#!/usr/bin/env perl
# select_review.pl - events.tsv (+ discarded.tsv) -> review.tsv: events for manual review, with what to check
use strict;
use warnings;
use FindBin qw($RealBin);
use lib "$RealBin/../lib";
use Getopt::Long qw(GetOptions);
use List::Util qw(shuffle);
use Strongnose::Common qw(open_out open_in read_tsv tsv_row need_columns read_fasta);

sub usage {
    return <<'USAGE';
usage: select_review.pl --events events.tsv [options] > review.tsv

required:
  --events FILE          from classify_events.pl

options (default):
  --discarded FILE       discarded.tsv; adds the false-negative checks (insertions that were dropped)
  --out FILE             review table (STDOUT)
  --per-category N       events sampled in a sampled category (10); "all" categories are capped at --max-all
  --max-all N            cap of the categories that list every event (40)
  --seed N               seed of the sampling; the same seed gives the same list (1)
  --window N             bp on each side of the position in the IGV locus (150)
  --extra FILE           TSV with a header whose first column is the insertion ID; its other columns are appended
                         (read support, allele fraction, QUAL...), for example:
                         bcftools query -H -f '%ID\t%INFO/SUPPORT\t%INFO/VAF\t%QUAL\n' calls.vcf.gz > extra.tsv
  --fasta FILE           insertion FASTA (the one RepeatMasker was run on); with --fasta-out the selected sequences are
                         written there, to paste into BLAT / Dfam
  --fasta-out FILE       FASTA of the selected insertions
  --checklist-out FILE   general checklist of what the pipeline cannot see (text)
  --help                 this message

Categories, in priority order (an event appears once, under the first category that picks it):
  validation of calls   high_alu, high_l1, high_sva, probable_alu, probable_l1, probable_sva, hervk_ok
  questions on rules    sva_element_after_3end, sva_no3end_with_evidence, l1_no3end_with_evidence, gate_rescue_alu,
                        gate_rescue_l1, tandem_dup_young
  calibration           age_near_ceiling (high Alus with the highest corrected divergence), old_alu_with_evidence
                        (AluS and older with poly-A and TSD), short_tsd_high (high with a TSD of 9 bp or less)
  false negatives       dropped insertions by reason (no_candidate_element, no_primary_element, no_repeat_hit,
                        tandem_repeat, satellite)

review.tsv columns: review_id, category, priority, query, chrom, pos, locus, ins_len, profile, repeat, tier, confidence,
  flags, cov, div_corr_w, polya_len, tsd_len, tsd_identity, gates_failed, other_elements, note, why_review, to_check,
  [extra columns], verdict, reviewer, comment
  verdict: for the reviewers (true_insertion | false_positive | uncertain); the three last columns are left empty
USAGE
}

my %o = ('per-category' => 10, 'max-all' => 40, seed => 1, window => 150);
GetOptions(\%o, qw(events=s discarded=s out=s per-category=i max-all=i seed=i window=i extra=s fasta=s fasta-out=s
                   checklist-out=s help)) or die usage();
if ($o{help}) { print usage(); exit 0 }
die usage() unless defined $o{events};
die "--fasta-out needs --fasta\n" if defined $o{'fasta-out'} && !defined $o{fasta};
srand($o{seed});

my ($eh, $ev) = read_tsv($o{events});
need_columns($o{events}, $eh, qw(query chrom pos ins_len profile repeat tier rm_family strand cov div_corr_w young polya_found
                                 polya_len tsd_found tsd_len tsd_identity tsd_capped gates_failed confidence other_elements
                                 note flags));
my $dis = [];
if (defined $o{discarded}) {
    my ($dh, $drows) = read_tsv($o{discarded});
    need_columns($o{discarded}, $dh, qw(query chrom pos ins_len reason detail));
    $dis = $drows;
}
my (@extra_cols, %extra);
if (defined $o{extra}) {
    my $fh = open_in($o{extra});
    my $hdr = <$fh>;
    chomp $hdr;
    $hdr =~ s/^#//;
    my ($k, @cols) = split /\t/, $hdr, -1;
    @extra_cols = map { my $c = $_; $c =~ s/^\[\d+\]//; $c =~ s/^(?:INFO|FORMAT)\///; $c } @cols;
    while (my $line = <$fh>) {
        chomp $line;
        my ($id, @v) = split /\t/, $line, -1;
        $extra{$id} = \@v;
    }
    close $fh;
}

# ---------------------------------------------------------------- categories
sub has_flag { my ($e, $f) = @_; return scalar grep { $_ eq $f } split /,/, $e->{flags} }
my $both = sub { $_[0]{polya_found} eq '1' && $_[0]{tsd_found} eq '1' };

my @cat = (
    ['high_alu',  'validation', 'Alu called high',
        sub { $_[0]{profile} eq 'Alu' && $_[0]{confidence} eq 'high' }, 'sample',
        'real insertion at the locus; TSD at the junctions; poly-A at the 3\' end; Alu orientation'],
    ['high_l1',   'validation', 'L1 called high',
        sub { $_[0]{profile} eq 'L1' && $_[0]{confidence} eq 'high' }, 'sample',
        'real insertion; 5\' truncation and 3\' end; TSD; poly-A; 5\' inversion'],
    ['high_sva',  'validation', 'SVA called high',
        sub { $_[0]{profile} eq 'SVA' && $_[0]{confidence} eq 'high' }, 'all',
        'real insertion; VNTR size; is it inside a reference SVA (length change of an existing element)?; TSD; poly-A'],
    ['probable_alu', 'validation', 'Alu called probable',
        sub { $_[0]{profile} eq 'Alu' && $_[0]{confidence} eq 'probable' }, 'sample',
        'which evidence is missing and is it really absent (poly-A or TSD)?'],
    ['probable_l1',  'validation', 'L1 called probable',
        sub { $_[0]{profile} eq 'L1' && $_[0]{confidence} eq 'probable' }, 'sample',
        'which evidence is missing and is it really absent (poly-A or TSD)?'],
    ['probable_sva', 'validation', 'SVA called probable',
        sub { $_[0]{profile} eq 'SVA' && $_[0]{confidence} eq 'probable' }, 'all',
        'which evidence is missing and is it really absent (poly-A or TSD)?'],
    ['hervk_ok', 'validation', 'HERVK with a valid structure or a solo LTR candidate',
        sub { $_[0]{profile} =~ /^HERVK/ && $_[0]{confidence} ne 'low' }, 'all',
        'LTR-internal-LTR layout; LTR identity; for a solo LTR: could it come from recombination of a reference provirus?'],
    ['sva_element_after_3end', 'rule question', 'SVA with another element after its 3\' end',
        sub { $_[0]{profile} eq 'SVA' && has_flag($_[0], 'element_after_3end') }, 'all',
        'sequence after the SVA 3\' end: 3\' transduction, composite insertion or a second event?; is the last tail the true end?'],
    ['sva_no3end_with_evidence', 'rule question', 'SVA whose 3\' end gate failed although it has poly-A and TSD',
        sub { $_[0]{profile} eq 'SVA' && has_flag($_[0], 'fragment_no_3end') && $both->($_[0]) }, 'all',
        'what covers the end of the insertion (no repeat hit)?; should the 3\' gate (max_3end_left) be waived with both evidence?'],
    ['l1_no3end_with_evidence', 'rule question', 'young L1 whose 3\' end gate failed although it has poly-A and TSD',
        sub { $_[0]{profile} eq 'L1' && $_[0]{young} eq '1' && has_flag($_[0], 'fragment_no_3end') && $both->($_[0]) }, 'sample',
        '3\' end of the L1 present or not; transduction?'],
    ['gate_rescue_alu', 'rule question', 'young Alu rejected by a gate although it has poly-A and TSD',
        sub { $_[0]{profile} eq 'Alu' && $_[0]{young} eq '1' && $_[0]{gates_failed} ne '.' && $both->($_[0]) }, 'sample',
        'does the failed gate reject a real insertion (length range, coverage, strand switch)?'],
    ['gate_rescue_l1', 'rule question', 'young L1 rejected by a gate although it has poly-A and TSD',
        sub { $_[0]{profile} eq 'L1' && $_[0]{young} eq '1' && $_[0]{gates_failed} ne '.' && $both->($_[0]) }, 'sample',
        'does the failed gate reject a real insertion?'],
    ['tandem_dup_young', 'rule question', 'young element whose TSD run reaches the scan limit',
        sub { $_[0]{young} eq '1' && has_flag($_[0], 'tandem_dup_of_flank') }, 'sample',
        'does the insertion repeat the neighbouring reference (length change of an existing repeat) or is it a real new copy next to a repeat?'],
    ['age_near_ceiling', 'calibration', 'high Alu with the highest corrected divergence',
        sub { $_[0]{profile} eq 'Alu' && $_[0]{confidence} eq 'high' }, 'top_div',
        'is the element really young? subfamily diagnostic positions (Dfam / BLAT); calibrates the divergence ceilings'],
    ['old_alu_with_evidence', 'calibration', 'AluS or older with poly-A and TSD (excluded by age)',
        sub { $_[0]{profile} eq 'Alu' && $_[0]{young} eq '0' && $_[0]{tier} =~ /^(intermediate|old)$/ && $both->($_[0]) }, 'sample',
        'is it a real polymorphic insertion of an older subfamily that the age filter wrongly excludes?'],
    ['short_tsd_high', 'calibration', 'called high with a TSD of 9 bp or less',
        sub { $_[0]{confidence} eq 'high' && $_[0]{tsd_found} eq '1' && $_[0]{tsd_len} ne '.' && $_[0]{tsd_len} <= 9 }, 'sample',
        'is the short TSD real (decoy: about 14% of 5-9 bp TSDs are chance)?'],
);

my %chosen;
my @rows;
my $prio = 0;
for my $c (@cat) {
    my ($id, $group, $desc, $sel, $mode, $check) = @$c;
    $prio++;
    my @cand = grep { !$chosen{ $_->{query} } && $sel->($_) } @$ev;
    my $n_total = scalar @cand;
    my @pick;
    if    ($mode eq 'all')     { @pick = splice(@{[ sort { $a->{chrom} cmp $b->{chrom} || $a->{pos} <=> $b->{pos} } @cand ]}, 0, $o{'max-all'}) }
    elsif ($mode eq 'top_div') { @pick = splice(@{[ sort { $b->{div_corr_w} <=> $a->{div_corr_w} } @cand ]}, 0, $o{'per-category'}) }
    else                       { @pick = splice(@{[ shuffle(@cand) ]}, 0, $o{'per-category'}) }
    for my $e (@pick) {
        $chosen{ $e->{query} } = 1;
        push @rows, { src => 'event', cat => $id, prio => $prio, e => $e, why => "$desc ($n_total in the table)", check => $check };
    }
}

# dropped insertions
my %dchecks = (
    no_candidate_element => 'does the insertion contain a mobile element that no profile recognised (unlisted subfamily, divergent copy)? BLAT/Dfam the sequence',
    no_primary_element   => 'is there a real element split into pieces too short to be primary?',
    no_repeat_hit        => 'is it a real insertion of a mobile element RepeatMasker missed (very divergent, or a short element)?',
    tandem_repeat        => 'tandem repeat expansion or a real element inside a repeat?',
    satellite            => 'satellite expansion (centromeric on T2T) or an element?',
);
for my $reason (qw(no_candidate_element no_primary_element no_repeat_hit tandem_repeat satellite)) {
    $prio++;
    my @cand = grep { $_->{reason} eq $reason } @$dis;
    next unless @cand;
    my $mode_all = $reason eq 'tandem_repeat' || $reason eq 'satellite';
    my @pick = $mode_all ? splice(@{[ sort { $a->{chrom} cmp $b->{chrom} || $a->{pos} <=> $b->{pos} } @cand ]}, 0, $o{'max-all'})
                         : splice(@{[ shuffle(@cand) ]}, 0, $o{'per-category'});
    for my $d (@pick) {
        my $e = { query => $d->{query}, chrom => $d->{chrom}, pos => $d->{pos}, ins_len => $d->{ins_len}, profile => '.',
                  repeat => '.', tier => '.', confidence => 'dropped', flags => '.', cov => '.', div_corr_w => '.',
                  polya_len => '.', tsd_len => '.', tsd_identity => '.', gates_failed => '.', other_elements => '.',
                  note => "$d->{reason}:$d->{detail}" };
        push @rows, { src => 'dropped', cat => "dropped_$reason", prio => $prio, e => $e,
                      why => "dropped insertion ($reason, " . scalar(@cand) . ' in the table)', check => $dchecks{$reason} };
    }
}

# ---------------------------------------------------------------- hints of one event
sub hints {
    my ($r) = @_;
    my $e = $r->{e};
    my @h;
    if ($r->{src} eq 'event') {
        push @h, 'TSD not detected: look for duplicated sequence at both junctions' if $e->{tsd_found} ne '1';
        push @h, 'TSD run reaches the scan limit: insertion repeats the neighbouring reference' if $e->{tsd_capped} eq '1';
        push @h, "short TSD ($e->{tsd_len} bp): could be chance" if $e->{tsd_found} eq '1' && $e->{tsd_len} ne '.' && $e->{tsd_len} <= 9 && $e->{tsd_capped} ne '1';
        push @h, 'no poly-A found: inspect the 3\' end of the insertion' if $e->{polya_found} ne '1';
        for my $g (split /,/, $e->{gates_failed}) {
            push @h, { length => 'insertion length outside the profile range',
                       coverage => 'element covers little of the insertion: what is the rest?',
                       three_prime_end_missing => 'element stops before the 3\' end of the consensus: what is after it?',
                       strand_switches => 'pieces on both strands: inversion or chimeric insertion?' }->{$g} if $g ne '.' && $g =~ /^(length|coverage|three_prime_end_missing|strand_switches)$/;
        }
        push @h, 'another element lies after the 3\' end: transduction, composite or a second event' if has_flag($e, 'element_after_3end');
        push @h, "note: $e->{note}" if $e->{note} ne '.';
        push @h, "other elements: $e->{other_elements}" if $e->{other_elements} ne '.';
    }
    push @h, 'IGV: insertion inside a reference repeat of the same family (length change)?';
    push @h, 'reference context: segmental duplication, gap, centromere or low mappability (hg38)?';
    return @h;
}

# ---------------------------------------------------------------- output
my $out = open_out($o{out});
print {$out} tsv_row(qw(review_id category priority query chrom pos locus ins_len profile repeat tier confidence flags cov
                        div_corr_w polya_len tsd_len tsd_identity gates_failed other_elements note why_review to_check),
                     @extra_cols, qw(verdict reviewer comment));
my $n = 0;
@rows = sort { $a->{prio} <=> $b->{prio} || $a->{e}{chrom} cmp $b->{e}{chrom} || $a->{e}{pos} <=> $b->{e}{pos} } @rows;
my %want;
for my $r (@rows) {
    my $e = $r->{e};
    $n++;
    $want{ $e->{query} } = 1;
    my $lo = $e->{pos} - $o{window};
    $lo = 1 if $lo < 1;
    my $x = $extra{ $e->{query} } || [];
    print {$out} tsv_row(sprintf('R%03d', $n), $r->{cat}, $r->{prio}, @{$e}{qw(query chrom pos)},
                         "$e->{chrom}:$lo-" . ($e->{pos} + $o{window}), @{$e}{qw(ins_len profile repeat tier confidence flags cov
                         div_corr_w polya_len tsd_len tsd_identity gates_failed other_elements note)},
                         $r->{why}, join('; ', $r->{check}, hints($r)),
                         (map { defined $x->[$_] ? $x->[$_] : '.' } 0 .. $#extra_cols), '', '', '');
}
close $out if defined $o{out};

if (defined $o{'fasta-out'}) {
    my $seqs = read_fasta($o{fasta});
    my $fo = open_out($o{'fasta-out'});
    my $miss = 0;
    for my $r (@rows) {
        my $q = $r->{e}{query};
        if (!defined $seqs->{$q}) { $miss++; next }
        print {$fo} ">$q $r->{cat}\n", join("\n", unpack('(A80)*', $seqs->{$q})), "\n";
    }
    close $fo;
    warn "select_review: $miss selected insertions are not in $o{fasta}\n" if $miss;
}

if (defined $o{'checklist-out'}) {
    my $cl = open_out($o{'checklist-out'});
    print {$cl} <<'TXT';
Checklist for the manual review (what the pipeline cannot see)

Per insertion, in IGV (locus column of review.tsv):
 1. Read support: how many reads show the insertion, are they long and consistent, what is the allele fraction
    (extra columns, if given). Low support and a clipped or soft-clipped pattern suggest a false call.
 2. Insertion sequence: paste it into BLAT or Dfam (review.fa). Does the subfamily call agree? Check diagnostic
    positions of the subfamily, not only the RepeatMasker name.
 3. TSD: do the reference bases on the two sides of the insertion point repeat the same short sequence, in the
    right length for the family (Alu 7-20, L1 7-30, SVA 7-20)?
 4. Poly-A: is there an A-rich tail at the 3' end of the element (T-rich at the 5' end on the minus strand)?
 5. Reference context: is the insertion point inside a reference repeat of the same family (a length change of an
    existing element, VNTR or poly-A length) and not a new insertion? Is the region a segmental duplication, a gap, a
    centromere model or low mappability (hg38)? On T2T the region may be resolved.
 6. Known polymorphism: is the insertion in a database of polymorphic mobile elements (germline), or only in this sample?
 7. Orientation: do the element strand and the poly-A side agree?

Verdict (column verdict): true_insertion | false_positive | uncertain; add the reason in comment. The verdicts give
precision per confidence level and per category, and show which rule to change.
TXT
    close $cl;
}

my %by;
$by{ $_->{cat} }++ for @rows;
warn sprintf("select_review: %d events selected: %s\n", $n, join(', ', map { "$_ $by{$_}" } sort { $by{$b} <=> $by{$a} || $a cmp $b } keys %by));
