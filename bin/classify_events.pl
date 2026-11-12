#!/usr/bin/env perl
# classify_events.pl - elements + polya + tsd + class profiles -> events.tsv (+ discarded.tsv)
# Every insertion of insertions.tsv ends up in exactly one of the two tables.
use strict;
use warnings;
use FindBin qw($RealBin);
use lib "$RealBin/../lib";
use Getopt::Long qw(GetOptions);
use Strongnose::Common qw(open_out read_tsv tsv_row need_columns parse_params read_bed4 intervals_len);

sub usage {
    return <<'USAGE';
usage: classify_events.pl --profiles class_profiles.tsv --elements elements.tsv --insertions insertions.tsv
                          --polya polya.tsv --tsd tsd.tsv [--structure ltr_structure.tsv] [options] > events.tsv

required:
  --profiles FILE               class_profiles.tsv
  --elements FILE               from merge_hits.pl
  --insertions FILE             from merge_hits.pl (--insertions-out)
  --polya FILE                  from find_polya.pl
  --tsd FILE                    from find_tsd.pl

options (default):
  --structure FILE              from check_ltr_structure.pl; needed to classify profiles with rule ltr_int_ltr or
                                solo_ltr (HERVK). Without it they are written as low with note rule_pending
  --out FILE                    events table (STDOUT)
  --discarded-out FILE          insertions that are not events, with the reason (none written)
  --bed FILE                    BED4 of all insertions; those without any RepeatMasker hit are listed as discarded
                                (reason no_repeat_hit)
  --artifact-max-tandem-frac F  discard when Simple_repeat + Low_complexity cover more than this fraction (0.50)
  --artifact-max-satellite-frac F  same for Satellite (0.30)
  --artifact-max-unknown-frac F same for Unknown (0.50)
  --polya-max-after-tail N      a tail only counts when at most N bp follow it in the insertion; the TSD copy
                                and a short transduction fit in this (50)
  --min-confidence low|probable|high   lowest level kept as an event, lower ones are discarded (low)
  --help                        this message

Per insertion:
  1. discarded when it has no candidate element, is mostly tandem repeat / satellite / unknown sequence,
     or no element covers enough bp to be primary
  2. the primary element is judged by each of its enabled profiles; the best result is kept
  3. gates of the profile: insertion length, element coverage, 3' end present (max_3end_left),
     strand switches (max_strand_switch). A failed gate caps the confidence at low
  4. evidence of the profile (the evidence column): polya = tail found, 3' end present and at most
     --polya-max-after-tail bp after it; tsd = TSD found, inside tsd_min..tsd_max, identity >= tsd_min_identity,
     and not a long duplication (tsd_capped)
  5. confidence: high = young + gates ok + all evidence; probable = young + gates ok + some evidence; low = the rest
  Profiles with rule ltr_int_ltr or solo_ltr use the table of check_ltr_structure.pl: young = the LTRs are young,
  the structure evidence is structure_ok, coverage is the fraction of the insertion spanned by the LTRs, and the
  poly-A and 3' end checks do not apply. Rules without an implementation give low and note rule_pending.
  Profiles with join_group_elements=yes (SVA) count elements of the same group and strand, and adjacent pieces of the
  classes of absorb_inside_span (within join_gap bp, 100), towards the coverage gate. Adjacent pieces count only
  inside the group span or on its 5' side; a piece after the 3' end is not part of the element. When the primary
  (largest) piece is not young or lacks its 3' end and a young piece of the same group and strand has it, that piece
  is the one judged and reported (note three_prime_piece).

events.tsv columns: query, chrom, pos, ins_len, profile, elem_id, repeat, tier, rm_family, strand, cov,
  div_corr_w, young (1/0), polya_found, polya_len, polya_purity, tsd_found, tsd_len, tsd_identity, tsd_pvalue,
  tsd_capped, tsd_seq, evidence_ok, evidence_total, gates_failed, confidence, other_elements, note,
  flags, tandem_frac, satellite_frac
  flags: element_after_3end (another element of 100 bp or more lies after the 3' end of the primary: 3' transduction,
  composite insertion or a second event; review by hand), fragment_no_3end (the 3' end gate failed), tandem_dup_of_flank (the TSD run reached the scan limit: the
  insertion repeats the neighbouring reference, typical of a length change of an existing repeat); '.' = none
discarded.tsv columns: query, chrom, pos, ins_len, reason, detail
USAGE
}

my %o = (
    'artifact-max-tandem-frac'    => 0.50,
    'artifact-max-satellite-frac' => 0.30,
    'artifact-max-unknown-frac'   => 0.50,
    'polya-max-after-tail'        => 50,
    'min-confidence'              => 'low',
);
GetOptions(\%o, qw(profiles=s elements=s insertions=s polya=s tsd=s structure=s bed=s out=s discarded-out=s
                   artifact-max-tandem-frac=f artifact-max-satellite-frac=f artifact-max-unknown-frac=f
                   polya-max-after-tail=i min-confidence=s help)) or die usage();
if ($o{help}) { print usage(); exit 0 }
die usage() unless defined $o{profiles} && defined $o{elements} && defined $o{insertions} && defined $o{polya} && defined $o{tsd};
my %RANK = (low => 0, probable => 1, high => 2);
die "--min-confidence must be low, probable or high\n" unless exists $RANK{ $o{'min-confidence'} };

# ---- profiles
my ($ph, $prows) = read_tsv($o{profiles});
need_columns($o{profiles}, $ph, qw(profile enabled rule evidence tsd_min tsd_max tsd_min_identity min_cov len_min len_max rule_params));
my %prof;
for my $r (@$prows) {
    next if $r->{enabled} eq 'no';
    $prof{ $r->{profile} } = { %$r, rp => parse_params($r->{rule_params}), ev => [ grep { length } split /;/, $r->{evidence} ] };
}

# ---- tables
my ($eh, $erows) = read_tsv($o{elements});
need_columns($o{elements}, $eh, qw(query elem_id rank is_primary repeat tier rm_family profiles strand cov cov_bases
                                   div_corr_w c_left young_profiles));
my %el_by_q;
push @{ $el_by_q{ $_->{query} } }, $_ for @$erows;

my ($ih, $irows) = read_tsv($o{insertions});
need_columns($o{insertions}, $ih, qw(query chrom bed_end q_len n_elements tandem_frac satellite_frac unknown_frac));

my ($ph2, $parows) = read_tsv($o{polya});
need_columns($o{polya}, $ph2, qw(elem_id polya_found polya_len polya_purity bases_after_tail));
my %polya = map { $_->{elem_id} => $_ } @$parows;

my ($th, $trows) = read_tsv($o{tsd});
need_columns($o{tsd}, $th, qw(query status tsd_found tsd_len tsd_identity tsd_pvalue tsd_capped tsd_seq));
my %tsd = map { $_->{query} => $_ } @$trows;

my (%struct, %struct_any);
if (defined $o{structure}) {
    my ($sh, $srows) = read_tsv($o{structure});
    need_columns($o{structure}, $sh, qw(query profile status ltr1_id ltr_identity span_frac young structure_ok fail_reason));
    for my $r (@$srows) {
        $struct_any{ $r->{query} }{ $r->{profile} } = $r;
        $struct{ $r->{query} }{ $r->{profile} } = $r if $r->{status} eq 'ok';
    }
}
my %el_by_id = map { $_->{elem_id} => $_ } @$erows;

# ---- outputs
my $out = open_out($o{out});
print {$out} tsv_row(qw(query chrom pos ins_len profile elem_id repeat tier rm_family strand cov div_corr_w young
                        polya_found polya_len polya_purity tsd_found tsd_len tsd_identity tsd_pvalue tsd_capped tsd_seq
                        evidence_ok evidence_total gates_failed confidence other_elements note flags tandem_frac satellite_frac));
my $dis;
if (defined $o{'discarded-out'}) {
    $dis = open_out($o{'discarded-out'});
    print {$dis} tsv_row(qw(query chrom pos ins_len reason detail));
}

my (%reason_n, %conf_n, $n_missing_tsd);
sub discard {
    my ($ins, $reason, $detail) = @_;
    $reason_n{$reason}++;
    print {$dis} tsv_row(@{$ins}{qw(query chrom)}, $ins->{bed_end}, $ins->{q_len}, $reason, defined $detail ? $detail : '.') if $dis;
}

sub list_has { my ($csv, $x) = @_; return scalar grep { $_ eq $x } split /,/, $csv }

sub evaluate {
    my ($pid, $e, $ins, $els) = @_;
    my $p  = $prof{$pid};
    my $rp = $p->{rp};
    my $rep_changed = 0;
    if (($rp->{join_group_elements} || '') eq 'yes' && defined $rp->{max_3end_left}
        && (!list_has($e->{young_profiles}, $pid) || $e->{c_left} > $rp->{max_3end_left})) {
        # the largest piece is not the one that carries the 3' end (SVA_C piece in front of a young SVA_F piece):
        # judge the element by the young piece of the same group and strand whose 3' end is present
        my ($alt) = sort { $b->{cov_bases} <=> $a->{cov_bases} }
                    grep { $_ != $e && $_->{group} eq $p->{group} && $_->{strand} eq $e->{strand}
                           && list_has($_->{young_profiles}, $pid) && $_->{c_left} <= $rp->{max_3end_left} } @$els;
        if ($alt) { $e = $alt; $rep_changed = 1 }
    }
    my @gates;
    push @gates, 'length'   if $ins->{q_len} < $p->{len_min} || $ins->{q_len} > $p->{len_max};
    my $cov = $e->{cov};
    if (($rp->{join_group_elements} || '') eq 'yes') {
        # one element can be annotated as pieces of several subfamilies (SVA_F + SVA_B + AluSx3...): count the pieces
        my @same = grep { $_->{group} eq $p->{group} && $_->{strand} eq $e->{strand} } @$els;
        my $lo = (sort { $a <=> $b } map { $_->{q_start} } @same)[0];
        my $hi = (sort { $b <=> $a } map { $_->{q_end} } @same)[0];
        my $gap = defined $rp->{join_gap} ? $rp->{join_gap} : 100;
        my @cls = split /,/, ($rp->{absorb_inside_span} || '');
        my %is_same = map { $_ => 1 } @same;
        # the Alu-like and repeat pieces of an SVA are at its 5' end: a piece entirely on the 3' side of the group span
        # (after the poly-A) is another element, not a component
        my @extra = grep {
            my $x = $_;
            my $past_3 = $e->{strand} eq '+' ? $x->{q_start} > $hi : $x->{q_end} < $lo;
            !$is_same{$x} && $x->{strand} eq $e->{strand} && !$past_3 && $x->{q_end} >= $lo - $gap && $x->{q_start} <= $hi + $gap
                && grep { index($x->{rm_family}, $_) == 0 } @cls
        } @$els;
        my $u = intervals_len(map { [$_->{q_start}, $_->{q_end}] } @same, @extra) / $ins->{q_len};
        $cov = $u if $u > $cov;
        $cov = 1 if $cov > 1;
    }
    $cov = sprintf('%.3f', $cov);
    push @gates, 'coverage' if $cov < $p->{min_cov};
    push @gates, 'three_prime_end_missing' if defined $rp->{max_3end_left} && $e->{c_left} > $rp->{max_3end_left};
    if (defined $rp->{max_strand_switch}) {
        my $sw = grep { $_ != $e && list_has($_->{profiles}, $pid) && $_->{strand} ne $e->{strand} } @$els;
        push @gates, 'strand_switches' if $sw > $rp->{max_strand_switch};
    }
    my $three_ok = !grep { $_ eq 'three_prime_end_missing' } @gates;
    my $young    = list_has($e->{young_profiles}, $pid) ? 1 : 0;

    my $pa = $polya{ $e->{elem_id} };
    my $td = $tsd{ $e->{query} };
    my $polya_ok = ($pa && $pa->{polya_found} eq '1' && $three_ok
                    && $pa->{bases_after_tail} <= $o{'polya-max-after-tail'}) ? 1 : 0;
    my $tsd_ok = ($td && $td->{status} eq 'ok' && $td->{tsd_found} eq '1' && $td->{tsd_capped} ne '1'
                  && $td->{tsd_len} >= $p->{tsd_min} && $td->{tsd_len} <= $p->{tsd_max}
                  && $td->{tsd_identity} >= $p->{tsd_min_identity}) ? 1 : 0;

    my ($ok, $total, @pending) = (0, 0);
    for my $ev (@{ $p->{ev} }) {
        $total++;
        if    ($ev eq 'polya') { $ok += $polya_ok }
        elsif ($ev eq 'tsd')   { $ok += $tsd_ok }
        else                   { push @pending, $ev }
    }
    my $conf = ($young && !@gates && $total && $ok == $total) ? 'high'
             : ($young && !@gates && $ok >= 1)                ? 'probable'
             : 'low';
    my $note = '.';
    if ($p->{rule} ne 'single_element') {
        $conf = 'low';
        my $sa = $struct_any{ $ins->{query} }{$pid};
        $note = (is_structure_rule($p->{rule}) && defined $o{structure})
              ? 'no_ltr_structure:' . ($sa ? $sa->{fail_reason} : 'no_ltr_pair')
              : "rule_pending:$p->{rule}";
    }
    elsif (@pending) {
        $note = 'unchecked_evidence:' . join(',', @pending);
    }
    elsif ($rep_changed) {
        $note = 'three_prime_piece';
    }
    return { pid => $pid, young => $young, polya => $pa, tsd => $td, polya_ok => $polya_ok, tsd_ok => $tsd_ok,
             ok => $ok, total => $total, gates => \@gates, conf => $conf, note => $note, cov => $cov, rep => $e };
}

sub is_structure_rule { return $_[0] eq 'ltr_int_ltr' || $_[0] eq 'solo_ltr' }

sub evaluate_ltr {
    my ($pid, $st, $ins) = @_;
    my $p  = $prof{$pid};
    my $e  = $el_by_id{ $st->{ltr1_id} };
    my @gates;
    push @gates, 'length'   if $ins->{q_len} < $p->{len_min} || $ins->{q_len} > $p->{len_max};
    push @gates, 'coverage' if $st->{span_frac} eq '.' || $st->{span_frac} < $p->{min_cov};
    my $young = $st->{young} eq '1' ? 1 : 0;
    my $td = $tsd{ $ins->{query} };
    my $tsd_ok = ($td && $td->{status} eq 'ok' && $td->{tsd_found} eq '1' && $td->{tsd_capped} ne '1'
                  && $td->{tsd_len} >= $p->{tsd_min} && $td->{tsd_len} <= $p->{tsd_max}
                  && $td->{tsd_identity} >= $p->{tsd_min_identity}) ? 1 : 0;
    my ($ok, $total, @pending) = (0, 0);
    for my $ev (@{ $p->{ev} }) {
        $total++;
        if    ($ev eq 'structure') { $ok += $st->{structure_ok} eq '1' ? 1 : 0 }
        elsif ($ev eq 'tsd')       { $ok += $tsd_ok }
        elsif ($ev eq 'polya')     { }
        else                       { push @pending, $ev }
    }
    my $conf = ($young && !@gates && $total && $ok == $total) ? 'high'
             : ($young && !@gates && $ok >= 1)                ? 'probable'
             : 'low';
    my $note = "ltr_identity=$st->{ltr_identity}";
    $note .= ';structure_fail=' . $st->{fail_reason} if $st->{fail_reason} ne '.';
    $note .= ';unchecked_evidence:' . join(',', @pending) if @pending;
    return { pid => $pid, young => $young, polya => undef, tsd => $td, polya_ok => 0, tsd_ok => $tsd_ok,
             ok => $ok, total => $total, gates => \@gates, conf => $conf, note => $note,
             elem => $e, cov => $st->{span_frac} };
}

# ---- main
my %seen;
for my $ins (@$irows) {
    $seen{ $ins->{query} } = 1;
    my @els = @{ $el_by_q{ $ins->{query} } || [] };
    my $fr  = sprintf('tandem=%s;satellite=%s;unknown=%s', @{$ins}{qw(tandem_frac satellite_frac unknown_frac)});

    if (!@els) { discard($ins, 'no_candidate_element', $fr); next }
    if ($ins->{tandem_frac} > $o{'artifact-max-tandem-frac'})       { discard($ins, 'tandem_repeat', $fr);     next }
    if ($ins->{satellite_frac} > $o{'artifact-max-satellite-frac'}) { discard($ins, 'satellite', $fr);         next }
    if ($ins->{unknown_frac} > $o{'artifact-max-unknown-frac'})     { discard($ins, 'unknown_sequence', $fr);  next }

    my ($prim) = grep { $_->{is_primary} eq '1' } @els;
    if (!$prim) {
        my ($top) = sort { $b->{cov_bases} <=> $a->{cov_bases} } @els;
        discard($ins, 'no_primary_element', "largest element $top->{repeat} covers $top->{cov_bases} bp");
        next;
    }
    my @pids = grep { $prof{$_} && !is_structure_rule($prof{$_}{rule}) } split /,/, $prim->{profiles};
    my %ltr_pids = %{ $struct{ $ins->{query} } || {} };
    for my $pid (keys %ltr_pids) { delete $ltr_pids{$pid} unless $prof{$pid} && $el_by_id{ $ltr_pids{$pid}{ltr1_id} } }
    if (!@pids && !%ltr_pids) {
        # the pending rules (solo_ltr, or ltr_int_ltr without --structure) keep the old behaviour
        @pids = grep { $prof{$_} } split /,/, $prim->{profiles};
    }
    if (!@pids && !%ltr_pids) { discard($ins, 'no_enabled_profile', $prim->{profiles}); next }

    my @res = ((map { evaluate($_, $prim, $ins, \@els) } @pids), (map { evaluate_ltr($_, $ltr_pids{$_}, $ins) } sort keys %ltr_pids));
    my ($best) = sort { $RANK{ $b->{conf} } <=> $RANK{ $a->{conf} } || $b->{ok} <=> $a->{ok} } @res;

    if ($RANK{ $best->{conf} } < $RANK{ $o{'min-confidence'} }) {
        discard($ins, 'below_min_confidence', "$best->{pid}: $best->{conf}");
        next;
    }
    $n_missing_tsd++ unless $best->{tsd};

    my ($pa, $td) = @{$best}{qw(polya tsd)};
    my $pe = $best->{elem} || $best->{rep} || $prim;
    my $cov = defined $best->{cov} ? $best->{cov} : $pe->{cov};
    my $others = join(';', map { "$_->{repeat}/$_->{strand}/$_->{cov}" } grep { $_ != $pe } @els) || '.';
    my @flags;
    push @flags, 'fragment_no_3end'  if grep { $_ eq 'three_prime_end_missing' } @{ $best->{gates} };
    push @flags, 'tandem_dup_of_flank' if $td && $td->{tsd_capped} eq '1';
    if (!$best->{elem}) {
        # another element of at least 100 bp after the 3' end of the primary: transduction, composite or a second event
        my $grp = $prof{ $best->{pid} }{group};
        my $after = grep { $_ != $pe && $_->{span} >= 100 && !($_->{group} eq $grp && $_->{strand} eq $pe->{strand})
                           && ($pe->{strand} eq '+' ? $_->{q_start} > $pe->{q_end} : $_->{q_end} < $pe->{q_start}) } @els;
        push @flags, 'element_after_3end' if $after;
    }
    $conf_n{ $best->{pid} }{ $best->{conf} }++;
    print {$out} tsv_row(
        @{$ins}{qw(query chrom bed_end q_len)}, $best->{pid}, @{$pe}{qw(elem_id repeat tier rm_family strand)}, $cov, $pe->{div_corr_w},
        $best->{young},
        $pa ? @{$pa}{qw(polya_found polya_len polya_purity)} : ('.', '.', '.'),
        $td ? @{$td}{qw(tsd_found tsd_len tsd_identity tsd_pvalue tsd_capped tsd_seq)} : ('.') x 6,
        $best->{ok}, $best->{total}, (@{ $best->{gates} } ? join(',', @{ $best->{gates} }) : '.'),
        $best->{conf}, $others, $best->{note}, (@flags ? join(',', @flags) : '.'), @{$ins}{qw(tandem_frac satellite_frac)},
    );
}

# insertions that RepeatMasker did not annotate at all
if (defined $o{bed}) {
    my $bed = read_bed4($o{bed});
    for my $id (sort keys %$bed) {
        next if $seen{$id};
        discard({ query => $id, chrom => $bed->{$id}[0], bed_end => $bed->{$id}[2], q_len => '.' }, 'no_repeat_hit', '.');
    }
}
close $out if defined $o{out};
close $dis if $dis;

# ---- summary
my $n_events = 0;
for my $pid (sort keys %conf_n) {
    my %c = %{ $conf_n{$pid} };
    $n_events += $_ for values %c;
    warn sprintf("classify_events: %-15s high %d, probable %d, low %d\n", $pid, $c{high} || 0, $c{probable} || 0, $c{low} || 0);
}
my $n_dis = 0;
$n_dis += $_ for values %reason_n;
warn sprintf("classify_events: %d events, %d discarded (%s)\n", $n_events, $n_dis,
             join(', ', map { "$_ $reason_n{$_}" } sort { $reason_n{$b} <=> $reason_n{$a} } keys %reason_n));
warn sprintf("classify_events: %d events have no row in %s (evidence counted as absent)\n", $n_missing_tsd, $o{tsd}) if $n_missing_tsd;
warn "classify_events: --discarded-out not given, so the discarded insertions are only counted\n" unless $dis;
