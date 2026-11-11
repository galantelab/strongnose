#!/usr/bin/env perl
# check_ltr_structure.pl - LTR - internal region - LTR structure of provirus profiles -> ltr_structure.tsv
# Used by every enabled profile whose rule is ltr_int_ltr (HERVK_provirus) or solo_ltr (HERVK_solo_ltr). Reads the parameters of the
# profile (ltr_names, internal_names, ltr_len, internal_min, same_orientation, ltr_min_identity) from rule_params.
use strict;
use warnings;
use FindBin qw($RealBin);
use lib "$RealBin/../lib";
use Getopt::Long qw(GetOptions);
use Strongnose::Common qw(open_out read_tsv tsv_row need_columns parse_params read_fasta revcomp intervals_len);

sub usage {
    return <<'USAGE';
usage: check_ltr_structure.pl --profiles class_profiles.tsv --elements elements.tsv --fasta insertions.fa [options] > ltr_structure.tsv

required:
  --profiles FILE     class_profiles.tsv (profiles with rule ltr_int_ltr and enabled=yes are checked)
  --elements FILE     from merge_hits.pl
  --fasta FILE        the insertion sequences RepeatMasker was run on (the coordinates of elements.tsv refer to them)

options (default):
  --out FILE          output table (STDOUT)
  --band N            half width, in bp, of the alignment band around the diagonal (80)
  --help              this message

Per insertion and provirus profile:
  Solo LTR profiles (rule solo_ltr): the longest LTR named in ltr_names; structure_ok = its length is inside the
  length range (ltr_len, or len_min-len_max of the profile), no internal-region element covers more than 20% of the
  insertion, and no second LTR is present. Written with status ok and the ltr2 columns '.'; status not_solo when it
  has an internal region or a second LTR (not a solo candidate, ignored by classify_events.pl).

  Provirus profiles (rule ltr_int_ltr):
  1. LTR candidates: elements named in ltr_names whose length is inside ltr_len
  2. pairs: LTR 1 before LTR 2, no overlap, same strand when same_orientation=yes, and the stretch between
     them at least internal_min bp (the internal region)
  3. the pair with the largest outer span is kept
  4. ltr_identity: the two LTRs, both in element orientation, are aligned (banded global alignment) and
     identity = matching columns / alignment columns
  5. internal region: internal_cov is the fraction of the stretch covered by elements named in
     internal_names (or with tier internal) of the same group; unannotated internal sequence is accepted,
     because the library may not name the internal region. Internal elements on the opposite strand of
     the LTRs over more than 20% of the stretch fail the orientation check
  6. structure_ok = pair found and stretch >= internal_min and orientation ok and ltr_identity >= ltr_min_identity
  young = 1 when both LTRs are young for the profile (young_profiles column of elements.tsv)

Columns: query, profile, status, ltr1_id, ltr2_id, ltr1_repeat, ltr2_repeat, strand, ltr1_start, ltr1_end, ltr2_start,
  ltr2_end, ltr1_len, ltr2_len, internal_len, internal_cov, internal_orient_ok, span_frac, ltr_identity, young,
  structure_ok, fail_reason
status: ok (pair found) | one_ltr (a single LTR: solo LTR) | no_pair (two LTRs, but none satisfies the pair conditions)
(insertions without any LTR of the profile are not listed)
fail_reason: comma list of internal_too_short, internal_orientation, ltr_identity, single_ltr, ltr_length, no_pair, or '.'
USAGE
}

my %o = (band => 80);
GetOptions(\%o, qw(profiles=s elements=s fasta=s out=s band=i help)) or die usage();
if ($o{help}) { print usage(); exit 0 }
die usage() unless defined $o{profiles} && defined $o{elements} && defined $o{fasta};

# ---- profiles
my ($ph, $prows) = read_tsv($o{profiles});
need_columns($o{profiles}, $ph, qw(profile group enabled rule min_cov len_min len_max rule_params));
my @prof;
for my $r (@$prows) {
    next unless $r->{enabled} ne 'no' && ($r->{rule} eq 'ltr_int_ltr' || $r->{rule} eq 'solo_ltr');
    my $rp = parse_params($r->{rule_params});
    my $solo = $r->{rule} eq 'solo_ltr';
    for my $k ($solo ? qw(ltr_names) : qw(ltr_names internal_names ltr_len internal_min ltr_min_identity)) {
        die "$o{profiles}: profile $r->{profile}: rule_params needs $k\n" unless defined $rp->{$k};
    }
    my $len = defined $rp->{ltr_len} ? $rp->{ltr_len} : "$r->{len_min}-$r->{len_max}";
    my ($lmin, $lmax) = $len =~ /^(\d+)-(\d+)$/ or die "$o{profiles}: profile $r->{profile}: ltr_len must be MIN-MAX\n";
    push @prof, {
        id => $r->{profile}, group => $r->{group}, min_cov => $r->{min_cov}, solo => $solo,
        ltr_names => { map { lc($_) => 1 } split /,/, $rp->{ltr_names} },
        int_names => { map { lc($_) => 1 } split /,/, (defined $rp->{internal_names} ? $rp->{internal_names} : '') },
        lmin => $lmin, lmax => $lmax, int_min => $rp->{internal_min},
        same => (($rp->{same_orientation} || 'yes') eq 'yes'), min_id => $rp->{ltr_min_identity},
    };
}
die "no enabled profile with rule ltr_int_ltr in $o{profiles}\n" unless @prof;

# ---- elements
my ($eh, $erows) = read_tsv($o{elements});
need_columns($o{elements}, $eh, qw(query elem_id repeat group tier profiles strand q_start q_end span young_profiles));
my %el_by_q;
push @{ $el_by_q{ $_->{query} } }, $_ for @$erows;
my $seqs = read_fasta($o{fasta});

# ---- banded global alignment; returns (matching columns, alignment columns)
sub align_identity {
    my ($sa, $sb, $band) = @_;
    my @a = split //, $sa;
    my @b = split //, $sb;
    my ($n, $m) = (scalar @a, scalar @b);
    return (0, 1) unless $n && $m;
    my ($MATCH, $MIS, $GAP, $NEG) = (1, -1, -2, -1e9);
    my $dmin = ($m - $n < 0 ? $m - $n : 0) - $band;
    my $dmax = ($m - $n > 0 ? $m - $n : 0) + $band;
    my $w = $dmax - $dmin + 1;
    my (@S, @P);
    for my $i (0 .. $n) {
        $S[$i] = [ ($NEG) x $w ];
        $P[$i] = [ (0) x $w ];
        my $jlo = $i + $dmin < 0 ? 0 : $i + $dmin;
        my $jhi = $i + $dmax > $m ? $m : $i + $dmax;
        for my $j ($jlo .. $jhi) {
            my $k = $j - $i - $dmin;
            if ($i == 0 && $j == 0) { $S[0][$k] = 0; next }
            my ($best, $ptr) = ($NEG, 0);
            if ($i > 0 && $j > 0) {                                   # diagonal: same k
                my $s = $S[$i - 1][$k];
                if ($s > $NEG) {
                    $s += ($a[$i - 1] eq $b[$j - 1] && $a[$i - 1] ne 'N') ? $MATCH : $MIS;
                    ($best, $ptr) = ($s, 1) if $s > $best;
                }
            }
            if ($i > 0 && $k + 1 < $w) {                              # up: (i-1, j) has k+1
                my $s = $S[$i - 1][$k + 1];
                if ($s > $NEG) { $s += $GAP; ($best, $ptr) = ($s, 2) if $s > $best }
            }
            if ($j > 0 && $k - 1 >= 0) {                              # left: (i, j-1) has k-1
                my $s = $S[$i][$k - 1];
                if ($s > $NEG) { $s += $GAP; ($best, $ptr) = ($s, 3) if $s > $best }
            }
            $S[$i][$k] = $best;
            $P[$i][$k] = $ptr;
        }
    }
    my ($i, $j, $match, $cols) = ($n, $m, 0, 0);
    while ($i > 0 || $j > 0) {
        my $k = $j - $i - $dmin;
        my $p = $P[$i][$k];
        if    ($p == 1) { $match++ if $a[$i - 1] eq $b[$j - 1] && $a[$i - 1] ne 'N'; $i--; $j-- }
        elsif ($p == 2) { $i-- }
        elsif ($p == 3) { $j-- }
        else            { last }
        $cols++;
    }
    return ($match, $cols || 1);
}

my $out = open_out($o{out});
print {$out} tsv_row(qw(query profile status ltr1_id ltr2_id ltr1_repeat ltr2_repeat strand ltr1_start ltr1_end ltr2_start
                        ltr2_end ltr1_len ltr2_len internal_len internal_cov internal_orient_ok span_frac ltr_identity
                        young structure_ok fail_reason));

my (%status, %n_ok);
for my $q (sort keys %el_by_q) {
    my @els = @{ $el_by_q{$q} };
    my $qseq = $seqs->{$q};
    for my $p (@prof) {
        my @mine = grep { grep({ $_ eq $p->{id} } split /,/, $_->{profiles}) || $_->{group} eq $p->{group} } @els;
        my @ltr = grep { $p->{ltr_names}{ lc $_->{repeat} } } @mine;
        next unless @ltr;                      # no LTR of this profile at all: nothing to report
        if ($p->{solo}) {
            my ($x) = sort { $b->{span} <=> $a->{span} } @ltr;
            my $q_len = defined $qseq ? length $qseq : 0;
            my (@iv, @fail);
            for my $e (@mine) {
                next unless $p->{int_names}{ lc $e->{repeat} } || $e->{tier} eq 'internal';
                push @iv, [$e->{q_start}, $e->{q_end}];
            }
            my $icov = $q_len ? intervals_len(@iv) / $q_len : 0;
            push @fail, 'ltr_length'       if $x->{span} < $p->{lmin} || $x->{span} > $p->{lmax};
            push @fail, 'internal_present' if $icov > 0.20;
            push @fail, 'multiple_ltr'     if @ltr > 1;
            my $young = (grep { $_ eq $p->{id} } split /,/, $x->{young_profiles}) ? 1 : 0;
            my $ok = @fail ? 0 : 1;
            # a provirus-shaped insertion (internal region or a second LTR) is not a solo LTR candidate at all
            my $stat = (grep { $_ eq 'internal_present' || $_ eq 'multiple_ltr' } @fail) ? 'not_solo' : 'ok';
            $status{$stat}++;
            $n_ok{ $p->{id} }++ if $ok;
            print {$out} tsv_row($q, $p->{id}, $stat, $x->{elem_id}, '.', $x->{repeat}, '.', $x->{strand}, $x->{q_start}, $x->{q_end},
                                 '.', '.', $x->{span}, '.', '.', sprintf('%.3f', $icov), 1,
                                 $q_len ? sprintf('%.3f', $x->{span} / $q_len) : '.', '.', $young, $ok,
                                 @fail ? join(',', @fail) : '.');
            next;
        }
        my $fail = sub {
            my ($status, $reason) = @_;
            $status{$status}++;
            print {$out} tsv_row($q, $p->{id}, $status, ('.') x 16, 0, 0, $reason);
        };

        my @cand = grep { $_->{span} >= $p->{lmin} && $_->{span} <= $p->{lmax} } @ltr;
        if (@ltr < 2) { $fail->('one_ltr', 'single_ltr'); next }
        if (@cand < 2) { $fail->('no_pair', 'ltr_length'); next }
        # pairs
        my ($best, $best_span);
        for my $x (@cand) {
            for my $y (@cand) {
                next if $x == $y || $x->{q_start} >= $y->{q_start} || $x->{q_end} >= $y->{q_start};
                next if $p->{same} && $x->{strand} ne $y->{strand};
                my $gap = $y->{q_start} - $x->{q_end} - 1;
                next if $gap < $p->{int_min};
                my $span = $y->{q_end} - $x->{q_start} + 1;
                if (!defined $best_span || $span > $best_span) { ($best, $best_span) = ([$x, $y, $gap], $span) }
            }
        }
        if (!$best) { $fail->('no_pair', 'no_pair'); next }

        my ($x, $y, $gap) = @$best;
        my $strand = $x->{strand} eq $y->{strand} ? $x->{strand} : '.';
        my $q_len = defined $qseq ? length $qseq : 0;

        # internal region
        my ($g1, $g2) = ($x->{q_end} + 1, $y->{q_start} - 1);
        my (@iv, @opp);
        for my $e (@mine) {
            next if $e == $x || $e == $y;
            next unless $p->{int_names}{ lc $e->{repeat} } || $e->{tier} eq 'internal';
            my ($s, $t) = ($e->{q_start} < $g1 ? $g1 : $e->{q_start}, $e->{q_end} > $g2 ? $g2 : $e->{q_end});
            next if $t < $s;
            push @iv, [$s, $t];
            push @opp, [$s, $t] if $strand ne '.' && $e->{strand} ne $strand;
        }
        my $icov = $gap ? intervals_len(@iv) / $gap : 0;
        my $oppf = $gap ? intervals_len(@opp) / $gap : 0;
        my $orient_ok = ($strand ne '.' && $oppf <= 0.20) ? 1 : 0;
        $orient_ok = 1 if $strand eq '.' && !$p->{same} && $oppf <= 0.20;

        # identity
        my $ident = '.';
        if (defined $qseq) {
            my $s1 = substr($qseq, $x->{q_start} - 1, $x->{span});
            my $s2 = substr($qseq, $y->{q_start} - 1, $y->{span});
            if ($y->{strand} eq '-') { $s2 = revcomp($s2) }
            if ($x->{strand} eq '-') { $s1 = revcomp($s1) }
            my ($m, $c) = align_identity($s1, $s2, $o{band});
            $ident = sprintf('%.3f', $m / $c);
        }
        my @fail;
        push @fail, 'internal_too_short' if $gap < $p->{int_min};
        push @fail, 'internal_orientation' unless $orient_ok;
        push @fail, 'ltr_identity' if $ident eq '.' || $ident < $p->{min_id};
        my $ok = @fail ? 0 : 1;
        my $young = (grep({ $_ eq $p->{id} } split /,/, $x->{young_profiles}) && grep({ $_ eq $p->{id} } split /,/, $y->{young_profiles})) ? 1 : 0;
        $n_ok{ $p->{id} }++ if $ok;
        $status{ok}++;
        print {$out} tsv_row($q, $p->{id}, 'ok', $x->{elem_id}, $y->{elem_id}, $x->{repeat}, $y->{repeat}, $strand,
                             @{$x}{qw(q_start q_end)}, @{$y}{qw(q_start q_end)}, $x->{span}, $y->{span}, $gap,
                             sprintf('%.3f', $icov), $orient_ok, $q_len ? sprintf('%.3f', $best_span / $q_len) : '.',
                             $ident, $young, $ok, @fail ? join(',', @fail) : '.');
    }
}
close $out if defined $o{out};
warn sprintf("check_ltr_structure: %s\n", join(', ', map { "$_ $status{$_}" } sort keys %status) || 'no candidates');
warn sprintf("check_ltr_structure: structure ok: %s\n", join(', ', map { "$_ $n_ok{$_}" } sort keys %n_ok) || 'none');
