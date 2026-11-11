#!/usr/bin/env perl
# find_tsd.pl - insertions.tsv + insertion FASTA + flanks.tsv -> tsd.tsv
#
# A target site duplication makes  reference-left + insertion + reference-right  periodic with a period equal to the
# insertion length: the two copies of the TSD are exactly that far apart, however the aligner split the TSD between the
# insertion and the flanks. So the TSD is a run of positions where S[i] = S[i + n] (n = insertion length) that touches
# the junction. No element boundary is used, which makes the test independent of truncations and transductions.
# The run is chosen by significance (chance of such a run in unrelated sequence), not by score, so a few ONT errors
# inside the copy do not split it.
use strict;
use warnings;
use FindBin qw($RealBin);
use lib "$RealBin/../lib";
use Getopt::Long qw(GetOptions);
use Strongnose::Common qw(open_out read_tsv tsv_row need_columns read_fasta);

sub usage {
    return <<'USAGE';
usage: find_tsd.pl --fasta insertions.fa --flanks flanks.tsv [options] insertions.tsv > tsd.tsv

required:
  --fasta FILE                  insertion sequences (>Sniffles ID ...)
  --flanks FILE                 output of extract_flanks.pl
  insertions.tsv                from merge_hits.pl; needs the columns query and n_elements

options (default):
  --out FILE                    output table (STDOUT)
  --min-elements N              only insertions with at least N elements are analysed (1; 0 = all)
  --tsd-min-len N               shortest run called a TSD (4); each profile applies its own range later
  --tsd-min-identity F          minimum identity of the two copies, 0-1 (0.85); ONT errors in the copy are tolerated
  --tsd-max-pvalue F            maximum adjusted p-value for tsd_found (0.01; decoy false-positive rate 1.6% with window 2)
  --tsd-pos-window N            the junction is also tried N bp to each side, since the Sniffles position is imprecise (2; real TSDs sit at shift 0 or +-1)
  --tsd-scan-max N              longest run scanned on each side of the junction (40); a run that reaches it is capped
  --tsd-null-match F            chance that two unrelated bases match, for the p-value (0.30)
  --anchor-base auto|yes|no     does the insertion FASTA start with the reference base at the position (a VCF ALT anchor)?
                                auto compares the first base with the reference over all insertions (auto)
  --decoy                       control run: every insertion is compared with the flank of the NEXT insertion (unrelated
                                sequence), so the fraction found estimates the false-positive rate of the thresholds
  --help                        this message

Columns of tsd.tsv, one row per analysed insertion:
  query, chrom, pos, ins_len
  status           ok | no_flank | no_sequence | too_short
  tsd_found        1 / 0 ('.' unless status is ok); the best run is reported even when 0
  tsd_len, tsd_matches, tsd_identity
  tsd_pvalue       chance of a run this good at the junction in unrelated sequence: binomial tail x (len+1) x (2*pos-window+1)
  tsd_delta        junction shift (bp) that gave the best run; 0 means the Sniffles position was exact
  tsd_start_offset first base of the first copy relative to the first insertion base (negative: it lies in the reference flank)
  tsd_capped       1 when a found run reaches the scan limit: a long duplication of the neighbouring reference (tandem
                   duplication), not a TSD. The limit is --tsd-scan-max, or half the insertion length for short insertions
  tsd_seq, tsd_seq2   the two copies, first 60 bases; differences between them are ONT errors or the edge of the call

Among runs that are about equally significant (within a factor 2) the one with the highest identity, then the shortest,
is reported. tsd_len can still be 1-2 bp too long when the bases next to the TSD match by chance, notably next to an A-rich tail.
USAGE
}

my %o = (
    'min-elements'     => 1,
    'tsd-min-len'      => 4,
    'tsd-min-identity' => 0.85,
    'tsd-max-pvalue'   => 0.01,
    'tsd-pos-window'   => 2,
    'tsd-scan-max'     => 40,
    'tsd-null-match'   => 0.30,
    'anchor-base'      => 'auto',
);
GetOptions(\%o, qw(fasta=s flanks=s out=s min-elements=i tsd-min-len=i tsd-min-identity=f tsd-max-pvalue=f
                   tsd-pos-window=i tsd-scan-max=i tsd-null-match=f anchor-base=s decoy help)) or die usage();
if ($o{help}) { print usage(); exit 0 }
die usage() unless defined $o{fasta} && defined $o{flanks} && @ARGV == 1;
die "--anchor-base must be auto, yes or no\n" unless $o{'anchor-base'} =~ /^(?:auto|yes|no)\z/;
die "--tsd-min-identity must be above 0 and at most 1\n" unless $o{'tsd-min-identity'} > 0 && $o{'tsd-min-identity'} <= 1;
die "--tsd-null-match must be between 0 and 1\n" unless $o{'tsd-null-match'} > 0 && $o{'tsd-null-match'} < 1;
die "--tsd-scan-max must be at least 8\n" unless $o{'tsd-scan-max'} >= 8;

my $w = $o{'tsd-pos-window'};

# P(X >= m) for X ~ Binomial(k, p)
sub binom_tail {
    my ($k, $m, $p) = @_;
    my $q = 1 - $p;
    my ($t, $sum) = ($q ** $k, 0);
    for my $j (0 .. $k) {
        $sum += $t if $j >= $m;
        $t *= ($k - $j) / ($j + 1) * $p / $q;
    }
    return $sum > 1 ? 1 : $sum;
}

# adjusted p-value of a run with $len positions and $mat matches: [len][mat]
my @pt;
for my $len (1 .. 2 * $o{'tsd-scan-max'}) {
    for my $mat (0 .. $len) {
        $pt[$len][$mat] = binom_tail($len, $mat, $o{'tsd-null-match'}) * ($len + 1) * (2 * $w + 1);   # not clipped: used to order runs, clipped when printed
    }
}

# ---- inputs
my ($fh, $frows) = read_tsv($o{flanks});
need_columns($o{flanks}, $fh, qw(query chrom pos flank_size status window_seq));
my %flank = map { $_->{query} => $_ } @$frows;
my $seqs  = read_fasta($o{fasta});
my ($ih, $irows) = read_tsv($ARGV[0]);
need_columns($ARGV[0], $ih, qw(query n_elements));
my @todo = grep { $_->{n_elements} >= $o{'min-elements'} } @$irows;

# ---- anchor base: does the FASTA sequence start with the reference base at the position?
my ($n_check, $n_same) = (0, 0);
for my $r (@todo) {
    my $f = $flank{ $r->{query} };
    my $s = $seqs->{ $r->{query} };
    next unless $f && $f->{status} eq 'ok' && defined $s && length $s;
    $n_check++;
    $n_same++ if substr($s, 0, 1) eq substr($f->{window_seq}, $f->{flank_size} - 1, 1);
}
my $same_frac = $n_check ? $n_same / $n_check : 0;
my $anchor = $o{'anchor-base'} eq 'yes' ? 1
           : $o{'anchor-base'} eq 'no'  ? 0
           : ($same_frac >= 0.6 ? 1 : 0);
warn sprintf("find_tsd: the first insertion base equals the reference base at the position in %.0f%% of %d insertions: %s\n",
             100 * $same_frac, $n_check,
             $o{'anchor-base'} eq 'auto' ? ($anchor ? 'treated as an anchor base and removed' : 'no anchor base')
                                         : "--anchor-base $o{'anchor-base'} given")
    if $n_check;
warn "find_tsd: between 40% and 60%: the anchor test is inconclusive, set --anchor-base yes or no\n"
    if $o{'anchor-base'} eq 'auto' && $n_check && $same_frac > 0.4 && $same_frac < 0.6;

# ---- main
my $out = open_out($o{out});
print {$out} tsv_row(qw(query chrom pos ins_len status tsd_found tsd_len tsd_matches tsd_identity tsd_pvalue
                        tsd_delta tsd_start_offset tsd_capped tsd_seq tsd_seq2));

my (%status_n, %len_hist, %delta_hist);
my ($n_found, $n_capped, $n_long) = (0, 0, 0);

for my $idx (0 .. $#todo) {
    my $r  = $todo[$idx];
    my $q  = $r->{query};
    my $f  = $flank{$q};                                  # chrom, pos and status of this insertion
    my $fw = $o{decoy} ? $flank{ $todo[ ($idx + 1) % @todo ]{query} } : $f;    # the flank whose sequence is compared
    my $s = $seqs->{$q};
    my ($chrom, $pos) = $f ? @{$f}{qw(chrom pos)} : ('.', '.');
    my $n = defined $s ? length($s) - $anchor : 0;
    my $status = !defined $s || !length $s ? 'no_sequence'
               : (!$f || $f->{status} ne 'ok' || !$fw || $fw->{status} ne 'ok') ? 'no_flank'
               : 'ok';
    my $h  = $status eq 'ok' ? $fw->{flank_size} : 0;
    my $wl = $h - $w;                                     # flank length used on each side after shifting the junction
    my $k  = $o{'tsd-scan-max'};
    $k = int($n / 2) if $k > int($n / 2);
    $k = $wl         if $k > $wl;
    $status = 'too_short' if $status eq 'ok' && ($k < 8 || length($fw->{window_seq}) != 2 * $h);

    if ($status ne 'ok') {
        $status_n{$status}++;
        print {$out} tsv_row($q, $chrom, $pos, $n, $status, ('.') x 10);
        next;
    }
    $status_n{ok}++;

    my $ins = substr($s, $anchor);
    my $win = $fw->{window_seq};
    my (@cands, @best_any);
    for my $delta (-$w .. $w) {
        my $left  = substr($win, $h + $delta - $wl, $wl);
        my $right = substr($win, $h + $delta, $wl);
        my @a = split //, substr($left, $wl - $k, $k) . substr($ins, 0, $k);       # copy 1; index k is the first insertion base
        my @b = split //, substr($ins, $n - $k, $k) . substr($right, 0, $k);       # copy 2, one insertion length later
        my @pre = (0);
        for my $i (0 .. 2 * $k - 1) {
            push @pre, $pre[-1] + (($a[$i] eq $b[$i] && $a[$i] ne 'N') ? 1 : 0);
        }
        # runs [x, y] that touch the junction: x <= k (first insertion base) and y >= k-1 (last flank base)
        for my $x (0 .. $k) {
            for my $y (($x > $k - 1 ? $x : $k - 1) .. 2 * $k - 1) {
                my $len = $y - $x + 1;
                my $mat = $pre[$y + 1] - $pre[$x];
                next unless $mat;
                my $p = $pt[$len][$mat];
                my $keep = $p <= 0.2;
                my $best = !@best_any || $p < $best_any[0] || ($p == $best_any[0] && $mat > $best_any[7]);
                next unless $keep || $best;
                my $cand = [$p, $mat / $len, $len, abs($delta), $delta, $x, $y, $mat, \@a, \@b, $k];
                push @cands, $cand if $keep;
                @best_any = @$cand if $best;
            }
        }
    }

    my @pick;
    if (@cands) {
        my ($minp) = sort { $a <=> $b } map { $_->[0] } @cands;
        my ($chosen) = sort { $b->[1] <=> $a->[1] || $a->[2] <=> $b->[2] || $a->[3] <=> $b->[3] }
                       grep { $_->[0] <= 2 * $minp } @cands;
        @pick = @$chosen;
    }
    else {
        @pick = @best_any;
    }
    if (!@pick) {
        print {$out} tsv_row($q, $chrom, $pos, $n, 'ok', 0, ('.') x 9);
        next;
    }

    my ($p, $ident, $len, undef, $delta, $x, $y, $mat, $ar, $br, $kk) = @pick;
    my $found  = ($len >= $o{'tsd-min-len'} && $ident >= $o{'tsd-min-identity'} && $p <= $o{'tsd-max-pvalue'}) ? 1 : 0;
    my $capped = ($found && ($x == 0 || $y == 2 * $kk - 1)) ? 1 : 0;
    my $trim   = sub { my $t = join('', @_); length($t) > 60 ? substr($t, 0, 60) . '...' : $t };
    $n_capped++ if $capped;
    if ($found) {
        $n_found++;
        $n_long++ if $len >= 30;
        $len_hist{$len}++;
        $delta_hist{$delta}++;
    }
    print {$out} tsv_row($q, $chrom, $pos, $n, 'ok', $found, $len, $mat, sprintf('%.3f', $ident),
                         sprintf('%.2e', $p > 1 ? 1 : $p), $delta, $x - $kk, $capped,
                         $trim->(@{$ar}[$x .. $y]), $trim->(@{$br}[$x .. $y]));
}
close $out if defined $o{out};

warn sprintf("find_tsd:%s %d insertions analysed, TSD found in %d, %d of them 30 bp or longer (status: %s)\n",
             $o{decoy} ? ' DECOY RUN (unrelated flanks, the found fraction is the false-positive rate):' : '',
             scalar(@todo), $n_found, $n_long, join(', ', map { "$_ $status_n{$_}" } sort keys %status_n));
warn sprintf("find_tsd: %d runs reach the scan limit (long duplication of the neighbouring reference)\n", $n_capped)
    if $n_capped;
if (%len_hist) {
    my %bin;
    $bin{ $_ >= 30 ? '30+' : sprintf('%02d-%02d', 5 * int($_ / 5), 5 * int($_ / 5) + 4) } += $len_hist{$_} for keys %len_hist;
    warn "find_tsd: TSD length: " . join(', ', map { "$_ bp: $bin{$_}" } sort keys %bin) . "\n";
    warn "find_tsd: junction shift (delta): " . join(', ', map { "$_: $delta_hist{$_}" } sort { $a <=> $b } keys %delta_hist) . "\n";
}
