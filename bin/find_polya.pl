#!/usr/bin/env perl
# find_polya.pl - elements.tsv + insertion FASTA -> polya.tsv
# Looks for an A-rich tail at the 3' end of every element, strand-aware: after q_end for "+" elements,
# before q_start (as a T-rich stretch on the insertion strand) for "-" elements.
use strict;
use warnings;
use FindBin qw($RealBin);
use lib "$RealBin/../lib";
use Getopt::Long qw(GetOptions);
use Strongnose::Common qw(open_out read_tsv tsv_row need_columns read_fasta revcomp);

sub usage {
    return <<'USAGE';
usage: find_polya.pl --fasta insertions.fa [options] elements.tsv > polya.tsv

required:
  --fasta FILE                  insertion sequences (>Sniffles ID ...), the FASTA RepeatMasker was run on

options (default):
  --out FILE                    output table (STDOUT)
  --polya-search-window N       bp searched past the element 3' end (40); the window grows while the tail reaches its end
  --polya-inside-margin N       bp before the element 3' end also searched: the hit often swallows the first
                                bases of the tail (20)
  --polya-max-gap N             the tail must start at most this many bp after the element end (10)
  --polya-min-len N             minimum tail length in bp (10)
  --polya-min-purity F          minimum A fraction of the tail, 0-1 (0.85); ONT homopolymer errors are tolerated
  --help                        this message

Method: the region around the 3' end is read in the orientation of the element (reverse-complemented for "-"
elements, so the tail always shows as A). The best A-rich segment is found with a maximum-scoring-segment scan
(A = +1, any other base = -purity/(1-purity)); it may start up to --polya-max-gap bp after the element end.
It counts as a tail when its length and purity pass the minimums. Tails further away (3' transductions) are not found.

Columns of polya.tsv, one row per element of elements.tsv:
  query, elem_id, strand
  tail_side        right (+ elements) or left (- elements)
  search_start, search_end   region searched, in insertion coordinates
  polya_found      1 / 0; '.' when the sequence is missing from the FASTA
  polya_start, polya_end, polya_len, polya_purity   best segment, insertion coordinates; reported also when found=0
  polya_gap        bases between the element end and the tail start (negative: the tail starts inside the element)
  bases_after_tail bases left in the insertion beyond the outer end of the tail (0: the tail ends the insertion)
  tail_seq         the segment read as A-rich, first 50 bases
USAGE
}

my %o = (
    'polya-search-window' => 40,
    'polya-inside-margin' => 20,
    'polya-max-gap'       => 10,
    'polya-min-len'       => 10,
    'polya-min-purity'    => 0.85,
);
GetOptions(\%o, qw(fasta=s out=s polya-search-window=i polya-inside-margin=i polya-max-gap=i
                   polya-min-len=i polya-min-purity=f help)) or die usage();
if ($o{help}) { print usage(); exit 0 }
die usage() unless defined $o{fasta} && @ARGV == 1;
die "--polya-min-purity must be above 0 and at most 1\n" unless $o{'polya-min-purity'} > 0 && $o{'polya-min-purity'} <= 1;

my $penalty = $o{'polya-min-purity'} >= 1 ? 1000 : $o{'polya-min-purity'} / (1 - $o{'polya-min-purity'});

# Best A-rich segment of $r whose first base is at an index <= $limit. Returns (start, end) or ().
sub best_tail {
    my ($r, $limit) = @_;
    my ($best, $best_s, $best_e) = (0);
    my ($cur, $cs) = (0, undef);
    for my $i (0 .. length($r) - 1) {
        my $is_a = substr($r, $i, 1) eq 'A';
        if (!defined $cs) {
            next unless $is_a && $i <= $limit;
            ($cs, $cur) = ($i, 1);
        }
        else {
            $cur += $is_a ? 1 : -$penalty;
            if ($cur <= 0) { ($cs, $cur) = (undef, 0); next }
        }
        ($best, $best_s, $best_e) = ($cur, $cs, $i) if $cur > $best;
    }
    return defined $best_s ? ($best_s, $best_e) : ();
}

my ($eh, $elems) = read_tsv($ARGV[0]);
need_columns($ARGV[0], $eh, qw(query elem_id strand q_start q_end q_len));
my $seqs = read_fasta($o{fasta});

my $out = open_out($o{out});
print {$out} tsv_row(qw(query elem_id strand tail_side search_start search_end polya_found polya_start polya_end
                        polya_len polya_purity polya_gap bases_after_tail tail_seq));

my ($n, $n_found, $n_missing, $n_len_diff) = (0, 0, 0, 0);
my (%found_by_strand, %n_by_strand, %len_warned);

for my $e (@$elems) {
    $n++;
    $n_by_strand{ $e->{strand} }++;
    my $side = $e->{strand} eq '+' ? 'right' : 'left';
    my $seq  = $seqs->{ $e->{query} };
    if (!defined $seq) {
        $n_missing++;
        print {$out} tsv_row(@{$e}{qw(query elem_id strand)}, $side, ('.') x 10);
        next;
    }
    my $len = length $seq;
    if ($len != $e->{q_len} && !$len_warned{ $e->{query} }++) {
        $n_len_diff++;
        warn "find_polya: $e->{query}: FASTA length $len differs from q_len $e->{q_len}; using the FASTA length\n"
            if $n_len_diff <= 5;
    }

    # The window grows while the tail touches its far end, so a tail longer than the window is measured whole
    # (otherwise the rest of the tail would be counted as bases after it).
    my ($rs, $re, $r, $b, $s, $t);                 # region start/end (insertion coordinates), the string, index of the element's last base
    my $win = $o{'polya-search-window'};
    while (1) {
        my $can_grow;
        if ($e->{strand} eq '+') {
            $rs = $e->{q_end} - $o{'polya-inside-margin'} + 1;
            $rs = 1 if $rs < 1;
            $re = $e->{q_end} + $win;
            $can_grow = $re < $len;
            $re = $len if $re > $len;
            $r  = substr($seq, $rs - 1, $re - $rs + 1);
            $b  = $e->{q_end} - $rs;
        }
        else {
            $rs = $e->{q_start} - $win;
            $can_grow = $rs > 1;
            $rs = 1 if $rs < 1;
            $re = $e->{q_start} + $o{'polya-inside-margin'} - 1;
            $re = $len if $re > $len;
            $r  = revcomp(substr($seq, $rs - 1, $re - $rs + 1));
            $b  = $re - $e->{q_start};
        }
        ($s, $t) = best_tail($r, $b + 1 + $o{'polya-max-gap'});
        last unless defined $s && $t == length($r) - 1 && $can_grow;
        $win *= 2;
    }

    if (!defined $s) {
        print {$out} tsv_row(@{$e}{qw(query elem_id strand)}, $side, $rs, $re, 0, ('.') x 7);
        next;
    }

    my $tail   = substr($r, $s, $t - $s + 1);
    my $tlen   = length $tail;
    my $purity = ($tail =~ tr/A//) / $tlen;
    my ($lo, $hi) = $e->{strand} eq '+' ? ($rs + $s, $rs + $t) : ($re - $t, $re - $s);
    my $found  = ($tlen >= $o{'polya-min-len'} && $purity >= $o{'polya-min-purity'}) ? 1 : 0;
    my $after  = $e->{strand} eq '+' ? $len - $hi : $lo - 1;
    my $shown  = $tlen > 50 ? substr($tail, 0, 50) . '...' : $tail;

    if ($found) { $n_found++; $found_by_strand{ $e->{strand} }++ }
    print {$out} tsv_row(@{$e}{qw(query elem_id strand)}, $side, $rs, $re, $found, $lo, $hi, $tlen,
                         sprintf('%.3f', $purity), $s - $b - 1, $after, $shown);
}
close $out if defined $o{out};

warn sprintf("find_polya: %d elements, tail found for %d (+ %d of %d, - %d of %d)\n",
             $n, $n_found, $found_by_strand{'+'} || 0, $n_by_strand{'+'} || 0,
             $found_by_strand{'-'} || 0, $n_by_strand{'-'} || 0);
warn sprintf("find_polya: %d elements have no sequence in %s\n", $n_missing, $o{fasta}) if $n_missing;
warn sprintf("find_polya: %d insertions have a FASTA length different from q_len\n", $n_len_diff) if $n_len_diff;
