#!/usr/bin/env perl
# extract_flanks.pl - BED4 + indexed reference -> flanks.tsv
# One window of reference sequence centred on each insertion point; find_tsd.pl cuts the flanks it needs from it.
# Reads the reference through its .fai index, so samtools is only needed to create the index.
use strict;
use warnings;
use FindBin qw($RealBin);
use lib "$RealBin/../lib";
use Getopt::Long qw(GetOptions);
use Strongnose::Common qw(open_out tsv_row read_bed4);

sub usage {
    return <<'USAGE';
usage: extract_flanks.pl --reference ref.fa [options] insertions.bed4 > flanks.tsv

required:
  --reference FILE     reference FASTA with a .fai index next to it (samtools faidx ref.fa)

options (default):
  --out FILE           output table (STDOUT)
  --flank-size N       bp on each side of the insertion point (80); must cover
                       tsd-scan-max + tsd-pos-window of find_tsd.pl, 80 is enough for the defaults
  --help               this message

The insertion point P is bed_end of the BED4 (1-based position of the reference base just before the insertion,
the VCF POS). The window is the reference from P-N+1 to P+N, upper-cased.

Columns of flanks.tsv:
  query, chrom, pos (=P), flank_size
  status        ok | chrom_missing (not in the .fai) | edge (window falls outside the contig)
  window_seq    2*N bases; '.' when status is not ok
USAGE
}

my %o = ('flank-size' => 80);
GetOptions(\%o, qw(reference=s out=s flank-size=i help)) or die usage();
if ($o{help}) { print usage(); exit 0 }
die usage() unless defined $o{reference} && @ARGV == 1;
die "--flank-size must be a positive integer\n" unless $o{'flank-size'} > 0;

my $fai_file = "$o{reference}.fai";
die "$fai_file not found; create it with: samtools faidx $o{reference}\n" unless -e $fai_file;
my %fai;
open(my $fi, '<', $fai_file) or die "cannot open $fai_file: $!\n";
while (my $line = <$fi>) {
    chomp $line;
    my ($name, $len, $offset, $linebases, $linewidth) = split /\t/, $line;
    next unless defined $linewidth;
    $fai{$name} = [$len, $offset, $linebases, $linewidth];
}
close $fi;
open(my $ref, '<:raw', $o{reference}) or die "cannot open $o{reference}: $!\n";

# 1-based inclusive region
sub fetch {
    my ($c, $s, $e) = @_;
    my ($len, $offset, $lb, $lw) = @$c;
    my $start = $offset + int(($s - 1) / $lb) * $lw + (($s - 1) % $lb);
    my $nl    = int(($e - 1) / $lb) - int(($s - 1) / $lb);
    my $nbytes = ($e - $s + 1) + $nl * ($lw - $lb);
    seek($ref, $start, 0) or die "seek failed in $o{reference}\n";
    my $got = read($ref, my $buf, $nbytes);
    die "short read in $o{reference} at byte $start\n" unless defined $got && $got == $nbytes;
    $buf =~ tr/\r\n//d;
    return uc $buf;
}

my $bed = read_bed4($ARGV[0]);
my $out = open_out($o{out});
print {$out} tsv_row(qw(query chrom pos flank_size status window_seq));

my $h = $o{'flank-size'};
my %n;
for my $id (sort { $bed->{$a}[0] cmp $bed->{$b}[0] || $bed->{$a}[2] <=> $bed->{$b}[2] || $a cmp $b } keys %$bed) {
    my ($chrom, $start, $end) = @{ $bed->{$id} };
    my $c = $fai{$chrom};
    if (!$c) {
        $n{chrom_missing}++;
        print {$out} tsv_row($id, $chrom, $end, $h, 'chrom_missing', '.');
        next;
    }
    my ($s, $e) = ($end - $h + 1, $end + $h);
    if ($s < 1 || $e > $c->[0]) {
        $n{edge}++;
        print {$out} tsv_row($id, $chrom, $end, $h, 'edge', '.');
        next;
    }
    $n{ok}++;
    print {$out} tsv_row($id, $chrom, $end, $h, 'ok', fetch($c, $s, $e));
}
close $out if defined $o{out};

warn sprintf("extract_flanks: %d insertions, %d ok, %d at a contig edge, %d on contigs missing from the index\n",
             scalar(keys %$bed), $n{ok} || 0, $n{edge} || 0, $n{chrom_missing} || 0);
if ($n{chrom_missing}) {
    warn "extract_flanks: the BED and the reference must use the same contig names (chr1 vs 1)\n";
}
