#!/usr/bin/env perl
# parse_rm_out.pl - RepeatMasker .out -> one TSV row per hit (hits.tsv)
use strict;
use warnings;
use FindBin qw($RealBin);
use lib "$RealBin/../lib";
use Getopt::Long qw(GetOptions);
use Strongnose::Common qw(open_in open_out parse_rm_line read_bed4 tsv_row hit_columns);

sub usage {
    return <<'USAGE';
usage: parse_rm_out.pl [--bed insertions.bed4] [--out hits.tsv] sample.fa.out [more.out ...]

  --bed FILE   BED4 with the Sniffles IDs; adds chrom, bed_start, bed_end to every hit
  --out FILE   output table (default: STDOUT)
  --help       this message

Input files may be "-" (STDIN) or gzipped (*.gz). Several .out files are concatenated.

Columns of hits.tsv:
  hit_row        order of the hit in the input (unique)
  query          Sniffles ID
  chrom, bed_start, bed_end   from the BED ('.' without --bed or when the ID is missing from it)
  q_start, q_end, q_len       position in the insertion; q_len = q_end + (left)
  strand         + or -  (RepeatMasker "C" becomes "-")
  score, div, del, ins        SW score and the three percentages
  repeat, rm_family           repeat name and class/family, e.g. AluYa5, SINE/Alu
  c_start, c_end, c_len       position in the repeat consensus; always c_start < c_end; c_len = c_end + c_left
  c_left         consensus bases missing after c_end (3' side), the same for + and - hits
  rm_id          RepeatMasker element ID (fragments of one element share it)
  overlapped     1 when the hit carries the "*" flag (overlapped by a higher-scoring hit)
USAGE
}

my ($bed_file, $out_file, $help);
GetOptions('bed=s' => \$bed_file, 'out=s' => \$out_file, 'help' => \$help)
    or die usage();
if ($help) { print usage(); exit 0 }

my $bed = defined $bed_file ? read_bed4($bed_file) : undef;
my $out = open_out($out_file);
print {$out} tsv_row(hit_columns());

my ($n_hits, $n_bad) = (0, 0);
my (%seen, %no_bed, %clen);

for my $file (@ARGV ? @ARGV : '-') {
    my $fh = open_in($file);
    while (my $line = <$fh>) {
        chomp $line;
        my $h = eval { parse_rm_line($line) };
        if ($@) {
            chomp(my $err = $@);
            warn "$file:$.: skipped ($err)\n";
            $n_bad++;
            next;
        }
        next unless $h;                            # header or blank line

        $n_hits++;
        $seen{ $h->{query} } = 1;

        my ($chrom, $bed_start, $bed_end) = ('.', '.', '.');
        if ($bed) {
            if (my $coords = $bed->{ $h->{query} }) { ($chrom, $bed_start, $bed_end) = @$coords }
            else                                    { $no_bed{ $h->{query} } = 1 }
        }
        # Simple repeats and low complexity have no fixed consensus: their "length" is the matched array
        $clen{ $h->{repeat} }{ $h->{c_len} }++ unless $h->{rm_family} =~ /^(?:Simple_repeat|Low_complexity)/;

        print {$out} tsv_row(
            $n_hits, $h->{query}, $chrom, $bed_start, $bed_end,
            @{$h}{qw(q_start q_end q_len strand score div del ins repeat rm_family
                     c_start c_end c_len c_left rm_id overlapped)}
        );
    }
    close $fh unless $file eq '-';
}
close $out if defined $out_file && $out_file ne '-';

# Summary on STDERR
warn sprintf("parse_rm_out: %d hits, %d insertions, %d lines skipped\n",
             $n_hits, scalar(keys %seen), $n_bad);
warn sprintf("parse_rm_out: %d insertions are not in the BED\n", scalar(keys %no_bed))
    if $bed && %no_bed;

# Sanity check of the consensus columns: one repeat name must have one consensus length.
my @inconsistent = sort grep { keys %{ $clen{$_} } > 1 } keys %clen;
if (@inconsistent) {
    warn sprintf("parse_rm_out: %d repeats show more than one consensus length (check the strand columns):\n",
                 scalar @inconsistent);
    for my $name (@inconsistent[0 .. ($#inconsistent < 9 ? $#inconsistent : 9)]) {
        warn "  $name: " . join(', ', sort { $a <=> $b } keys %{ $clen{$name} }) . "\n";
    }
}

if ($n_hits == 0) {
    warn "parse_rm_out: no hit lines found\n";
    exit 1;
}
