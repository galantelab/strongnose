package Strongnose::Common;

# Strongnose::Common - helpers shared by the Strongnose scripts
# (Sniffles insertions -> RepeatMasker -> mobile DNA events).
#
#   open_in($path)        read handle; "-" is STDIN, "*.gz" is read through gzip
#   open_out($path)       write handle; undef or "-" is STDOUT
#   tsv_row(@fields)      one tab-separated line, newline included
#   read_tsv($path)       (\@header, \@rows) with each row a hashref; "#" comment lines and blank lines are skipped
#   read_bed4($path)      { sniffles_id => [chrom, start, end] }
#   parse_rm_line($line)  hashref for a RepeatMasker .out hit line; undef for header/blank lines; dies on a malformed hit line
#   hit_columns()         column order of hits.tsv

use strict;
use warnings;
use Exporter 'import';

our $VERSION   = '0.01';
our @EXPORT_OK = qw(open_in open_out tsv_row read_tsv read_bed4 parse_rm_line hit_columns);

# ---------------------------------------------------------------- files

sub open_in {
    my ($path) = @_;
    die "open_in: no path given\n" unless defined $path;
    return \*STDIN if $path eq '-';
    my $fh;
    if ($path =~ /\.gz\z/) {
        open($fh, '-|', 'gzip', '-dc', '--', $path)
            or die "cannot read $path through gzip: $!\n";
    }
    else {
        open($fh, '<', $path) or die "cannot open $path: $!\n";
    }
    return $fh;
}

sub open_out {
    my ($path) = @_;
    return \*STDOUT if !defined $path || $path eq '-';
    open(my $fh, '>', $path) or die "cannot write $path: $!\n";
    return $fh;
}

# ----------------------------------------------------------------- tables

sub tsv_row {
    return join("\t", @_) . "\n";
}

sub read_tsv {
    my ($path) = @_;
    my $fh = open_in($path);
    my (@header, @rows);
    while (my $line = <$fh>) {
        $line =~ s/\r?\n\z//;
        next if $line eq '' || $line =~ /^#/;
        my @f = split /\t/, $line, -1;          # -1 keeps empty trailing fields
        if (!@header) {
            @header = @f;
            next;
        }
        die "$path:$.: expected " . scalar(@header) . " columns, got " . scalar(@f) . "\n"
            if @f != @header;
        my %row;
        @row{@header} = @f;
        push @rows, \%row;
    }
    close $fh unless $path eq '-';
    return (\@header, \@rows);
}

sub read_bed4 {
    my ($path) = @_;
    my $fh = open_in($path);
    my %bed;
    while (my $line = <$fh>) {
        $line =~ s/\r?\n\z//;
        next if $line eq '' || $line =~ /^(?:#|track|browser)/;
        my ($chrom, $start, $end, $id) = split /\t/, $line;
        die "$path:$.: BED4 needs 4 tab-separated columns\n"
            unless defined $id && $id ne '';
        if (exists $bed{$id}) {
            warn "$path:$.: duplicate ID $id, keeping the first\n";
            next;
        }
        $bed{$id} = [$chrom, $start, $end];
    }
    close $fh unless $path eq '-';
    return \%bed;
}

# ------------------------------------------------------------ RepeatMasker

sub hit_columns {
    return qw(hit_row query chrom bed_start bed_end q_start q_end q_len strand
              score div del ins repeat rm_family c_start c_end c_len c_left
              rm_id overlapped);
}

# .out hit line:
#   score div del ins query q_begin q_end (q_left) strand repeat class/family  <3 repeat columns>  ID [*]
# The repeat columns are "begin end (left)" for "+" hits and "(left) end begin" for "C" hits.
# In both cases (left) is the number of consensus bases after the end, i.e. missing at the 3' side.
sub parse_rm_line {
    my ($line) = @_;
    return unless defined $line;
    $line =~ s/^\s+//;
    return unless $line =~ /^\d/;               # header and blank lines
    $line =~ s/\s+\z//;

    my @f = split ' ', $line;
    die "expected at least 15 fields, got " . scalar(@f) . ": $line\n" if @f < 15;

    my ($score, $div, $del, $ins, $query, $q_start, $q_end, $q_left, $strand,
        $repeat, $rm_family, $r1, $r2, $r3, $rm_id, $star) = @f;

    tr/()//d for $q_left, $r1, $r2, $r3;

    die "unknown strand '$strand': $line\n" unless $strand eq '+' || $strand eq 'C';
    for my $n ($score, $q_start, $q_end, $q_left, $r1, $r2, $r3, $rm_id) {
        die "non-integer field '$n': $line\n" unless $n =~ /^\d+\z/;
    }
    for my $n ($div, $del, $ins) {
        die "non-numeric field '$n': $line\n" unless $n =~ /^\d+(?:\.\d+)?\z/;
    }

    my ($c_start, $c_end, $c_left) = $strand eq '+' ? ($r1, $r2, $r3) : ($r3, $r2, $r1);

    die "query start after end: $line\n"         if $q_start > $q_end;
    die "consensus start after end: $line\n"     if $c_start > $c_end;

    return {
        query      => $query,
        q_start    => $q_start,
        q_end      => $q_end,
        q_len      => $q_end + $q_left,
        strand     => ($strand eq 'C' ? '-' : '+'),
        score      => $score,
        div        => $div,
        del        => $del,
        ins        => $ins,
        repeat     => $repeat,
        rm_family  => $rm_family,
        c_start    => $c_start,
        c_end      => $c_end,
        c_len      => $c_end + $c_left,
        c_left     => $c_left,
        rm_id      => $rm_id,
        overlapped => (defined $star && $star eq '*') ? 1 : 0,
    };
}

1;
