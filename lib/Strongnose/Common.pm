package Strongnose::Common;

# Strongnose::Common - helpers shared by the Strongnose scripts
# (Sniffles insertions -> RepeatMasker -> mobile DNA events).
#
#   open_in($path)        read handle; "-" is STDIN, "*.gz" is read through gzip
#   open_out($path)       write handle; undef or "-" is STDOUT
#   tsv_row(@fields)      one tab-separated line, newline included
#   read_tsv($path)       (\@header, \@rows) with each row a hashref; "#" comment lines and blank lines are skipped;
#                         rows shorter than the header are padded with empty fields, longer rows are an error
#   read_bed4($path)      { sniffles_id => [chrom, start, end] }
#   parse_rm_line($line)  hashref for a RepeatMasker .out hit line; undef for header/blank lines; dies on a malformed hit line
#   hit_columns()         column order of hits.tsv
#   name_normalizer(%o)   coderef mapping a repeat name to the key used for table lookups (see below)
#   need_columns($file, \@header, @names)   die unless every name is a column of the table
#   parse_params($str)    "k=v;k2=v2" (the rule_params column) -> hashref; a bare "k" is 1
#   merge_intervals(@iv)  sort and join [start,end] pairs (inclusive, 1-based) that overlap or touch
#   intervals_len(@iv)    number of bases covered by the union of the pairs
#   overlap_len(s1,e1,s2,e2)  bases shared by two inclusive intervals (0 when disjoint)
#   read_fasta($path)     { id => SEQUENCE }; id is the first word of the header, sequence upper-cased and unwrapped
#   revcomp($seq)         reverse complement (ACGTN, any case)

use strict;
use warnings;
use Exporter 'import';

our $VERSION   = '0.01';
our @EXPORT_OK = qw(open_in open_out tsv_row read_tsv read_bed4 parse_rm_line hit_columns name_normalizer
                    need_columns parse_params merge_intervals intervals_len overlap_len read_fasta revcomp);

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
        die "$path:$.: expected at most " . scalar(@header) . " columns, got " . scalar(@f) . "\n"
            if @f > @header;
        push @f, ('') x (@header - @f);          # short rows: missing trailing fields are empty (hand-edited files)
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

# ------------------------------------------------------------------ names

# name_normalizer(case_insensitive => 1, strip_suffixes => '_5,_5end,_3end', int_equivalent => 1)
# Returns a coderef that maps a repeat name to the key used to look it up, so that library variants
# of one name meet: L1MEc / L1MEC_5, HERVK9-int / HERVK9I.
#   case_insensitive  fold case
#   strip_suffixes    comma list; the first one the name ends with is removed (positional consensus fragments)
#   int_equivalent    a trailing "-int" becomes "I"
sub name_normalizer {
    my (%o) = @_;
    my $ci    = exists $o{case_insensitive} ? $o{case_insensitive} : 1;
    my $int   = exists $o{int_equivalent}   ? $o{int_equivalent}   : 1;
    my $strip = defined $o{strip_suffixes}  ? $o{strip_suffixes}   : '_5,_5end,_3end';
    my @suffix = sort { length($b) <=> length($a) }
                 grep { length }
                 map  { my $s = $_; $s =~ s/^\s+|\s+\z//g; $ci ? lc $s : $s }
                 split /,/, $strip;
    return sub {
        my ($name) = @_;
        my $n = $ci ? lc $name : $name;
        for my $s (@suffix) {
            if (length($n) > length($s) && substr($n, -length($s)) eq $s) {
                substr($n, -length($s)) = '';
                last;
            }
        }
        $n =~ s/-int\z/$ci ? 'i' : 'I'/e if $int;
        return $n;
    };
}

# ----------------------------------------------------------- tables, params

sub need_columns {
    my ($file, $header, @need) = @_;
    my %have = map { $_ => 1 } @$header;
    my @missing = grep { !$have{$_} } @need;
    die "$file: missing column(s): " . join(', ', @missing) . "\n" if @missing;
    return 1;
}

sub parse_params {
    my ($str) = @_;
    my %p;
    return \%p unless defined $str && length $str;
    for my $kv (split /;/, $str) {
        $kv =~ s/^\s+|\s+\z//g;
        next unless length $kv;
        my ($k, $v) = split /=/, $kv, 2;
        $p{$k} = defined $v ? $v : 1;
    }
    return \%p;
}

# ------------------------------------------------------------- intervals

sub merge_intervals {
    my @iv = sort { $a->[0] <=> $b->[0] || $a->[1] <=> $b->[1] }
             grep { $_->[1] >= $_->[0] } @_;
    my @out;
    for my $i (@iv) {
        if (@out && $i->[0] <= $out[-1][1] + 1) {
            $out[-1][1] = $i->[1] if $i->[1] > $out[-1][1];
        }
        else {
            push @out, [@$i];
        }
    }
    return @out;
}

sub intervals_len {
    my $n = 0;
    $n += $_->[1] - $_->[0] + 1 for merge_intervals(@_);
    return $n;
}

sub overlap_len {
    my ($s1, $e1, $s2, $e2) = @_;
    my $s = $s1 > $s2 ? $s1 : $s2;
    my $e = $e1 < $e2 ? $e1 : $e2;
    return $e >= $s ? $e - $s + 1 : 0;
}

# ----------------------------------------------------------------- sequences

sub read_fasta {
    my ($path) = @_;
    my $fh = open_in($path);
    my (%seq, $id);
    while (my $line = <$fh>) {
        $line =~ s/\r?\n\z//;
        if ($line =~ /^>(\S+)/) {
            $id = $1;
            die "$path:$.: duplicate sequence ID $id\n" if exists $seq{$id};
            $seq{$id} = '';
        }
        elsif (defined $id) {
            $line =~ s/\s+//g;
            $seq{$id} .= uc $line;
        }
    }
    close $fh unless $path eq '-';
    return \%seq;
}

sub revcomp {
    my ($s) = @_;
    $s = reverse $s;
    $s =~ tr/ACGTNacgtn/TGCANtgcan/;
    return $s;
}

1;
