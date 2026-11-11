#!/usr/bin/env perl
# summarize_events.pl - events.tsv (+ elements.tsv, discarded.tsv) -> plain-text summary used to calibrate the profiles
use strict;
use warnings;
use FindBin qw($RealBin);
use lib "$RealBin/../lib";
use Getopt::Long qw(GetOptions);
use Strongnose::Common qw(open_out read_tsv need_columns);

sub usage {
    return <<'USAGE';
usage: summarize_events.pl --events events.tsv [options] > summary.txt

required:
  --events FILE        from classify_events.pl

options (default):
  --elements FILE      elements.tsv; adds c_left (missing 3' consensus bases) to the quantile tables
  --discarded FILE     discarded.tsv; adds the reasons of the dropped insertions
  --out FILE           summary (STDOUT)
  --div-bin N          width of the divergence histogram bins, in percentage points (1)
  --top N              repeat names listed per profile (12)
  --help               this message

Sections:
  1. confidence per profile; dropped insertions by reason
  2. gates failed and flags per profile (young events)
  3. evidence cross table per profile: poly-A x TSD for young events that passed every gate
  4. quantiles per profile of the REFERENCE SET = young events with both poly-A and TSD found (TSD not capped),
     whatever their gates: length, coverage, divergence, c_left, TSD length and identity, poly-A length and purity.
     Poly-A and TSD do not depend on those values, so the set is a fair guide to the limits of the profile
     (len_min, len_max, min_cov, max_3end_left, tsd range, tsd_min_identity)
  5. divergence histograms (corrected divergence, div_corr_w) per tier for events with both poly-A and TSD found:
     the peak of the active tier shows how much ONT inflates the divergence (ont_div_correction) and where
     max_div_active / max_div_young can sit
  6. repeat names per profile with count, median divergence and share of high / probable
USAGE
}

my %o = ('div-bin' => 1, top => 12);
GetOptions(\%o, qw(events=s elements=s discarded=s out=s div-bin=f top=i help)) or die usage();
if ($o{help}) { print usage(); exit 0 }
die usage() unless defined $o{events};

my ($eh, $ev) = read_tsv($o{events});
need_columns($o{events}, $eh, qw(query profile elem_id repeat tier young ins_len cov div_corr_w polya_found polya_len polya_purity
                                 tsd_found tsd_len tsd_identity tsd_capped gates_failed confidence flags));
my %cleft;
if (defined $o{elements}) {
    my ($lh, $lr) = read_tsv($o{elements});
    need_columns($o{elements}, $lh, qw(elem_id c_left));
    $cleft{ $_->{elem_id} } = $_->{c_left} for @$lr;
}
my @dis;
if (defined $o{discarded}) {
    my ($dh, $dr) = read_tsv($o{discarded});
    need_columns($o{discarded}, $dh, qw(reason));
    @dis = @$dr;
}

my $out = open_out($o{out});
sub say_out { print {$out} @_, "\n" }
sub num { return defined $_[0] && $_[0] =~ /^-?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?$/ ? $_[0] : undef }
sub quant {
    my ($s, $q) = @_;
    return undef unless @$s;
    my $pos = $q * $#$s;
    my ($lo, $hi) = (int $pos, int($pos) + ($pos > int $pos ? 1 : 0));
    return $s->[$lo] + ($s->[$hi] - $s->[$lo]) * ($pos - $lo);
}
sub fmt { defined $_[0] ? ($_[0] == int($_[0]) ? sprintf('%d', $_[0]) : sprintf('%.3g', $_[0])) : 'NA' }
sub section { say_out(''); say_out("== $_[0]"); say_out('') }
sub table {
    my ($rows) = @_;
    my @w;
    for my $r (@$rows) { for my $i (0 .. $#$r) { my $l = length $r->[$i]; $w[$i] = $l if !defined $w[$i] || $l > $w[$i] } }
    for my $r (@$rows) { say_out(join('  ', map { sprintf('%-*s', $w[$_], $r->[$_]) } 0 .. $#$r)) }
}

my %by_prof;
push @{ $by_prof{ $_->{profile} } }, $_ for @$ev;
my @profiles = sort keys %by_prof;

# ---- 1
section('1. confidence per profile');
my @t = (['profile', 'events', 'high', 'probable', 'low', 'young']);
for my $p (@profiles) {
    my %c;
    my $young = 0;
    for my $e (@{ $by_prof{$p} }) { $c{ $e->{confidence} }++; $young++ if $e->{young} eq '1' }
    push @t, [$p, scalar @{ $by_prof{$p} }, map({ $c{$_} || 0 } qw(high probable low)), $young];
}
table(\@t);
if (@dis) {
    say_out('');
    my %r;
    $r{ $_->{reason} }++ for @dis;
    say_out('dropped insertions: ' . join(', ', map { "$_ $r{$_}" } sort { $r{$b} <=> $r{$a} } keys %r));
}

# ---- 2
section('2. gates failed and flags (young events)');
for my $p (@profiles) {
    my (%g, %f);
    my $n = 0;
    for my $e (grep { $_->{young} eq '1' } @{ $by_prof{$p} }) {
        $n++;
        $g{ $e->{gates_failed} }++;
        $f{$_}++ for grep { $_ ne '.' } split /,/, $e->{flags};
    }
    say_out("$p ($n young events)");
    say_out('  gates : ' . (join('; ', map { ($_ eq '.' ? 'none' : $_) . " $g{$_}" } sort { $g{$b} <=> $g{$a} } keys %g) || '-'));
    say_out('  flags : ' . (join('; ', map { "$_ $f{$_}" } sort { $f{$b} <=> $f{$a} } keys %f) || '-'));
}

# ---- 3
section('3. evidence of young events that passed every gate (poly-A x TSD, raw detection)');
@t = (['profile', 'n', 'polya+tsd', 'polya only', 'tsd only', 'neither', 'polya %', 'tsd %']);
for my $p (@profiles) {
    my @s = grep { $_->{young} eq '1' && $_->{gates_failed} eq '.' } @{ $by_prof{$p} };
    next unless @s;
    my ($b, $pa, $ts, $ne) = (0, 0, 0, 0);
    for my $e (@s) {
        my ($a, $c) = ($e->{polya_found} eq '1', $e->{tsd_found} eq '1' && $e->{tsd_capped} ne '1');
        if ($a && $c) { $b++ } elsif ($a) { $pa++ } elsif ($c) { $ts++ } else { $ne++ }
    }
    my $n = @s;
    push @t, [$p, $n, $b, $pa, $ts, $ne, sprintf('%.0f', 100 * ($b + $pa) / $n), sprintf('%.0f', 100 * ($b + $ts) / $n)];
}
table(\@t);

# ---- 4
section('4. reference set: young events with poly-A and TSD found (TSD not capped), whatever their gates');
my @metrics = (
    ['ins_len', 'ins_len'], ['cov', 'cov'], ['div_corr_w', 'div_corr_w'], ['c_left', undef],
    ['tsd_len', 'tsd_len'], ['tsd_identity', 'tsd_identity'], ['polya_len', 'polya_len'], ['polya_purity', 'polya_purity'],
);
for my $p (@profiles) {
    my @s = grep { $_->{young} eq '1' && $_->{polya_found} eq '1' && $_->{tsd_found} eq '1' && $_->{tsd_capped} ne '1' } @{ $by_prof{$p} };
    say_out("$p (n = " . scalar(@s) . ')');
    next unless @s;
    my @rows = (['  metric', 'min', 'p5', 'p25', 'median', 'p75', 'p95', 'max']);
    for my $m (@metrics) {
        my ($name, $col) = @$m;
        next if $name eq 'c_left' && !%cleft;
        my @v = sort { $a <=> $b } grep { defined } map { num($name eq 'c_left' ? $cleft{ $_->{elem_id} } : $_->{$col}) } @s;
        push @rows, ["  $name", map { fmt(quant(\@v, $_)) } (0, 0.05, 0.25, 0.5, 0.75, 0.95, 1)] if @v;
    }
    table(\@rows);
    say_out('');
}

# ---- 5
section("5. corrected divergence (div_corr_w), events with poly-A and TSD found, bins of $o{'div-bin'} percentage point(s)");
my %hist;
my $maxbin = 0;
for my $e (grep { $_->{polya_found} eq '1' && $_->{tsd_found} eq '1' && $_->{tsd_capped} ne '1' } @$ev) {
    my $d = num($e->{div_corr_w});
    next unless defined $d;
    my $b = int($d / $o{'div-bin'});
    $hist{ $e->{tier} }{$b}++;
    $maxbin = $b if $b > $maxbin;
}
my @tiers = grep { $hist{$_} } qw(active young generic intermediate old internal);
my @h = (['bin (%)', @tiers]);
for my $b (0 .. $maxbin) {
    push @h, [sprintf('%g-%g', $b * $o{'div-bin'}, ($b + 1) * $o{'div-bin'}), map { $hist{$_}{$b} || 0 } @tiers];
}
table(\@h);
for my $t (@tiers) {
    my @all = sort { $a <=> $b } map { num($_->{div_corr_w}) } grep { $_->{tier} eq $t && $_->{polya_found} eq '1' && $_->{tsd_found} eq '1' && $_->{tsd_capped} ne '1' && defined num($_->{div_corr_w}) } @$ev;
    say_out(sprintf("%-12s n=%d median %s p95 %s", $t, scalar @all, fmt(quant(\@all, 0.5)), fmt(quant(\@all, 0.95))));
}

# ---- 6
section('6. repeat names per profile');
for my $p (@profiles) {
    my %n;
    for my $e (@{ $by_prof{$p} }) {
        my $r = $n{ $e->{repeat} } ||= { n => 0, div => [], high => 0, probable => 0, tier => $e->{tier} };
        $r->{n}++;
        push @{ $r->{div} }, num($e->{div_corr_w}) if defined num($e->{div_corr_w});
        $r->{ $e->{confidence} }++ if $e->{confidence} eq 'high' || $e->{confidence} eq 'probable';
    }
    say_out($p);
    my @rows = (['  repeat', 'tier', 'n', 'median div', 'high', 'probable']);
    my @names = sort { $n{$b}{n} <=> $n{$a}{n} || $a cmp $b } keys %n;
    splice(@names, $o{top}) if @names > $o{top};
    for my $r (@names) {
        my @d = sort { $a <=> $b } @{ $n{$r}{div} };
        push @rows, ["  $r", $n{$r}{tier}, $n{$r}{n}, fmt(quant(\@d, 0.5)), $n{$r}{high}, $n{$r}{probable}];
    }
    table(\@rows);
    say_out('');
}
close $out if defined $o{out};
