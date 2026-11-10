#!/usr/bin/env perl
# annotate_hits.pl - hits.tsv -> hits_annotated.tsv
# Adds the repeat group and tier, the candidate/artifact class, the matching profiles,
# the ONT-corrected divergence and the profiles for which each hit is young.
# Nothing is dropped: filters are flags, so every later step can still report why an insertion was lost.
use strict;
use warnings;
use FindBin qw($RealBin);
use lib "$RealBin/../lib";
use Getopt::Long qw(GetOptions);
use Scalar::Util qw(looks_like_number);
use Strongnose::Common qw(open_out read_tsv tsv_row name_normalizer need_columns);

sub usage {
    return <<'USAGE';
usage: annotate_hits.pl --subfamilies subfamilies.tsv --profiles class_profiles.tsv [options] hits.tsv > hits_annotated.tsv

required:
  --subfamilies FILE            repeat name -> group and age tier
  --profiles FILE               event profiles (class match, young tiers, ...)

options (default):
  --out FILE                    output table (STDOUT)
  --ont-div-correction N        percentage points subtracted from "div" before the age test (0.0)
  --max-div-active N            corrected div ceiling for tier "active" (8.0)
  --max-div-young N             corrected div ceiling for tier "young" (10.0)
  --max-div-generic N           corrected div ceiling for tier "generic" (6.0)
  --[no-]name-case-insensitive  fold case when matching repeat names (on)
  --name-strip-suffixes LIST    comma list of suffixes removed before matching ("_5,_5end,_3end")
  --[no-]name-int-equivalent    treat a trailing "-int" and "I" as the same (on)
  --artifact-classes LIST       class prefixes flagged as artifacts ("Simple_repeat,Low_complexity,Satellite,Unknown")
  --help                        this message

Columns added to hits.tsv:
  name_norm        normalized repeat name used for the lookup
  group, tier      from subfamilies.tsv; '.' and "old" when the name is not listed
  tier_source      table | default
  hit_kind         candidate = class matches an enabled profile (and the group agrees); artifact = artifact class; other
  profiles         enabled profiles the hit belongs to ('.' when none)
  div_corr         max(0, div - ont-div-correction)
  young_profiles   profiles for which the tier counts as young and div_corr is within the ceiling ('.' when none)
  div_limit        the ceiling that applies to the hit (tier ceiling or per-name max_div); '.' when there is none

Ceilings come from --max-div-* by tier; a max_div in subfamilies.tsv overrides the ceiling for that name.
Tiers without a ceiling (intermediate, old, internal) are never young. A profile's rm_family_match is
matched against the whole rm_family string.
USAGE
}

my %o = (
    'ont-div-correction'    => 0.0,
    'max-div-active'        => 8.0,
    'max-div-young'         => 10.0,
    'max-div-generic'       => 6.0,
    'name-case-insensitive' => 1,
    'name-strip-suffixes'   => '_5,_5end,_3end',
    'name-int-equivalent'   => 1,
    'artifact-classes'      => 'Simple_repeat,Low_complexity,Satellite,Unknown',
);
GetOptions(\%o, qw(subfamilies=s profiles=s out=s ont-div-correction=f max-div-active=f
                   max-div-young=f max-div-generic=f name-case-insensitive! name-strip-suffixes=s
                   name-int-equivalent! artifact-classes=s help))
    or die usage();
if ($o{help}) { print usage(); exit 0 }
die usage() unless defined $o{subfamilies} && defined $o{profiles} && @ARGV == 1;

my %TIER = map { $_ => 1 } qw(active young intermediate old internal generic);
my %ceiling = (
    active  => $o{'max-div-active'},
    young   => $o{'max-div-young'},
    generic => $o{'max-div-generic'},
);

my $norm = name_normalizer(
    case_insensitive => $o{'name-case-insensitive'},
    strip_suffixes   => $o{'name-strip-suffixes'},
    int_equivalent   => $o{'name-int-equivalent'},
);

# ---- subfamilies.tsv
my ($sh, $srows) = read_tsv($o{subfamilies});
need_columns($o{subfamilies}, $sh, qw(repeat_name group tier rm_family max_div));
my %table;
for my $r (@$srows) {
    my $name = $r->{repeat_name};
    die "$o{subfamilies}: unknown tier '$r->{tier}' for $name (use: " . join(' ', sort keys %TIER) . ")\n"
        unless $TIER{ $r->{tier} };
    die "$o{subfamilies}: max_div '$r->{max_div}' of $name is not a number\n"
        if $r->{max_div} ne '' && !looks_like_number($r->{max_div});
    my $key = $norm->($name);
    die "$o{subfamilies}: $name and $table{$key}{repeat_name} become the same name '$key'\n"
        if exists $table{$key};
    $table{$key} = $r;
}

# ---- class_profiles.tsv
my ($ph, $prows) = read_tsv($o{profiles});
need_columns($o{profiles}, $ph, qw(profile group rm_family_match enabled young_tiers));
my @profiles;
for my $r (@$prows) {
    die "$o{profiles}: enabled of $r->{profile} must be yes or no\n" unless $r->{enabled} =~ /^(?:yes|no)\z/;
    next if $r->{enabled} eq 'no';
    my $re = eval { qr/\A(?:$r->{rm_family_match})\z/ }
        or die "$o{profiles}: bad rm_family_match '$r->{rm_family_match}' in $r->{profile}\n";
    my %young = map { $_ => 1 } grep { length } split /,/, $r->{young_tiers};
    for my $t (keys %young) {
        die "$o{profiles}: unknown tier '$t' in young_tiers of $r->{profile}\n" unless $TIER{$t};
    }
    push @profiles, { id => $r->{profile}, group => $r->{group}, re => $re, young => \%young };
}
warn "annotate_hits: no enabled profile in $o{profiles}\n" unless @profiles;

my @artifact = grep { length } map { my $s = $_; $s =~ s/^\s+|\s+\z//g; $s } split /,/, $o{'artifact-classes'};

# ---- hits
my ($hh, $hits) = read_tsv($ARGV[0]);
need_columns($ARGV[0], $hh, qw(query repeat rm_family div));
my @new_cols = qw(name_norm group tier tier_source hit_kind profiles div_corr young_profiles div_limit);

my $out = open_out($o{out});
print {$out} tsv_row(@$hh, @new_cols);

my (%kind_n, %tier_n, %unlisted, %mismatch, %young_query, %seen_query);
my $n_young = 0;

for my $r (@$hits) {
    my $key = $norm->($r->{repeat});
    my $sf  = $table{$key};
    my ($group, $tier, $source, $override) = $sf
        ? ($sf->{group}, $sf->{tier}, 'table', $sf->{max_div})
        : ('.', 'old', 'default', '');

    my $family = $r->{rm_family};
    my @belong = grep { $family =~ $_->{re} && ($group eq '.' || $group eq $_->{group}) } @profiles;

    my $kind = @belong ? 'candidate'
             : (grep { index($family, $_) == 0 } @artifact) ? 'artifact'
             : 'other';

    my $corr = $r->{div} - $o{'ont-div-correction'};
    $corr = 0 if $corr < 0;

    my $limit = $kind eq 'candidate' ? ($override ne '' ? $override : $ceiling{$tier}) : undef;
    my @young;
    if (defined $limit && $corr <= $limit) {
        @young = map { $_->{id} } grep { $_->{young}{$tier} } @belong;
    }

    $seen_query{ $r->{query} } = 1;
    $kind_n{$kind}++;
    if ($kind eq 'candidate') {
        $tier_n{$tier}++;
        $unlisted{ $r->{repeat} }++ unless $sf;
    }
    if ($sf && $sf->{rm_family} ne '' && $sf->{rm_family} ne $family) {
        $mismatch{"$r->{repeat}: table says $sf->{rm_family}, hit says $family"}++;
    }
    if (@young) { $n_young++; $young_query{ $r->{query} } = 1 }

    print {$out} tsv_row(
        @{$r}{@$hh}, $key, $group, $tier, $source, $kind,
        (@belong ? join(',', map { $_->{id} } @belong) : '.'),
        sprintf('%.2f', $corr),
        (@young ? join(',', @young) : '.'),
        (defined $limit ? $limit : '.'),
    );
}
close $out if defined $o{out};

# ---- summary on STDERR
warn sprintf("annotate_hits: %d hits in %d insertions (candidate %d, artifact %d, other %d)\n",
             scalar(@$hits), scalar(keys %seen_query),
             map { $kind_n{$_} || 0 } qw(candidate artifact other));
warn sprintf("annotate_hits: candidate hits by tier: %s\n",
             join(', ', map { "$_ $tier_n{$_}" } sort keys %tier_n)) if %tier_n;
warn sprintf("annotate_hits: %d young hits in %d insertions\n", $n_young, scalar(keys %young_query));

if (%unlisted) {
    my @top = (sort { $unlisted{$b} <=> $unlisted{$a} || $a cmp $b } keys %unlisted)[0 .. (keys(%unlisted) > 10 ? 9 : keys(%unlisted) - 1)];
    warn sprintf("annotate_hits: %d candidate names are not in subfamilies.tsv and were treated as old: %s%s\n",
                 scalar(keys %unlisted), join(', ', map { "$_ ($unlisted{$_})" } @top),
                 keys(%unlisted) > 10 ? ', ...' : '');
}
if (%mismatch) {
    warn sprintf("annotate_hits: %d names whose class differs from subfamilies.tsv:\n", scalar(keys %mismatch));
    warn "  $_\n" for (sort keys %mismatch)[0 .. (keys(%mismatch) > 5 ? 4 : keys(%mismatch) - 1)];
}
