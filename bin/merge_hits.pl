#!/usr/bin/env perl
# merge_hits.pl - hits_annotated.tsv -> elements.tsv + insertions.tsv
# Per insertion: resolves overlapping hits, joins fragments into elements, applies the absorb rules
# of the profiles, ranks the elements and measures how much of the insertion is artifact-like.
use strict;
use warnings;
use FindBin qw($RealBin);
use lib "$RealBin/../lib";
use Getopt::Long qw(GetOptions);
use Strongnose::Common qw(open_out read_tsv tsv_row parse_params intervals_len overlap_len);

sub usage {
    return <<'USAGE';
usage: merge_hits.pl --profiles class_profiles.tsv --insertions-out insertions.tsv [options] hits_annotated.tsv > elements.tsv

required:
  --profiles FILE               class_profiles.tsv (absorb_inside_span rules)
  --insertions-out FILE         per-insertion summary (every insertion, also those without elements)

options (default):
  --out FILE                    elements table (STDOUT)
  --alternatives-out FILE       candidate hits that were not used, with the reason
  --[no-]merge-use-rm-id        RepeatMasker ID is the first joining pass (on)
  --merge-max-gap N             max query gap in bp between fragments of one element (100)
  --merge-cons-tolerance N      bp of consensus overlap tolerated between consecutive fragments (20)
  --overlap-max-frac F          a hit overlapping accepted hits by more than this fraction of its length is an alternative (0.20)
  --[no-]ignore-overlapped-hits hits flagged "*" are alternatives, not used for coverage (on)
  --min-primary-element-len N   an element needs this many covered bp to be primary (100)
  --help                        this message

Steps, per insertion:
  1. only hit_kind=candidate hits are used; "*" hits become alternatives (reason: star)
  2. by score, a hit that overlaps accepted hits by more than --overlap-max-frac is an alternative (overlap)
  3. hits of one element are joined: same RepeatMasker ID, then same repeat name + strand + profile with a
     query gap <= --merge-max-gap and consensus coordinates that advance in the direction of the strand
  4. elements whose profile has absorb_inside_span take the hits inside their span that belong to those classes
     (SVA: Simple_repeat, Low_complexity, SINE/Alu); absorbed hits are alternatives (absorbed_by) and add to coverage
  5. elements are ranked by covered bases, then score; rank 1 is primary when it covers --min-primary-element-len bp

elements.tsv, one row per element:
  query, chrom, bed_start, bed_end, q_len   insertion
  elem_id, rank, is_primary                 query#rank; is_primary 1/0
  repeat, group, tier, rm_family, profiles  from the hits
  strand, q_start, q_end, span              span of the member hits in the insertion
  cov_bases, cov                            bases covered by members and absorbed hits; cov = cov_bases / q_len
  n_hits, score_sum
  div_w, div_corr_w                         length-weighted divergence, raw and ONT-corrected
  c_start, c_end, c_len, c_left             consensus span; c_left = consensus bases after c_end (3' side)
  div_limit, young_profiles                 young when tier is young for the profile and div_corr_w <= div_limit ('.' otherwise)
  absorbed_n, absorbed_bases, hit_rows

insertions.tsv, one row per insertion:
  query, chrom, bed_start, bed_end, q_len, n_hits, n_candidate_hits, n_elements
  elem_frac       fraction covered by elements (members and absorbed hits)
  tandem_frac     Simple_repeat + Low_complexity;  satellite_frac Satellite*;  unknown_frac Unknown
                  (hits absorbed into an element do not count; a poly-A tail annotated as (A)n counts as tandem)
USAGE
}

my %o = (
    'merge-use-rm-id'         => 1,
    'merge-max-gap'           => 100,
    'merge-cons-tolerance'    => 20,
    'overlap-max-frac'        => 0.20,
    'ignore-overlapped-hits'  => 1,
    'min-primary-element-len' => 100,
);
GetOptions(\%o, qw(profiles=s insertions-out=s out=s alternatives-out=s merge-use-rm-id! merge-max-gap=i
                   merge-cons-tolerance=i overlap-max-frac=f ignore-overlapped-hits! min-primary-element-len=i help))
    or die usage();
if ($o{help}) { print usage(); exit 0 }
die usage() unless defined $o{profiles} && defined $o{'insertions-out'} && @ARGV == 1;

# ---- absorb rules from class_profiles.tsv
my %absorb;                                       # profile id => [class prefixes]
my ($ph, $prows) = read_tsv($o{profiles});
for my $r (@$prows) {
    next if ($r->{enabled} // 'yes') eq 'no';
    my $p = parse_params($r->{rule_params});
    $absorb{ $r->{profile} } = [ grep { length } split /,/, $p->{absorb_inside_span} ] if defined $p->{absorb_inside_span};
}

# ---- hits
my ($hh, $hits) = read_tsv($ARGV[0]);
my %have = map { $_ => 1 } @$hh;
my @need = qw(hit_row query chrom bed_start bed_end q_start q_end q_len strand score div rm_family c_start c_end
              c_len rm_id overlapped repeat name_norm group tier hit_kind profiles div_corr div_limit tier_young_profiles);
my @missing = grep { !$have{$_} } @need;
die "$ARGV[0]: missing column(s): " . join(', ', @missing) . "\n(run annotate_hits.pl again, it now writes div_limit and tier_young_profiles)\n"
    if @missing;

my (@order, %by_q);
for my $h (@$hits) {
    push @order, $h->{query} unless $by_q{ $h->{query} };
    push @{ $by_q{ $h->{query} } }, $h;
}

# ---- helpers

sub hit_len { return $_[0]{q_end} - $_[0]{q_start} + 1 }

sub cons_ok {                                     # do consensus coordinates advance in the direction of the strand?
    my ($prev, $next) = @_;
    my $tol = $o{'merge-cons-tolerance'};
    return $next->{strand} eq '+'
        ? $next->{c_start} >= $prev->{c_end} - $tol
        : $next->{c_end}   <= $prev->{c_start} + $tol;
}

sub by_q { $a->{q_start} <=> $b->{q_start} || $a->{hit_row} <=> $b->{hit_row} }

sub element_stats {
    my ($e) = @_;
    my @h = @{ $e->{hits} };
    my ($qs, $qe, $cs, $ce);
    my ($wsum, $dsum, $csum, $score) = (0, 0, 0, 0);
    for my $h (@h) {
        $qs = $h->{q_start} if !defined $qs || $h->{q_start} < $qs;
        $qe = $h->{q_end}   if !defined $qe || $h->{q_end}   > $qe;
        $cs = $h->{c_start} if !defined $cs || $h->{c_start} < $cs;
        $ce = $h->{c_end}   if !defined $ce || $h->{c_end}   > $ce;
        my $len = hit_len($h);
        $wsum += $len;
        $dsum += $h->{div} * $len;
        $csum += $h->{div_corr} * $len;
        $score += $h->{score};
    }
    my $f = $h[0];
    my $divc = $csum / $wsum;
    my $young = ($f->{div_limit} ne '.' && $f->{tier_young_profiles} ne '.' && $divc <= $f->{div_limit})
              ? $f->{tier_young_profiles} : '.';
    my %repeats;
    my @names = grep { !$repeats{$_}++ } map { $_->{repeat} } @h;
    my @absorbed = @{ $e->{absorbed} };
    $e->{s} = {
        repeat => join(',', @names), group => $f->{group}, tier => $f->{tier}, rm_family => $f->{rm_family},
        profiles => $f->{profiles}, strand => $f->{strand}, q_start => $qs, q_end => $qe, span => $qe - $qs + 1,
        cov_bases => intervals_len(map { [$_->{q_start}, $_->{q_end}] } @h, @absorbed),
        n_hits => scalar(@h), score_sum => $score, div_w => $dsum / $wsum, div_corr_w => $divc,
        c_start => $cs, c_end => $ce, c_len => $f->{c_len}, c_left => $f->{c_len} - $ce,
        div_limit => $f->{div_limit}, young => $young, absorbed_n => scalar(@absorbed),
        absorbed_bases => intervals_len(map { [$_->{q_start}, $_->{q_end}] } @absorbed),
        hit_rows => join(',', map { $_->{hit_row} } @h),
    };
}

# ---- outputs
my $out = open_out($o{out});
print {$out} tsv_row(qw(query chrom bed_start bed_end q_len elem_id rank is_primary repeat group tier rm_family
                        profiles strand q_start q_end span cov_bases cov n_hits score_sum div_w div_corr_w
                        c_start c_end c_len c_left div_limit young_profiles absorbed_n absorbed_bases hit_rows));
my $ins = open_out($o{'insertions-out'});
print {$ins} tsv_row(qw(query chrom bed_start bed_end q_len n_hits n_candidate_hits n_elements elem_frac
                        tandem_frac satellite_frac unknown_frac));
my $alt;
if (defined $o{'alternatives-out'}) {
    $alt = open_out($o{'alternatives-out'});
    print {$alt} tsv_row(qw(hit_row query repeat rm_family q_start q_end score reason));
}

my (%alt_n, $n_elem, $n_primary, $n_primary_young, $n_merged, $n_with_elem, $n_ins, $n_absorbed);

for my $query (@order) {
    my @hits = @{ $by_q{$query} };
    my $f0   = $hits[0];
    my $q_len = $f0->{q_len};
    my @cand = grep { $_->{hit_kind} eq 'candidate' } @hits;

    my @alts;                                      # [hit, reason]
    my @pool;
    for my $h (@cand) {
        if ($o{'ignore-overlapped-hits'} && $h->{overlapped}) { push @alts, [$h, 'star']; next }
        push @pool, $h;
    }

    # 2. overlap resolution, best score first
    my @accepted;
    for my $h (sort { $b->{score} <=> $a->{score} || $a->{hit_row} <=> $b->{hit_row} } @pool) {
        my $ov = 0;
        for my $acc (@accepted) { $ov += overlap_len(@{$h}{qw(q_start q_end)}, @{$acc}{qw(q_start q_end)}) }
        if ($ov / hit_len($h) > $o{'overlap-max-frac'}) { push @alts, [$h, 'overlap']; next }
        push @accepted, $h;
    }

    # 3. units (RepeatMasker ID), then elements (gap + consensus order)
    my @units;
    if ($o{'merge-use-rm-id'}) {
        my %u;
        push @{ $u{ join("\t", $_->{rm_id}, $_->{name_norm}, $_->{strand}) } }, $_ for @accepted;
        @units = values %u;
    }
    else {
        @units = map { [$_] } @accepted;
    }
    @$_ = sort by_q @$_ for @units;
    @units = sort { $a->[0]{q_start} <=> $b->[0]{q_start} || $a->[0]{hit_row} <=> $b->[0]{hit_row} } @units;

    my @elements;
    for my $u (@units) {
        my $first = $u->[0];
        my $key   = join("\t", $first->{name_norm}, $first->{strand}, $first->{profiles});
        my ($best, $best_gap);
        for my $e (@elements) {
            next unless $e->{key} eq $key;
            my $prev = $e->{hits}[-1];
            my $gap  = $first->{q_start} - $prev->{q_end} - 1;
            next if $gap > $o{'merge-max-gap'};
            next unless cons_ok($prev, $first);
            ($best, $best_gap) = ($e, $gap) if !defined $best_gap || $gap < $best_gap;
        }
        if ($best) {
            @{ $best->{hits} } = sort by_q (@{ $best->{hits} }, @$u);
        }
        else {
            push @elements, { key => $key, hits => [@$u], absorbed => [] };
        }
    }
    element_stats($_) for @elements;

    # 4. absorb rules, biggest elements first
    my %absorbed_row;
    for my $e (sort { $b->{s}{cov_bases} <=> $a->{s}{cov_bases} || $b->{s}{score_sum} <=> $a->{s}{score_sum} } @elements) {
        next if $e->{removed};
        my %prefix = map { $_ => 1 } map { @{ $absorb{$_} || [] } } split /,/, $e->{s}{profiles};
        next unless %prefix;
        my %member = map { $_->{hit_row} => 1 } @{ $e->{hits} };
        my ($qs, $qe) = @{ $e->{s} }{qw(q_start q_end)};
        for my $h (@hits) {
            next if $member{ $h->{hit_row} } || $absorbed_row{ $h->{hit_row} };
            next unless grep { $h->{rm_family} eq $_ || index($h->{rm_family}, "$_/") == 0 } keys %prefix;
            next if overlap_len($h->{q_start}, $h->{q_end}, $qs, $qe) / hit_len($h) < 0.8;
            $absorbed_row{ $h->{hit_row} } = 1;
            push @{ $e->{absorbed} }, $h;
            for my $e2 (@elements) {
                next if $e2 == $e || $e2->{removed};
                my $n = @{ $e2->{hits} };
                @{ $e2->{hits} } = grep { $_->{hit_row} != $h->{hit_row} } @{ $e2->{hits} };
                next if @{ $e2->{hits} } == $n;
                if (@{ $e2->{hits} }) { element_stats($e2) } else { $e2->{removed} = 1 }
            }
        }
        element_stats($e);
    }
    @elements = grep { !$_->{removed} } @elements;

    # 5. ranking
    @elements = sort { $b->{s}{cov_bases} <=> $a->{s}{cov_bases} || $b->{s}{score_sum} <=> $a->{s}{score_sum}
                       || $a->{s}{q_start} <=> $b->{s}{q_start} } @elements;
    my $rank = 0;
    for my $e (@elements) {
        $rank++;
        $e->{rank} = $rank;
        $e->{id}   = "$query#$rank";
        $e->{primary} = ($rank == 1 && $e->{s}{cov_bases} >= $o{'min-primary-element-len'}) ? 1 : 0;
    }
    for my $e (@elements) {
        for my $h (@{ $e->{absorbed} }) { push @alts, [$h, "absorbed_by:$e->{id}"] }
    }

    # outputs
    for my $e (@elements) {
        my $s = $e->{s};
        print {$out} tsv_row(
            $query, @{$f0}{qw(chrom bed_start bed_end)}, $q_len, $e->{id}, $e->{rank}, $e->{primary},
            @{$s}{qw(repeat group tier rm_family profiles strand q_start q_end span cov_bases)},
            sprintf('%.3f', $s->{cov_bases} / $q_len), $s->{n_hits}, $s->{score_sum},
            sprintf('%.2f', $s->{div_w}), sprintf('%.2f', $s->{div_corr_w}),
            @{$s}{qw(c_start c_end c_len c_left div_limit)}, $s->{young}, $s->{absorbed_n}, $s->{absorbed_bases},
            $s->{hit_rows},
        );
        $n_elem++;
        $n_merged++ if $s->{n_hits} > 1;
        if ($e->{primary}) {
            $n_primary++;
            $n_primary_young++ if $s->{young} ne '.';
        }
    }
    $n_with_elem++ if @elements;
    $n_ins++;

    my (@tandem, @sat, @unk, @elem_iv);
    for my $h (@hits) {
        next if $absorbed_row{ $h->{hit_row} };
        my $iv = [$h->{q_start}, $h->{q_end}];
        if    ($h->{rm_family} =~ /^(?:Simple_repeat|Low_complexity)/) { push @tandem, $iv }
        elsif ($h->{rm_family} =~ /^Satellite/)                        { push @sat, $iv }
        elsif ($h->{rm_family} =~ /^Unknown/)                          { push @unk, $iv }
    }
    for my $e (@elements) { push @elem_iv, map { [$_->{q_start}, $_->{q_end}] } @{ $e->{hits} }, @{ $e->{absorbed} } }
    my $frac = sub { my $x = intervals_len(@_) / $q_len; sprintf('%.3f', $x > 1 ? 1 : $x) };
    print {$ins} tsv_row($query, @{$f0}{qw(chrom bed_start bed_end)}, $q_len, scalar(@hits), scalar(@cand),
                         scalar(@elements), $frac->(@elem_iv), $frac->(@tandem), $frac->(@sat), $frac->(@unk));

    for my $alt_hit (@alts) {
        my ($h, $reason) = @$alt_hit;
        (my $short = $reason) =~ s/:.*//;
        next if $short ne 'absorbed_by' && $absorbed_row{ $h->{hit_row} };    # report an absorbed hit once
        $alt_n{$short}++;
        print {$alt} tsv_row(@{$h}{qw(hit_row query repeat rm_family q_start q_end score)}, $reason) if $alt;
    }
    $n_absorbed += scalar(keys %absorbed_row);
}
close $out if defined $o{out};
close $ins;
close $alt if $alt;

warn sprintf("merge_hits: %d insertions, %d with at least one element, %d elements (%d joined from several hits)\n",
             $n_ins || 0, $n_with_elem || 0, $n_elem || 0, $n_merged || 0);
warn sprintf("merge_hits: %d primary elements, %d of them young\n", $n_primary || 0, $n_primary_young || 0);
warn sprintf("merge_hits: alternatives: %s\n", join(', ', map { "$_ $alt_n{$_}" } sort keys %alt_n)) if %alt_n;
