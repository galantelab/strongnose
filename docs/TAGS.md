# Strongnose tag reference

Every fixed value ("tag") the scripts write, grouped by the table where you meet it. Taken from the code; when a
script changes a tag, change it here too. `.` always means "none" or "not applicable".

Reading order for a call in `events.tsv`: `confidence`, then `gates_failed` (why it is capped at low),
then the evidence columns (`polya_found`, `tsd_found`, `tsd_capped`), then `flags` and `note`.

## 1. Final classification (`events.tsv`)

### confidence

| Tag | Meaning |
|---|---|
| `high` | young element, no gate failed, every evidence item of the profile present |
| `probable` | young element, no gate failed, some but not all evidence items present |
| `low` | everything else: not young, or a gate failed, or no evidence at all |

A solo LTR (`HERVK_solo_ltr`) can reach `probable` at most, because its TSD (4-6 bp) cannot be detected reliably.

### gates_failed (comma list; any gate caps the confidence at low)

| Tag | Meaning | Set by |
|---|---|---|
| `length` | insertion length outside `len_min`-`len_max` of the profile | class_profiles.tsv |
| `coverage` | the element covers less than `min_cov` of the insertion (for SVA, same-group pieces and 5' pieces are counted) | class_profiles.tsv |
| `three_prime_end_missing` | more than `max_3end_left` consensus bases are missing at the 3' end (`c_left`) | rule_params |
| `strand_switches` | more elements of the profile on the opposite strand than `max_strand_switch` | rule_params |

### flags (comma list; informational, they do not change the confidence)

| Tag | Meaning | What it usually is |
|---|---|---|
| `fragment_no_3end` | the 3' end gate failed | a fragment from the middle of the consensus; for SVA typically a VNTR length change of an existing element |
| `tandem_dup_of_flank` | the TSD run reached the scan limit (`tsd_capped` = 1) | the insertion repeats the neighbouring reference: tandem duplication or length change of an existing repeat |
| `element_after_3end` | another element of 100 bp or more, not a piece of the same group, lies after the 3' end of the primary | 3' transduction, composite insertion or a second event; review by hand |

### note

| Tag | Meaning |
|---|---|
| `ltr_identity=X` | identity of the two LTRs (provirus profiles); `.` for a solo LTR |
| `structure_fail=REASONS` | the LTR structure was found but is not valid; reasons in section 5 |
| `no_ltr_structure:REASON` | a provirus or solo-LTR profile without a usable structure row (`single_ltr`, `no_pair`, `ltr_length`, `no_ltr_pair`) |
| `three_prime_piece` | the event is judged by a young piece of the same group that carries the 3' end, not by the largest piece (SVA split among subfamily names) |
| `rule_pending:RULE` | the rule of the profile has no implementation, or `--structure` was not given; confidence is low |
| `unchecked_evidence:ITEMS` | the profile lists an evidence item the classifier does not know |

### evidence columns

| Column | Values |
|---|---|
| `young` | `1` = young for the profile (tier in `young_tiers` and corrected divergence under the ceiling), `0` otherwise |
| `polya_found`, `tsd_found` | `1` / `0` as detected; `.` when the insertion has no row in that table. Raw detection: the classifier applies more conditions before it counts them (`evidence_ok`) |
| `tsd_capped` | `1` = the run reached `tsd_scan_max`; never counts as a TSD |
| `evidence_ok` / `evidence_total` | evidence items that count / items the profile asks for |
| `other_elements` | `repeat/strand/coverage` of the other elements of the insertion, `;` separated |

## 2. Dropped insertions (`discarded.tsv`, column `reason`)

| Tag | Meaning |
|---|---|
| `no_repeat_hit` | RepeatMasker reported nothing for the insertion (listed from the BED) |
| `no_candidate_element` | hits exist, but none belongs to an enabled profile |
| `tandem_repeat` | Simple_repeat + Low_complexity cover more than `artifact_max_tandem_frac` |
| `satellite` | Satellite covers more than `artifact_max_satellite_frac` |
| `unknown_sequence` | Unknown covers more than `artifact_max_unknown_frac` |
| `no_primary_element` | the largest element is shorter than `min_primary_element_len` |
| `no_enabled_profile` | the primary element belongs only to disabled profiles |
| `below_min_confidence` | classified, but under `min_confidence` |

## 3. Age and hit annotation (`subfamilies.tsv`, annotated hits, `elements.tsv`)

### tier

| Tag | Meaning |
|---|---|
| `active` | subfamily with known recent activity (AluYa5, AluYb8, AluYb9, L1HS, SVA_E/F/F1, LTR5_Hs) |
| `young` | young subfamily with weaker evidence of activity (other AluY, L1PA2-6, SVA_D) |
| `intermediate` | mostly fixed (AluS, L1PA7/8, SVA_A-C, LTR5A/B) |
| `old` | ancient; also the default for names that are not in the table |
| `internal` | internal region of a multi-part element (HERVK-int); its age comes from the LTRs |
| `generic` | the name gives no subfamily (plain `Alu`); judged by divergence, with a stricter ceiling |

### other annotation tags

| Column | Tag | Meaning |
|---|---|---|
| `tier_source` | `table` / `default` | tier read from subfamilies.tsv / name not listed, so `old` |
| `hit_kind` | `candidate` | class matches an enabled profile (and the group agrees) |
| | `artifact` | class in `artifact_classes` (Simple_repeat, Low_complexity, Satellite, Unknown) |
| | `other` | any other repeat |
| `overlapped` | `1` | RepeatMasker `*` flag: the hit is overlapped by a higher-scoring hit |
| `is_primary` | `1` | the element that covers the most bases of the insertion |

### profiles (`class_profiles.tsv`)

`Alu`, `L1`, `SVA`, `HERVK_provirus`, `HERVK_solo_ltr` (disabled by default). Rules: `single_element`, `ltr_int_ltr`,
`solo_ltr`. Evidence items: `polya`, `tsd`, `structure`.

### hits not used by the merge (`--alternatives-out`, column `reason`)

| Tag | Meaning |
|---|---|
| `star` | the hit carries the `*` flag |
| `overlap` | the hit overlaps accepted hits by more than `overlap_max_frac` of its length |
| `absorbed_by:ELEMENT` | the hit was counted as part of that element (SVA components) |

## 4. Poly-A and TSD tables

| Table | Column | Tags |
|---|---|---|
| `polya.tsv` | `tail_side` | `right` (plus-strand element) / `left` (minus-strand element): side of the insertion where the tail is |
| `flanks.tsv` | `status` | `ok` / `chrom_missing` (contig not in the `.fai`) / `edge` (window outside the contig) |
| `tsd.tsv` | `status` | `ok` / `no_flank` (no usable flank window) / `no_sequence` (insertion not in the FASTA) / `too_short` (insertion under 8 bp, or incomplete window) |

`tsd_delta` is the junction shift, in bp, at which the TSD was found (0 = the position reported by Sniffles).

## 5. LTR structure (`ltr_structure.tsv`)

### status

| Tag | Meaning |
|---|---|
| `ok` | provirus: a valid-looking pair was found (check `structure_ok`); solo LTR: a solo candidate |
| `one_ltr` | provirus profile, but the insertion has a single LTR |
| `no_pair` | two or more LTRs, none satisfies the pair conditions |
| `not_solo` | solo-LTR profile, but the insertion has an internal region or a second LTR; ignored by the classifier |

### fail_reason

| Tag | Profile | Meaning |
|---|---|---|
| `ltr_identity` | provirus | identity of the two LTRs under `ltr_min_identity`, or not computed |
| `internal_too_short` | provirus | fewer than `internal_min` bp between the LTRs |
| `internal_orientation` | provirus | internal region on the opposite strand of the LTRs over more than 20% |
| `single_ltr` | provirus | only one LTR |
| `no_pair` | provirus | no pair with the same orientation and enough internal sequence |
| `ltr_length` | both | LTR length outside `ltr_len` (or `len_min`-`len_max` for the solo profile) |
| `internal_present` | solo | internal-region elements cover more than 20% of the insertion |
| `multiple_ltr` | solo | more than one LTR |

## 6. Manual review (`review.tsv`)

### category

| Group | Tag | Events |
|---|---|---|
| validation | `high_alu`, `high_l1`, `high_sva` | calls at high (SVA: all of them) |
| validation | `probable_alu`, `probable_l1`, `probable_sva` | calls at probable |
| validation | `hervk_ok` | HERVK calls above low |
| rule question | `sva_element_after_3end` | SVA with `element_after_3end` |
| rule question | `sva_no3end_with_evidence`, `l1_no3end_with_evidence` | 3' gate failed although poly-A and TSD were found |
| rule question | `gate_rescue_alu`, `gate_rescue_l1` | young, a gate failed, poly-A and TSD found |
| rule question | `tandem_dup_young` | young with `tandem_dup_of_flank` |
| calibration | `age_near_ceiling` | high Alus with the highest corrected divergence |
| calibration | `old_alu_with_evidence` | AluS or older with poly-A and TSD |
| calibration | `short_tsd_high` | high with a TSD of 9 bp or less |
| false negatives | `dropped_REASON` | samples of the dropped insertions (reasons of section 2) |

### verdict (filled by the reviewers)

`true_insertion`, `false_positive`, `uncertain`.
