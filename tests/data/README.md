# Test data

Minimal dataset used by the `test` profile (`conf/test.config`) and by the
nf-test suite. It is deliberately tiny so that the pipeline runs in seconds.

## Files

| File | Description |
|------|-------------|
| `samplesheet.csv` | Samplesheet used by `-profile test` (columns: `sample,vcf,tbi`) |
| `test.vcf.gz` | Subset of a Sniffles2 multi-record VCF (insertions and other SV types) |
| `test.vcf.gz.tbi` | Tabix index of `test.vcf.gz` |

## Provenance

- **Source:** private long-read dataset (origin not disclosed). The sample
  name in the VCF was replaced by `toy`.
- **Organism / reference build:** <human, GRCh38 | T2T-CHM13 | ...>
- **Sequencing platform:** <ONT | PacBio HiFi>
- **Caller:** Sniffles <version>, run with: `<sniffles command / key options>`
- **Insertion sequence format:** <inline in ALT | symbolic <INS> with INFO/SEQ>

## How the subset was generated

```bash
# 1. Extract one region (contains a mix of passing and failing INS records)
bcftools view -r chr22:20000000-25000000 <full_sniffles.vcf.gz> -Oz -o subset.vcf.gz

# 2. Dump the header and drop the unwanted lines
bcftools view -h subset.vcf.gz | grep -v -E 'sniffles|command|VEP|bcftools' > header.txt

# 3. Replace the header and rename the sample in the same step
echo "<original_sample_name> toy" > rename.txt
bcftools reheader -h header.txt -s rename.txt subset.vcf.gz \
  -o test.vcf.gz

# 4. Re-index and replace the old file
mv test.vcf.gz tests/data/test.vcf.gz
tabix -p vcf tests/data/test.vcf.gz
```

## Content summary

- Records: <N> total, <N> `SVTYPE=INS`, <N> with `FILTER=PASS`
- Expected to be kept by the default filters: <N> (see `params.min_support`, etc.)
- Contains known retroelement insertions: <yes/no, which type if known>

## Regenerating or extending

Repeat the commands above with a different region. If the region changes,
update the counts in this file and refresh the nf-test snapshot:

```bash
nf-test test --update-snapshot
```

## License / data sharing

<State the conditions under which this subset may be redistributed.>
