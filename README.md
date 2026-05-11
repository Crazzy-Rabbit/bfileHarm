# bfileHarm

`bfileHarm` harmonizes two PLINK binary bfile datasets (`.bed/.bim/.fam`).

It keeps common SNPs, aligns A1/A2 allele order to `bfile1`, flips `bfile2`
genotypes when necessary, and writes two harmonized bfile outputs with identical
SNP order.

## Features

- Match common variants by SNP ID from BIM column 2.
- Keep only SNPs with matching chromosome and compatible allele pairs.
- Align both output BIM files to the A1/A2 order of `bfile1`.
- Flip `bfile2.bed` genotype encoding for swapped alleles.
- Stream BED processing block by block; genotype data are not loaded into memory.
- Use a memory-conscious BIM planning path for large WGS bfiles.
- Optionally write two BED outputs concurrently with `--bed-threads 2`.

## Build

Build with Rust/Cargo:

```bash
cargo build --release
```

The executable is:

```bash
target/release/bfileHarm
```

## Input

Both inputs must be PLINK binary bfile prefixes:

```text
prefix.bed
prefix.bim
prefix.fam
```

Pass prefixes without file extensions:

```bash
--bfile1 1kgp_EUR
--bfile2 1kgp_EAS
```

The `.bed` files must be standard PLINK SNP-major binary BED files.

## Output

Use explicit output prefixes:

```bash
bfileHarm \
  --bfile1 1kgp_EUR \
  --bfile2 1kgp_EAS \
  --out1 1kgp_EUR.harm \
  --out2 1kgp_EAS.harm
```

This writes:

```text
1kgp_EUR.harm.bed
1kgp_EUR.harm.bim
1kgp_EUR.harm.fam
1kgp_EAS.harm.bed
1kgp_EAS.harm.bim
1kgp_EAS.harm.fam
```

Or use a shared output prefix:

```bash
bfileHarm --bfile1 1kgp_EUR --bfile2 1kgp_EAS --out-prefix 1kgp_EUR_EAS.harm
```

This writes:

```text
1kgp_EUR_EAS.harm.bfile1.bed/.bim/.fam
1kgp_EUR_EAS.harm.bfile2.bed/.bim/.fam
```

## Usage

```text
bfileHarm --bfile1 <prefix1> --bfile2 <prefix2> --out1 <out1> --out2 <out2>
bfileHarm --bfile1 <prefix1> --bfile2 <prefix2> --out-prefix <out>
bfileHarm <prefix1> <prefix2> <out-prefix>
bfileHarm <prefix1> <prefix2> <out1> <out2>
```

Options:

| Option | Description |
|---|---|
| `--bfile1`, `--source` | First/source bfile prefix. Output allele order follows this file. |
| `--bfile2`, `--target` | Second/target bfile prefix. Genotypes are flipped when A1/A2 are swapped. |
| `--out1` | Output prefix for harmonized `bfile1`. |
| `--out2` | Output prefix for harmonized `bfile2`. |
| `--out-prefix`, `--out` | Shared output prefix. Writes `<prefix>.bfile1.*` and `<prefix>.bfile2.*`. |
| `--report` | Optional text report path. |
| `--bed-threads <1\|2>` | BED writer threads. Default: `1`. Use `2` to write both BED outputs concurrently. |
| `--strict-duplicates` | Stop on duplicated SNP IDs instead of skipping duplicated IDs. |
| `--force` | Overwrite existing output files. |
| `-h`, `--help` | Show help text. |

## Harmonization Rules

For each SNP in `bfile1.bim`, `bfileHarm` checks whether the same SNP ID exists
in `bfile2.bim`.

A SNP is retained only if:

- SNP ID is unique in both BIM files, unless `--strict-duplicates` is used, in
  which case duplicates stop the run.
- Chromosome is the same in both BIM files.
- Allele pair is compatible after uppercasing:
  - same order: `bfile1 A1/A2 == bfile2 A1/A2`
  - swapped order: `bfile1 A1/A2 == bfile2 A2/A1`

When alleles are swapped in `bfile2`:

- `bfile2.bed` genotype codes are flipped as `0 <-> 2`.
- Heterozygous and missing genotypes are unchanged.
- Output `bfile2.bim` receives `bfile1` A1/A2 allele order.

Both output bfiles are written in the retained SNP order of `bfile1`.
`.fam` files are copied unchanged.

## Performance

The implementation is designed for large reference bfiles such as WGS 1KGP
subsets.

### Memory

Peak memory is mainly driven by:

```text
bfile2 BIM SNP index
retained SNP match list
BED action map
64 MiB BED buffer per active BED writer
```

BED genotype data are streamed and are not kept in memory.

For about 500 samples:

```text
BED bytes per SNP = ceil(500 / 4) = 125 bytes
BED action map    = 1 byte per input SNP
match list        ≈ 12 bytes per retained SNP
```

Approximate peak RAM:

| Variant count | Expected RAM |
|---:|---:|
| 10M | 3-7 GB |
| 30M | 9-20 GB |
| 80M | 22-50 GB |

These are rough estimates. Very long SNP IDs or highly fragmented BIM strings
can increase memory use.

### Disk

For 500 samples and `K` retained SNPs:

```text
two output BED files ≈ 2 * K * 125 bytes
two output BIM files ≈ 2 * K * 40-100 bytes
```

For 80M retained SNPs, output size is roughly:

```text
BED: 20 GB
BIM: 6-16 GB
total: 26-36 GB
```

### Runtime

BED processing is fastest when both inputs have similar SNP order.

- `bfile1.bed` is always streamed in input order.
- `bfile2.bed` is streamed when its retained SNP indices are also increasing.
- If `bfile2` order differs strongly from `bfile1`, `bfile2.bed` falls back to
  ordered random-access extraction to preserve identical output SNP order.

`--bed-threads 2` may improve throughput on NVMe SSDs or separate disks. On HDDs
or saturated shared storage, it may be slower than the default `--bed-threads 1`.

## Recommended Workflow

For full WGS reference panels:

1. Keep both input bfiles sorted in the same genomic/SNP order whenever possible.
2. Prefer NVMe SSD storage for `.bed` I/O.
3. Use `--bed-threads 2` only when storage bandwidth is sufficient.
4. Consider running by chromosome for very large panels or memory-constrained
   machines.

Example:

```bash
bfileHarm \
  --bfile1 1kgp_EUR_chr1 \
  --bfile2 1kgp_EAS_chr1 \
  --out-prefix 1kgp_EUR_EAS_chr1.harm \
  --bed-threads 2 \
  --force
```

## Notes

- Strand-complement matching is not currently performed. Alleles must match
  exactly as the same pair or swapped pair after uppercasing.
- Multi-allelic representations with duplicated SNP IDs are skipped by default.
- Input prefixes should not include `.bed`, `.bim`, or `.fam`.
