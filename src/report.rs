use crate::harmonize::HarmonizeStats;
use crate::Config;

pub fn format_report(config: &Config, stats: &HarmonizeStats) -> String {
    format!(
        "\
bfileHarm report
  bfile1: {bfile1}
  bfile2: {bfile2}
  out1:   {out1}
  out2:   {out2}

Input SNPs
  bfile1: {source_snps}
  bfile2: {target_snps}

Filtering
  duplicated SNP IDs in bfile1: {dup_ids_source} IDs / {dup_records_source} records
  duplicated SNP IDs in bfile2: {dup_ids_target} IDs / {dup_records_target} records
  common unique SNP IDs:        {common_ids}
  chromosome mismatches:        {chr_mismatch}
  allele mismatches:            {allele_mismatch}
  bfile2 allele flips:          {target_flips}

Output
  retained harmonized SNPs:     {retained}
",
        bfile1 = config.bfile1,
        bfile2 = config.bfile2,
        out1 = config.out1,
        out2 = config.out2,
        source_snps = stats.source_snps,
        target_snps = stats.target_snps,
        dup_ids_source = stats.duplicate_ids_source,
        dup_records_source = stats.duplicate_records_source,
        dup_ids_target = stats.duplicate_ids_target,
        dup_records_target = stats.duplicate_records_target,
        common_ids = stats.common_ids,
        chr_mismatch = stats.chr_mismatch,
        allele_mismatch = stats.allele_mismatch,
        target_flips = stats.target_flips,
        retained = stats.retained,
    )
}
