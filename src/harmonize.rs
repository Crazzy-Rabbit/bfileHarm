use crate::plink::BimRecord;
use crate::plink::{create_output, scan_bim};

use std::collections::{HashMap, HashSet};
use std::io::{self, Write};
use std::path::Path;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct Match {
    pub source_idx: u32,
    pub target_idx: u32,
    pub flip_target: bool,
}

#[derive(Debug, Default)]
pub struct HarmonizeStats {
    pub source_snps: usize,
    pub target_snps: usize,
    pub duplicate_ids_source: usize,
    pub duplicate_ids_target: usize,
    pub duplicate_records_source: usize,
    pub duplicate_records_target: usize,
    pub common_ids: usize,
    pub chr_mismatch: usize,
    pub allele_mismatch: usize,
    pub target_flips: usize,
    pub retained: usize,
}

#[derive(Debug)]
pub struct HarmonizationPlan {
    pub matches: Vec<Match>,
    pub stats: HarmonizeStats,
    target_by_snp: HashMap<String, TargetRecord>,
}

#[derive(Clone, Debug)]
struct TargetRecord {
    idx: u32,
    chr: String,
    cm: String,
    bp: String,
    a1: String,
    a2: String,
    duplicate: bool,
}

#[derive(Debug, Default)]
struct DuplicateScan {
    total: usize,
    duplicate_ids: usize,
    duplicate_records: usize,
    duplicate_set: HashSet<String>,
}

#[derive(Debug, Default)]
struct TargetBuild {
    total: usize,
    duplicate_ids: usize,
    duplicate_records: usize,
    target_by_snp: HashMap<String, TargetRecord>,
}

pub fn plan_harmonization_from_bim_files(
    source_bim: &Path,
    target_bim: &Path,
    strict_duplicates: bool,
) -> Result<HarmonizationPlan, String> {
    let source_dups = scan_source_duplicates(source_bim)
        .map_err(|e| format!("failed to scan bfile1 BIM duplicates: {e}"))?;
    let target_build =
        build_target_map(target_bim).map_err(|e| format!("failed to index bfile2 BIM: {e}"))?;

    if strict_duplicates && (source_dups.duplicate_ids > 0 || target_build.duplicate_ids > 0) {
        return Err(format!(
            "duplicate SNP IDs found: {} IDs in bfile1, {} IDs in bfile2",
            source_dups.duplicate_ids, target_build.duplicate_ids
        ));
    }

    let mut stats = HarmonizeStats {
        source_snps: source_dups.total,
        target_snps: target_build.total,
        duplicate_ids_source: source_dups.duplicate_ids,
        duplicate_ids_target: target_build.duplicate_ids,
        duplicate_records_source: source_dups.duplicate_records,
        duplicate_records_target: target_build.duplicate_records,
        ..HarmonizeStats::default()
    };
    let mut matches = Vec::new();

    scan_bim(source_bim, |source_idx, source| {
        if source_dups.duplicate_set.contains(&source.snp) {
            return Ok(());
        }
        let Some(target) = target_build.target_by_snp.get(&source.snp) else {
            return Ok(());
        };
        if target.duplicate {
            return Ok(());
        }
        stats.common_ids += 1;

        if !same_chr(&source.chr, &target.chr) {
            stats.chr_mismatch += 1;
            return Ok(());
        }

        let flip_target = if same_alleles_to_target(&source, target) {
            false
        } else if swapped_alleles_to_target(&source, target) {
            true
        } else {
            stats.allele_mismatch += 1;
            return Ok(());
        };

        if flip_target {
            stats.target_flips += 1;
        }
        matches.push(Match {
            source_idx: u32::try_from(source_idx).map_err(|_| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    "bfile1 BIM has more than 4,294,967,295 SNPs",
                )
            })?,
            target_idx: target.idx,
            flip_target,
        });
        Ok(())
    })
    .map_err(|e| format!("failed to plan harmonization from bfile1 BIM: {e}"))?;

    stats.retained = matches.len();
    if matches.is_empty() {
        return Err("no SNPs left after common-SNP and allele harmonization".to_string());
    }

    Ok(HarmonizationPlan {
        matches,
        stats,
        target_by_snp: target_build.target_by_snp,
    })
}

pub fn write_harmonized_bims(
    source_bim: &Path,
    out1_bim: &Path,
    out2_bim: &Path,
    plan: &HarmonizationPlan,
    force: bool,
) -> io::Result<()> {
    let mut out1 = create_output(out1_bim, force)?;
    let mut out2 = create_output(out2_bim, force)?;
    let mut next_match = 0usize;

    scan_bim(source_bim, |source_idx, source| {
        if next_match >= plan.matches.len() {
            return Ok(());
        }
        let m = plan.matches[next_match];
        if usize::try_from(m.source_idx).unwrap_or(usize::MAX) != source_idx {
            return Ok(());
        }
        let target = plan.target_by_snp.get(&source.snp).ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "internal error: missing target BIM record for {}",
                    source.snp
                ),
            )
        })?;
        if target.idx != m.target_idx {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("internal error: target index changed for {}", source.snp),
            ));
        }

        writeln!(
            out1,
            "{}\t{}\t{}\t{}\t{}\t{}",
            source.chr, source.snp, source.cm, source.bp, source.a1, source.a2
        )?;
        writeln!(
            out2,
            "{}\t{}\t{}\t{}\t{}\t{}",
            target.chr, source.snp, target.cm, target.bp, source.a1, source.a2
        )?;
        next_match += 1;
        Ok(())
    })?;

    if next_match != plan.matches.len() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "internal error: not all retained SNPs were written to BIM outputs",
        ));
    }
    out1.flush()?;
    out2.flush()
}

fn scan_source_duplicates(path: &Path) -> io::Result<DuplicateScan> {
    let mut scan = DuplicateScan::default();
    let mut seen = HashSet::new();

    scan_bim(path, |_idx, rec| {
        scan.total += 1;
        if !seen.insert(rec.snp.clone()) {
            if scan.duplicate_set.insert(rec.snp) {
                scan.duplicate_ids += 1;
                scan.duplicate_records += 2;
            } else {
                scan.duplicate_records += 1;
            }
        }
        Ok(())
    })?;

    Ok(scan)
}

fn build_target_map(path: &Path) -> io::Result<TargetBuild> {
    let mut build = TargetBuild::default();

    scan_bim(path, |idx, rec| {
        build.total += 1;
        let idx = u32::try_from(idx).map_err(|_| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                "bfile2 BIM has more than 4,294,967,295 SNPs",
            )
        })?;

        if let Some(existing) = build.target_by_snp.get_mut(&rec.snp) {
            if !existing.duplicate {
                existing.duplicate = true;
                build.duplicate_ids += 1;
                build.duplicate_records += 2;
            } else {
                build.duplicate_records += 1;
            }
            return Ok(());
        }

        build.target_by_snp.insert(
            rec.snp,
            TargetRecord {
                idx,
                chr: rec.chr,
                cm: rec.cm,
                bp: rec.bp,
                a1: rec.a1,
                a2: rec.a2,
                duplicate: false,
            },
        );
        Ok(())
    })?;

    Ok(build)
}

#[cfg(test)]
fn plan_harmonization(
    source: &[BimRecord],
    target: &[BimRecord],
    strict_duplicates: bool,
) -> Result<(Vec<Match>, HarmonizeStats), String> {
    let source_dups = duplicate_ids(source);
    let target_dups = duplicate_ids(target);

    if strict_duplicates && (!source_dups.is_empty() || !target_dups.is_empty()) {
        return Err(format!(
            "duplicate SNP IDs found: {} IDs in bfile1, {} IDs in bfile2",
            source_dups.len(),
            target_dups.len()
        ));
    }

    let target_map = unique_snp_map(target, &target_dups);
    let mut matches = Vec::new();
    let mut stats = HarmonizeStats {
        source_snps: source.len(),
        target_snps: target.len(),
        duplicate_ids_source: source_dups.len(),
        duplicate_ids_target: target_dups.len(),
        duplicate_records_source: count_duplicate_records(source, &source_dups),
        duplicate_records_target: count_duplicate_records(target, &target_dups),
        ..HarmonizeStats::default()
    };

    for (source_idx, s) in source.iter().enumerate() {
        if source_dups.contains(&s.snp) {
            continue;
        }
        let Some(&target_idx) = target_map.get(&s.snp) else {
            continue;
        };
        stats.common_ids += 1;

        let t = &target[target_idx];
        if !same_chr(&s.chr, &t.chr) {
            stats.chr_mismatch += 1;
            continue;
        }

        let flip_target = if same_alleles(s, t) {
            false
        } else if swapped_alleles(s, t) {
            true
        } else {
            stats.allele_mismatch += 1;
            continue;
        };

        if flip_target {
            stats.target_flips += 1;
        }
        matches.push(Match {
            source_idx: source_idx.try_into().unwrap(),
            target_idx: target_idx.try_into().unwrap(),
            flip_target,
        });
    }

    stats.retained = matches.len();
    if matches.is_empty() {
        return Err("no SNPs left after common-SNP and allele harmonization".to_string());
    }
    Ok((matches, stats))
}

#[cfg(test)]
fn duplicate_ids(records: &[BimRecord]) -> HashSet<String> {
    let mut seen = HashSet::new();
    let mut duplicates = HashSet::new();
    for rec in records {
        if !seen.insert(rec.snp.clone()) {
            duplicates.insert(rec.snp.clone());
        }
    }
    duplicates
}

#[cfg(test)]
fn unique_snp_map(records: &[BimRecord], duplicates: &HashSet<String>) -> HashMap<String, usize> {
    let mut map = HashMap::new();
    for (idx, rec) in records.iter().enumerate() {
        if !duplicates.contains(&rec.snp) {
            map.insert(rec.snp.clone(), idx);
        }
    }
    map
}

#[cfg(test)]
fn count_duplicate_records(records: &[BimRecord], duplicates: &HashSet<String>) -> usize {
    records
        .iter()
        .filter(|rec| duplicates.contains(&rec.snp))
        .count()
}

fn same_chr(a: &str, b: &str) -> bool {
    normalize_chr(a) == normalize_chr(b)
}

fn normalize_chr(chr: &str) -> String {
    chr.trim()
        .strip_prefix("chr")
        .or_else(|| chr.trim().strip_prefix("CHR"))
        .unwrap_or(chr.trim())
        .to_ascii_uppercase()
}

#[cfg(test)]
fn same_alleles(source: &BimRecord, target: &BimRecord) -> bool {
    source.a1 == target.a1 && source.a2 == target.a2
}

#[cfg(test)]
fn swapped_alleles(source: &BimRecord, target: &BimRecord) -> bool {
    source.a1 == target.a2 && source.a2 == target.a1
}

fn same_alleles_to_target(source: &BimRecord, target: &TargetRecord) -> bool {
    source.a1 == target.a1 && source.a2 == target.a2
}

fn swapped_alleles_to_target(source: &BimRecord, target: &TargetRecord) -> bool {
    source.a1 == target.a2 && source.a2 == target.a1
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rec(chr: &str, snp: &str, a1: &str, a2: &str) -> BimRecord {
        BimRecord {
            chr: chr.to_string(),
            snp: snp.to_string(),
            cm: "0".to_string(),
            bp: "1".to_string(),
            a1: a1.to_string(),
            a2: a2.to_string(),
        }
    }

    #[test]
    fn plan_keeps_source_order_and_marks_target_flip() {
        let source = vec![rec("1", "rs1", "A", "C"), rec("1", "rs2", "G", "T")];
        let target = vec![rec("1", "rs2", "T", "G"), rec("1", "rs1", "A", "C")];

        let (matches, stats) = plan_harmonization(&source, &target, false).unwrap();

        assert_eq!(matches.len(), 2);
        assert_eq!(
            matches[0],
            Match {
                source_idx: 0,
                target_idx: 1,
                flip_target: false
            }
        );
        assert_eq!(
            matches[1],
            Match {
                source_idx: 1,
                target_idx: 0,
                flip_target: true
            }
        );
        assert_eq!(stats.target_flips, 1);
    }

    #[test]
    fn plan_skips_allele_mismatch() {
        let source = vec![rec("1", "rs1", "A", "C"), rec("1", "rs2", "G", "T")];
        let target = vec![rec("1", "rs1", "A", "G"), rec("1", "rs2", "G", "T")];

        let (matches, stats) = plan_harmonization(&source, &target, false).unwrap();

        assert_eq!(matches.len(), 1);
        assert_eq!(matches[0].source_idx, 1);
        assert_eq!(stats.allele_mismatch, 1);
    }
}
