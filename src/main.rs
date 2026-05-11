mod harmonize;
mod plink;
mod report;

use std::env;
use std::path::PathBuf;
use std::process;
use std::thread;

use harmonize::{plan_harmonization_from_bim_files, write_harmonized_bims};
use plink::{
    ensure_outputs_available, make_bed_action_map, path_with_ext, read_fam, write_bed_action_map,
    write_bed_subset, write_fam, write_text,
};
use report::format_report;

const HELP: &str = "\
bfileHarm

Harmonize two PLINK binary bfile datasets by common SNP ID and A1/A2 allele order.

Usage:
  bfileHarm --bfile1 <prefix1> --bfile2 <prefix2> --out1 <prefix1_out> --out2 <prefix2_out>
  bfileHarm --bfile1 <prefix1> --bfile2 <prefix2> --out-prefix <prefix>
  bfileHarm <prefix1> <prefix2> <out-prefix>
  bfileHarm <prefix1> <prefix2> <out1-prefix> <out2-prefix>

Options:
  --bfile1, --source <prefix>     First/source bfile prefix.
  --bfile2, --target <prefix>     Second/target bfile prefix.
  --out1 <prefix>                 Output prefix for harmonized bfile1.
  --out2 <prefix>                 Output prefix for harmonized bfile2.
  --out-prefix, --out <prefix>    Write <prefix>.bfile1.* and <prefix>.bfile2.*.
  --report <path>                 Optional text report path.
  --bed-threads <1|2>             BED writer threads. Default: 1.
  --strict-duplicates             Error on duplicate SNP IDs instead of skipping them.
  --force                         Overwrite existing output files.
  -h, --help                      Show this help.

Input prefixes should not include .bed/.bim/.fam.
";

#[derive(Debug)]
pub struct Config {
    pub bfile1: String,
    pub bfile2: String,
    pub out1: String,
    pub out2: String,
    pub report: Option<PathBuf>,
    pub bed_threads: usize,
    pub strict_duplicates: bool,
    pub force: bool,
}

fn main() {
    if let Err(err) = run() {
        eprintln!("Error: {err}");
        process::exit(1);
    }
}

fn run() -> Result<(), String> {
    let Some(config) = parse_args()? else {
        println!("{HELP}");
        return Ok(());
    };

    run_with_config(&config)
}

fn run_with_config(config: &Config) -> Result<(), String> {
    let fam1 = read_fam(&path_with_ext(&config.bfile1, "fam"))
        .map_err(|e| format!("failed to read bfile1 FAM: {e}"))?;
    let fam2 = read_fam(&path_with_ext(&config.bfile2, "fam"))
        .map_err(|e| format!("failed to read bfile2 FAM: {e}"))?;

    let bfile1_bim = path_with_ext(&config.bfile1, "bim");
    let bfile2_bim = path_with_ext(&config.bfile2, "bim");
    let plan =
        plan_harmonization_from_bim_files(&bfile1_bim, &bfile2_bim, config.strict_duplicates)?;

    ensure_outputs_available(
        &config.out1,
        &config.out2,
        config.report.as_deref(),
        config.force,
    )
    .map_err(|e| format!("output path check failed: {e}"))?;

    write_fam(
        &path_with_ext(&config.out1, "fam"),
        &fam1.content,
        config.force,
    )
    .map_err(|e| format!("failed to write out1 FAM: {e}"))?;
    write_fam(
        &path_with_ext(&config.out2, "fam"),
        &fam2.content,
        config.force,
    )
    .map_err(|e| format!("failed to write out2 FAM: {e}"))?;

    write_harmonized_bims(
        &bfile1_bim,
        &path_with_ext(&config.out1, "bim"),
        &path_with_ext(&config.out2, "bim"),
        &plan,
        config.force,
    )
    .map_err(|e| format!("failed to write harmonized BIM files: {e}"))?;

    if config.bed_threads >= 2 {
        thread::scope(|scope| {
            let out1_job = scope.spawn(|| {
                write_source_bed(config, fam1.n_samples, &plan)
                    .map_err(|e| format!("failed to write out1 BED: {e}"))
            });
            let out2_job = scope.spawn(|| {
                write_target_bed(config, fam2.n_samples, &plan)
                    .map_err(|e| format!("failed to write out2 BED: {e}"))
            });

            out1_job
                .join()
                .map_err(|_| "out1 BED writer thread panicked".to_string())??;
            out2_job
                .join()
                .map_err(|_| "out2 BED writer thread panicked".to_string())?
        })?;
    } else {
        write_source_bed(config, fam1.n_samples, &plan)
            .map_err(|e| format!("failed to write out1 BED: {e}"))?;
        write_target_bed(config, fam2.n_samples, &plan)
            .map_err(|e| format!("failed to write out2 BED: {e}"))?;
    }

    let report = format_report(config, &plan.stats);
    print!("{report}");
    if let Some(report_path) = &config.report {
        write_text(report_path, &report, config.force)
            .map_err(|e| format!("failed to write report: {e}"))?;
    }

    Ok(())
}

fn write_source_bed(
    config: &Config,
    n_samples: usize,
    plan: &harmonize::HarmonizationPlan,
) -> std::io::Result<()> {
    let actions = make_bed_action_map(
        plan.stats.source_snps,
        plan.matches.iter().map(|m| (m.source_idx as usize, false)),
    )?;
    write_bed_action_map(
        &path_with_ext(&config.bfile1, "bed"),
        &path_with_ext(&config.out1, "bed"),
        n_samples,
        &actions,
        config.force,
    )
}

fn write_target_bed(
    config: &Config,
    n_samples: usize,
    plan: &harmonize::HarmonizationPlan,
) -> std::io::Result<()> {
    if target_bed_can_stream(plan) {
        let actions = make_bed_action_map(
            plan.stats.target_snps,
            plan.matches
                .iter()
                .map(|m| (m.target_idx as usize, m.flip_target)),
        )?;
        write_bed_action_map(
            &path_with_ext(&config.bfile2, "bed"),
            &path_with_ext(&config.out2, "bed"),
            n_samples,
            &actions,
            config.force,
        )
    } else {
        write_bed_subset(
            &path_with_ext(&config.bfile2, "bed"),
            &path_with_ext(&config.out2, "bed"),
            n_samples,
            plan.stats.target_snps,
            plan.matches
                .iter()
                .map(|m| (m.target_idx as usize, m.flip_target)),
            config.force,
        )
    }
}

fn target_bed_can_stream(plan: &harmonize::HarmonizationPlan) -> bool {
    plan.matches
        .windows(2)
        .all(|pair| pair[0].target_idx < pair[1].target_idx)
}

fn parse_args() -> Result<Option<Config>, String> {
    let mut args = env::args().skip(1).peekable();
    if args.peek().is_none() {
        return Ok(None);
    }

    let mut bfile1 = None;
    let mut bfile2 = None;
    let mut out1 = None;
    let mut out2 = None;
    let mut out_prefix = None;
    let mut report = None;
    let mut bed_threads = 1usize;
    let mut strict_duplicates = false;
    let mut force = false;
    let mut positional = Vec::new();

    while let Some(arg) = args.next() {
        match arg.as_str() {
            "-h" | "--help" => return Ok(None),
            "--bfile1" | "--source" => bfile1 = Some(next_value(&mut args, &arg)?),
            "--bfile2" | "--target" => bfile2 = Some(next_value(&mut args, &arg)?),
            "--out1" => out1 = Some(next_value(&mut args, &arg)?),
            "--out2" => out2 = Some(next_value(&mut args, &arg)?),
            "--out-prefix" | "--out" => out_prefix = Some(next_value(&mut args, &arg)?),
            "--report" => report = Some(PathBuf::from(next_value(&mut args, &arg)?)),
            "--bed-threads" => {
                bed_threads = next_value(&mut args, &arg)?
                    .parse::<usize>()
                    .map_err(|_| "--bed-threads expects 1 or 2".to_string())?;
                if !(1..=2).contains(&bed_threads) {
                    return Err("--bed-threads expects 1 or 2".to_string());
                }
            }
            "--strict-duplicates" => strict_duplicates = true,
            "--force" => force = true,
            _ if arg.starts_with('-') => return Err(format!("unknown option: {arg}\n\n{HELP}")),
            _ => positional.push(arg),
        }
    }

    if !positional.is_empty() {
        match positional.len() {
            3 => {
                bfile1.get_or_insert_with(|| positional[0].clone());
                bfile2.get_or_insert_with(|| positional[1].clone());
                out_prefix.get_or_insert_with(|| positional[2].clone());
            }
            4 => {
                bfile1.get_or_insert_with(|| positional[0].clone());
                bfile2.get_or_insert_with(|| positional[1].clone());
                out1.get_or_insert_with(|| positional[2].clone());
                out2.get_or_insert_with(|| positional[3].clone());
            }
            _ => {
                return Err(format!(
                    "expected 3 or 4 positional arguments, got {}\n\n{HELP}",
                    positional.len()
                ));
            }
        }
    }

    let bfile1 = bfile1.ok_or_else(|| format!("missing --bfile1\n\n{HELP}"))?;
    let bfile2 = bfile2.ok_or_else(|| format!("missing --bfile2\n\n{HELP}"))?;

    let (out1, out2) = match (out1, out2, out_prefix) {
        (Some(a), Some(b), None) => (a, b),
        (None, None, Some(prefix)) => (format!("{prefix}.bfile1"), format!("{prefix}.bfile2")),
        (Some(_), Some(_), Some(_)) => {
            return Err("use either --out1/--out2 or --out-prefix, not both".to_string());
        }
        _ => {
            return Err(format!(
                "missing output prefix: provide --out-prefix or both --out1 and --out2\n\n{HELP}"
            ));
        }
    };

    Ok(Some(Config {
        bfile1,
        bfile2,
        out1,
        out2,
        report,
        bed_threads,
        strict_duplicates,
        force,
    }))
}

fn next_value<I>(args: &mut std::iter::Peekable<I>, flag: &str) -> Result<String, String>
where
    I: Iterator<Item = String>,
{
    args.next()
        .filter(|v| !v.starts_with('-'))
        .ok_or_else(|| format!("{flag} requires a value"))
}

#[cfg(test)]
mod tests {
    use super::*;
    use plink::BED_MAGIC;
    use std::fs;

    #[test]
    fn end_to_end_harmonizes_two_small_bfiles() {
        let unique = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let dir = env::temp_dir().join(format!("bfile_harm_e2e_{unique}"));
        fs::create_dir_all(&dir).unwrap();

        let b1 = dir.join("b1");
        let b2 = dir.join("b2");
        let out1 = dir.join("out1");
        let out2 = dir.join("out2");

        fs::write(
            path_with_ext(b1.to_str().unwrap(), "fam"),
            "F1 I1 0 0 1 -9\nF2 I2 0 0 2 -9\n",
        )
        .unwrap();
        fs::write(
            path_with_ext(b2.to_str().unwrap(), "fam"),
            "G1 J1 0 0 1 -9\nG2 J2 0 0 2 -9\n",
        )
        .unwrap();
        fs::write(
            path_with_ext(b1.to_str().unwrap(), "bim"),
            "1 rs1 0 101 A C\n1 rs2 0 102 G T\n1 rs3 0 103 C G\n",
        )
        .unwrap();
        fs::write(
            path_with_ext(b2.to_str().unwrap(), "bim"),
            "1 rs2 0 202 T G\n1 rs1 0 201 A C\n1 rs3 0 203 A G\n",
        )
        .unwrap();
        fs::write(
            path_with_ext(b1.to_str().unwrap(), "bed"),
            [BED_MAGIC[0], BED_MAGIC[1], BED_MAGIC[2], 0x0c, 0x03, 0x00],
        )
        .unwrap();
        fs::write(
            path_with_ext(b2.to_str().unwrap(), "bed"),
            [BED_MAGIC[0], BED_MAGIC[1], BED_MAGIC[2], 0x0c, 0x08, 0x00],
        )
        .unwrap();

        let config = Config {
            bfile1: b1.to_string_lossy().into_owned(),
            bfile2: b2.to_string_lossy().into_owned(),
            out1: out1.to_string_lossy().into_owned(),
            out2: out2.to_string_lossy().into_owned(),
            report: None,
            bed_threads: 1,
            strict_duplicates: false,
            force: false,
        };
        run_with_config(&config).unwrap();

        assert_eq!(
            fs::read_to_string(path_with_ext(out1.to_str().unwrap(), "bim")).unwrap(),
            "1\trs1\t0\t101\tA\tC\n1\trs2\t0\t102\tG\tT\n"
        );
        assert_eq!(
            fs::read_to_string(path_with_ext(out2.to_str().unwrap(), "bim")).unwrap(),
            "1\trs1\t0\t201\tA\tC\n1\trs2\t0\t202\tG\tT\n"
        );
        assert_eq!(
            fs::read(path_with_ext(out1.to_str().unwrap(), "bed")).unwrap(),
            vec![BED_MAGIC[0], BED_MAGIC[1], BED_MAGIC[2], 0x0c, 0x03]
        );
        assert_eq!(
            fs::read(path_with_ext(out2.to_str().unwrap(), "bed")).unwrap(),
            vec![BED_MAGIC[0], BED_MAGIC[1], BED_MAGIC[2], 0x08, 0x03]
        );

        fs::remove_dir_all(&dir).unwrap();
    }
}
