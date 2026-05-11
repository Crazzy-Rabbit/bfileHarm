use std::fs::{self, File};
use std::io::{self, BufRead, BufReader, BufWriter, Read, Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};

pub const BED_MAGIC: [u8; 3] = [0x6c, 0x1b, 0x01];
const BED_IO_BUFFER_BYTES: usize = 64 * 1024 * 1024;
const BED_SKIP: u8 = 0;
const BED_COPY: u8 = 1;
const BED_FLIP: u8 = 2;

#[derive(Clone, Debug)]
pub struct BimRecord {
    pub chr: String,
    pub snp: String,
    pub cm: String,
    pub bp: String,
    pub a1: String,
    pub a2: String,
}

#[derive(Debug)]
pub struct FamFile {
    pub content: String,
    pub n_samples: usize,
}

pub fn path_with_ext(prefix: &str, ext: &str) -> PathBuf {
    PathBuf::from(format!("{prefix}.{ext}"))
}

pub fn read_fam(path: &Path) -> io::Result<FamFile> {
    let content = fs::read_to_string(path)?;
    let mut n_samples = 0;
    for (line_no, line) in content.lines().enumerate() {
        if line.trim().is_empty() {
            continue;
        }
        let n_fields = line.split_whitespace().count();
        if n_fields < 6 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "{}:{} has fewer than 6 FAM columns",
                    path.display(),
                    line_no + 1
                ),
            ));
        }
        n_samples += 1;
    }
    Ok(FamFile { content, n_samples })
}

pub fn scan_bim<F>(path: &Path, mut f: F) -> io::Result<usize>
where
    F: FnMut(usize, BimRecord) -> io::Result<()>,
{
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut idx = 0usize;

    for (line_no, line) in reader.lines().enumerate() {
        let line = line?;
        if line.trim().is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split_whitespace().collect();
        if fields.len() < 6 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "{}:{} has fewer than 6 BIM columns",
                    path.display(),
                    line_no + 1
                ),
            ));
        }
        f(
            idx,
            BimRecord {
                chr: fields[0].to_string(),
                snp: fields[1].to_string(),
                cm: fields[2].to_string(),
                bp: fields[3].to_string(),
                a1: fields[4].to_ascii_uppercase(),
                a2: fields[5].to_ascii_uppercase(),
            },
        )?;
        idx += 1;
    }

    Ok(idx)
}

pub fn ensure_outputs_available(
    out1: &str,
    out2: &str,
    report: Option<&Path>,
    force: bool,
) -> io::Result<()> {
    for prefix in [out1, out2] {
        for ext in ["bed", "bim", "fam"] {
            let path = path_with_ext(prefix, ext);
            if !force && path.exists() {
                return Err(io::Error::new(
                    io::ErrorKind::AlreadyExists,
                    format!(
                        "{} already exists; use --force to overwrite",
                        path.display()
                    ),
                ));
            }
        }
    }
    if let Some(path) = report {
        if !force && path.exists() {
            return Err(io::Error::new(
                io::ErrorKind::AlreadyExists,
                format!(
                    "{} already exists; use --force to overwrite",
                    path.display()
                ),
            ));
        }
    }
    Ok(())
}

pub fn write_fam(path: &Path, content: &str, force: bool) -> io::Result<()> {
    write_text(path, content, force)
}

pub fn write_bed_subset<I>(
    input_bed: &Path,
    output_bed: &Path,
    n_samples: usize,
    input_snp_count: usize,
    selected: I,
    force: bool,
) -> io::Result<()>
where
    I: IntoIterator<Item = (usize, bool)>,
{
    let bytes_per_snp = n_samples.div_ceil(4);
    if bytes_per_snp == 0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "FAM file contains zero samples",
        ));
    }

    let expected_len = 3_u64 + (bytes_per_snp as u64) * (input_snp_count as u64);
    let actual_len = fs::metadata(input_bed)?.len();
    if actual_len < expected_len {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "{} is too short for {} SNPs and {} samples",
                input_bed.display(),
                input_snp_count,
                n_samples
            ),
        ));
    }

    let mut input = File::open(input_bed)?;
    validate_bed_magic(&mut input, input_bed)?;

    let mut output = create_output(output_bed, force)?;
    output.write_all(&BED_MAGIC)?;

    let max_snps_per_read = (BED_IO_BUFFER_BYTES / bytes_per_snp).max(1);
    let mut run_start = 0usize;
    let mut run_flips = Vec::with_capacity(max_snps_per_read.min(1 << 20));
    let mut buffer = Vec::new();

    for (snp_idx, flip) in selected {
        if snp_idx >= input_snp_count {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!(
                    "SNP index {snp_idx} is out of bounds for {}",
                    input_bed.display()
                ),
            ));
        }

        let expected_next = run_start + run_flips.len();
        if run_flips.is_empty() {
            run_start = snp_idx;
        } else if snp_idx != expected_next || run_flips.len() >= max_snps_per_read {
            flush_bed_run(
                &mut input,
                &mut output,
                bytes_per_snp,
                n_samples,
                run_start,
                &run_flips,
                &mut buffer,
            )?;
            run_flips.clear();
            run_start = snp_idx;
        }
        run_flips.push(flip);
    }

    if !run_flips.is_empty() {
        flush_bed_run(
            &mut input,
            &mut output,
            bytes_per_snp,
            n_samples,
            run_start,
            &run_flips,
            &mut buffer,
        )?;
    }
    output.flush()
}

pub fn make_bed_action_map<I>(input_snp_count: usize, selected: I) -> io::Result<Vec<u8>>
where
    I: IntoIterator<Item = (usize, bool)>,
{
    let mut actions = vec![BED_SKIP; input_snp_count];
    for (snp_idx, flip) in selected {
        if snp_idx >= input_snp_count {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("SNP index {snp_idx} is out of bounds"),
            ));
        }
        if actions[snp_idx] != BED_SKIP {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("SNP index {snp_idx} appears more than once"),
            ));
        }
        actions[snp_idx] = if flip { BED_FLIP } else { BED_COPY };
    }
    Ok(actions)
}

pub fn write_bed_action_map(
    input_bed: &Path,
    output_bed: &Path,
    n_samples: usize,
    actions: &[u8],
    force: bool,
) -> io::Result<()> {
    let bytes_per_snp = n_samples.div_ceil(4);
    if bytes_per_snp == 0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "FAM file contains zero samples",
        ));
    }

    let input_snp_count = actions.len();
    let expected_len = 3_u64 + (bytes_per_snp as u64) * (input_snp_count as u64);
    let actual_len = fs::metadata(input_bed)?.len();
    if actual_len < expected_len {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "{} is too short for {} SNPs and {} samples",
                input_bed.display(),
                input_snp_count,
                n_samples
            ),
        ));
    }

    let mut input = File::open(input_bed)?;
    validate_bed_magic(&mut input, input_bed)?;
    let mut output = create_output(output_bed, force)?;
    output.write_all(&BED_MAGIC)?;

    let max_snps_per_read = (BED_IO_BUFFER_BYTES / bytes_per_snp).max(1);
    let mut buffer = Vec::new();
    let mut raw_idx = 0usize;

    while raw_idx < input_snp_count {
        let run_snps = max_snps_per_read.min(input_snp_count - raw_idx);
        let run_bytes = run_snps * bytes_per_snp;
        buffer.resize(run_bytes, 0);
        input.read_exact(buffer.as_mut_slice())?;

        for (chunk_idx, chunk) in buffer.chunks_exact_mut(bytes_per_snp).enumerate() {
            match actions[raw_idx + chunk_idx] {
                BED_SKIP => {}
                BED_COPY => {
                    clear_bed_padding(chunk, n_samples);
                    output.write_all(chunk)?;
                }
                BED_FLIP => {
                    flip_bed_chunk(chunk, n_samples);
                    output.write_all(chunk)?;
                }
                action => {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        format!("invalid BED action code: {action}"),
                    ));
                }
            }
        }
        raw_idx += run_snps;
    }
    output.flush()
}

fn flush_bed_run(
    input: &mut File,
    output: &mut BufWriter<File>,
    bytes_per_snp: usize,
    n_samples: usize,
    run_start: usize,
    flips: &[bool],
    buffer: &mut Vec<u8>,
) -> io::Result<()> {
    let n_bytes = bytes_per_snp * flips.len();
    buffer.resize(n_bytes, 0);

    let offset = 3_u64 + (run_start as u64) * (bytes_per_snp as u64);
    input.seek(SeekFrom::Start(offset))?;
    input.read_exact(buffer.as_mut_slice())?;

    for (chunk, flip) in buffer.chunks_exact_mut(bytes_per_snp).zip(flips.iter()) {
        if *flip {
            flip_bed_chunk(chunk, n_samples);
        } else {
            clear_bed_padding(chunk, n_samples);
        }
    }

    output.write_all(buffer.as_slice())
}

fn validate_bed_magic(input: &mut File, path: &Path) -> io::Result<()> {
    let mut magic = [0_u8; 3];
    input.read_exact(&mut magic)?;
    if magic != BED_MAGIC {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("{} is not PLINK SNP-major BED format", path.display()),
        ));
    }
    Ok(())
}

fn flip_bed_chunk(chunk: &mut [u8], n_samples: usize) {
    for sample in 0..n_samples {
        let byte_idx = sample / 4;
        let shift = (sample % 4) * 2;
        let code = (chunk[byte_idx] >> shift) & 0b11;
        let flipped = match code {
            0b00 => 0b11,
            0b11 => 0b00,
            0b10 => 0b10,
            0b01 => 0b01,
            _ => unreachable!(),
        };
        chunk[byte_idx] = (chunk[byte_idx] & !(0b11 << shift)) | (flipped << shift);
    }
    clear_bed_padding(chunk, n_samples);
}

fn clear_bed_padding(chunk: &mut [u8], n_samples: usize) {
    for sample in n_samples..(chunk.len() * 4) {
        let byte_idx = sample / 4;
        let shift = (sample % 4) * 2;
        chunk[byte_idx] &= !(0b11 << shift);
    }
}

pub fn write_text(path: &Path, content: &str, force: bool) -> io::Result<()> {
    let mut out = create_output(path, force)?;
    out.write_all(content.as_bytes())?;
    out.flush()
}

pub fn create_output(path: &Path, force: bool) -> io::Result<BufWriter<File>> {
    if !force && path.exists() {
        return Err(io::Error::new(
            io::ErrorKind::AlreadyExists,
            format!(
                "{} already exists; use --force to overwrite",
                path.display()
            ),
        ));
    }
    if let Some(parent) = path.parent() {
        if !parent.as_os_str().is_empty() {
            fs::create_dir_all(parent)?;
        }
    }
    File::create(path).map(|file| BufWriter::with_capacity(BED_IO_BUFFER_BYTES, file))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn flip_chunk_recodes_plink_genotypes() {
        let original = (3_u8 << 2) | (2_u8 << 4) | (1_u8 << 6);
        let expected = 3_u8 | (2_u8 << 4) | (1_u8 << 6);
        let mut chunk = vec![original];
        flip_bed_chunk(&mut chunk, 4);
        assert_eq!(chunk[0], expected);
    }

    #[test]
    fn flip_chunk_clears_padding_bits() {
        let mut chunk = vec![0xff];
        flip_bed_chunk(&mut chunk, 3);
        assert_eq!(chunk[0] >> 6, 0);
    }

    #[test]
    fn write_bed_subset_copies_flips_and_clears_padding() {
        let unique = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let dir = std::env::temp_dir().join(format!("bfile_harm_test_{unique}"));
        fs::create_dir_all(&dir).unwrap();

        let input = dir.join("in.bed");
        let output = dir.join("out.bed");
        fs::write(
            &input,
            [
                BED_MAGIC[0],
                BED_MAGIC[1],
                BED_MAGIC[2],
                0b1111_1100,
                0b0000_0000,
                0b1111_1100,
            ],
        )
        .unwrap();

        write_bed_subset(&input, &output, 2, 3, vec![(0, false), (2, true)], false).unwrap();

        let written = fs::read(&output).unwrap();
        assert_eq!(
            written,
            vec![
                BED_MAGIC[0],
                BED_MAGIC[1],
                BED_MAGIC[2],
                0b0000_1100,
                0b0000_0011
            ]
        );

        fs::remove_dir_all(&dir).unwrap();
    }

    #[test]
    fn action_map_streams_and_flips_bed() {
        let unique = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let dir = std::env::temp_dir().join(format!("bfile_harm_action_test_{unique}"));
        fs::create_dir_all(&dir).unwrap();

        let input = dir.join("in.bed");
        let output = dir.join("out.bed");
        fs::write(
            &input,
            [
                BED_MAGIC[0],
                BED_MAGIC[1],
                BED_MAGIC[2],
                0b1111_1100,
                0b0000_0000,
                0b1111_1100,
            ],
        )
        .unwrap();

        let actions = make_bed_action_map(3, vec![(0, false), (2, true)]).unwrap();
        write_bed_action_map(&input, &output, 2, &actions, false).unwrap();

        let written = fs::read(&output).unwrap();
        assert_eq!(
            written,
            vec![
                BED_MAGIC[0],
                BED_MAGIC[1],
                BED_MAGIC[2],
                0b0000_1100,
                0b0000_0011
            ]
        );

        fs::remove_dir_all(&dir).unwrap();
    }
}
