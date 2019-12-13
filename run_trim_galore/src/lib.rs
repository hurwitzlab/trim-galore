extern crate clap;
extern crate regex;

use clap::{App, Arg};
use regex::Regex;
use std::collections::HashMap;
use std::error::Error;
use std::process::{Command, Stdio};
use std::{
    env,
    fs::{self, DirBuilder},
    io::Write,
    path::{Path, PathBuf},
};

#[derive(Debug)]
struct SplitPath {
    stem: String,
    ext: Option<String>,
}

#[derive(Debug)]
pub struct Config {
    query: Vec<String>,
    out_dir: PathBuf,
    num_concurrent_jobs: Option<u32>,
    num_halt: Option<u32>,
    phred_64: Option<bool>,
    quality: Option<u32>,
    adapter: Option<String>,
    adapter2: Option<String>,
    illumina: Option<bool>,
    nextera: Option<bool>,
    small_rna: Option<bool>,
    consider_already_trimmed: Option<bool>,
    max_length: Option<u32>,
    stringency: Option<u32>,
    error_rate: Option<f32>,
    gzip: Option<bool>,
    dont_gzip: Option<bool>,
    length: Option<u32>,
    max_n: Option<u32>,
    trim_n: Option<bool>,
    clip_r1: Option<u32>,
    clip_r2: Option<u32>,
    three_prime_clip_r1: Option<u32>,
    three_prime_clip_r2: Option<u32>,
    nextseq: Option<u32>,
    basename: Option<String>,
    hardtrim5: Option<u32>,
    hardtrim3: Option<u32>,
    rbbs: Option<bool>,
    non_directional: Option<bool>,
    keep: Option<bool>,
    trim: Option<bool>,
    retain_unpaired: Option<bool>,
    length_1: Option<u32>,
    length_2: Option<u32>,
}

#[derive(Debug, PartialEq, Eq, Hash)]
enum ReadDirection {
    Forward,
    Reverse,
}

type MyResult<T> = Result<T, Box<dyn Error>>;
type ReadPair = HashMap<ReadDirection, String>;
type ReadPairLookup = HashMap<String, ReadPair>;
type SingleReads = Vec<String>;

// --------------------------------------------------
pub fn get_args() -> MyResult<Config> {
    let matches = App::new("run_trim_galore")
        .version("0.1.0")
        .author("Ken Youens-Clark <kyclark@email.arizona.edu>")
        .about("Runs TrimGalore")
        .arg(
            Arg::with_name("query")
                .short("Q")
                .long("query")
                .value_name("FILE_OR_DIR")
                .help("File input or directory")
                .required(true)
                .min_values(1),
        )
        .arg(
            Arg::with_name("out_dir")
                .short("o")
                .long("out_dir")
                .value_name("DIR")
                .help("Output directory"),
        )
        .arg(
            Arg::with_name("num_concurrent_jobs")
                .short("J")
                .long("num_concurrent_jobs")
                .value_name("INT")
                .default_value("8")
                .help("Number of concurrent jobs for parallel"),
        )
        .arg(
            Arg::with_name("num_halt")
                .short("H")
                .long("num_halt")
                .value_name("INT")
                .default_value("1")
                .help("Halt after this many failing jobs"),
        )
        .arg(
            Arg::with_name("phred_64")
                .long("--phred64")
                .help("ASCII+64 quality scores as Phred scores"),
        )
        .arg(
            Arg::with_name("quality")
                .short("q")
                .long("quality")
                .value_name("INT")
                .help("Minimum Phred score"),
        )
        .arg(
            Arg::with_name("adapter")
                .short("a")
                .long("adapter")
                .value_name("STR")
                .help("Adapter sequence to be trimmed"),
        )
        .arg(
            Arg::with_name("adapter2")
                .long("adapter2")
                .value_name("STR")
                .help("Adapter sequence to be trimmed"),
        )
        .arg(
            Arg::with_name("illumina")
                .long("illumina")
                .help("Illumina universal adapter"),
        )
        .arg(
            Arg::with_name("nextera")
                .long("nextera")
                .help("Nextera adapter"),
        )
        .arg(
            Arg::with_name("small_rna")
                .long("small_rna")
                .help("Illumina Small RNA 3' Adapter"),
        )
        .arg(
            Arg::with_name("consider_already_trimmed")
                .long("--consider_already_trimmed")
                .help("threshold to consider already adapter-trimmed"),
        )
        .arg(
            Arg::with_name("max_length")
                .long("max_length")
                .value_name("INT")
                .help("Maximum length"),
        )
        .arg(
            Arg::with_name("stringency")
                .long("stringency")
                .value_name("INT")
                .help(
                    "Overlap with adapter sequence required to trim a sequence",
                ),
        )
        .arg(
            Arg::with_name("error_rate")
                .short("e")
                .long("error_rate")
                .value_name("FLOAT")
                .help("Maximum allowed error rate"),
        )
        .arg(
            Arg::with_name("gzip")
                .long("gzip")
                .help("Compress the output file with GZIP"),
        )
        .arg(
            Arg::with_name("dont_gzip")
                .long("dont_gzip")
                .conflicts_with("gzip")
                .help("Do not compress the output file with GZIP"),
        )
        .arg(
            Arg::with_name("length")
                .long("length")
                .value_name("INT")
                .help("Minimum length"),
        )
        .arg(
            Arg::with_name("max_n")
                .long("max_n")
                .value_name("INT")
                .help("Maximum allowed Ns before discarding"),
        )
        .arg(
            Arg::with_name("trim_n")
                .long("trim_n")
                .help("Remove Ns from either side of the read"),
        )
        .arg(
            Arg::with_name("clip_r1")
                .long("clip_r1")
                .value_name("INT")
                .help("Remove N bp from the 5' end of read 1"),
        )
        .arg(
            Arg::with_name("clip_r2")
                .long("clip_r2")
                .value_name("INT")
                .help("Remove N bp from the 5' end of read 2"),
        )
        .arg(
            Arg::with_name("three_prime_clip_r1")
                .long("three_prime_clip_r1")
                .value_name("INT")
                .help("Remove N bp from the 3' end of read 1"),
        )
        .arg(
            Arg::with_name("three_prime_clip_r2")
                .long("three_prime_clip_r2")
                .value_name("INT")
                .help("Remove N bp from the 3' end of read 2"),
        )
        .arg(
            Arg::with_name("nextseq")
                .long("nextseq")
                .value_name("INT")
                .help("Quality cutoff"),
        )
        .arg(
            Arg::with_name("basename")
                .long("basename")
                .value_name("STR")
                .help("Preferred basename"),
        )
        .arg(
            Arg::with_name("hardtrim5")
                .long("hardtrim5")
                .value_name("INT")
                .help("Hard-trim sequences to N bp at the 5'-end"),
        )
        .arg(
            Arg::with_name("hardtrim3")
                .long("hardtrim3")
                .value_name("INT")
                .help("Hard-trim sequences to N bp at the 3'-end"),
        )
        .arg(
            Arg::with_name("rbbs")
                .long("rbbs")
                .help("Input file was an MspI digested RRBS sample"),
        )
        .arg(
            Arg::with_name("non_directional")
                .long("non_directional")
                .help("Non-directional RRBS libraries"),
        )
        .arg(
            Arg::with_name("keep")
                .long("keep")
                .help("Keep the quality trimmed intermediate file"),
        )
        .arg(
            Arg::with_name("trim")
                .long("trim")
                .help("Trims 1 bp off every read from its 3' end"),
        )
        .arg(
            Arg::with_name("retain_unpaired")
                .long("retain_unpaired")
                .help("Retain unpaired reads"),
        )
        .arg(
            Arg::with_name("length_1")
                .long("length_1")
                .value_name("INT")
                .help(
                    "Unpaired single-end read length cutoff needed for read 1",
                ),
        )
        .arg(
            Arg::with_name("length_2")
                .long("length_2")
                .value_name("INT")
                .help(
                    "Unpaired single-end read length cutoff needed for read 2",
                ),
        )
        .get_matches();

    let out_dir = match matches.value_of("out_dir") {
        Some(x) => PathBuf::from(x),
        _ => {
            let cwd = env::current_dir()?;
            cwd.join(PathBuf::from("trim-galore-out"))
        }
    };

    let num_concurrent_jobs = matches
        .value_of("num_concurrent_jobs")
        .and_then(|x| x.trim().parse::<u32>().ok());

    let num_halt = matches
        .value_of("num_halt")
        .and_then(|x| x.trim().parse::<u32>().ok());

    let phred_64 = Some(matches.is_present("phred_64"));

    let quality = matches
        .value_of("quality")
        .and_then(|x| x.trim().parse::<u32>().ok());

    let adapter = matches.value_of("adapter").map(|a| a.to_string());

    let adapter2 = matches.value_of("adapter2").map(|a| a.to_string());

    let max_length = matches
        .value_of("max_length")
        .and_then(|x| x.trim().parse::<u32>().ok());

    let stringency = matches
        .value_of("stringency")
        .and_then(|x| x.trim().parse::<u32>().ok());

    let error_rate = matches
        .value_of("error_rate")
        .and_then(|x| x.trim().parse::<f32>().ok());

    let length = matches
        .value_of("length")
        .and_then(|x| x.trim().parse::<u32>().ok());

    let max_n = matches
        .value_of("max_n")
        .and_then(|x| x.trim().parse::<u32>().ok());

    let clip_r1 = matches
        .value_of("clip_r1")
        .and_then(|x| x.trim().parse::<u32>().ok());

    let clip_r2 = matches
        .value_of("clip_r1")
        .and_then(|x| x.trim().parse::<u32>().ok());

    let three_prime_clip_r1 = matches
        .value_of("three_prime_clip_r1")
        .and_then(|x| x.trim().parse::<u32>().ok());

    let three_prime_clip_r2 = matches
        .value_of("three_prime_clip_r2")
        .and_then(|x| x.trim().parse::<u32>().ok());

    let nextseq = matches
        .value_of("nextseq")
        .and_then(|x| x.trim().parse::<u32>().ok());

    let basename = matches.value_of("adapter").map(|a| a.to_string());

    let hardtrim5 = matches
        .value_of("hardtrim5")
        .and_then(|x| x.trim().parse::<u32>().ok());

    let hardtrim3 = matches
        .value_of("hardtrim3")
        .and_then(|x| x.trim().parse::<u32>().ok());

    let length_1 = matches
        .value_of("length_1")
        .and_then(|x| x.trim().parse::<u32>().ok());

    let length_2 = matches
        .value_of("length_2")
        .and_then(|x| x.trim().parse::<u32>().ok());

    Ok(Config {
        query: matches.values_of_lossy("query").unwrap(),
        out_dir,
        num_concurrent_jobs,
        num_halt,
        phred_64,
        quality,
        adapter,
        adapter2,
        illumina: Some(matches.is_present("illumina")),
        nextera: Some(matches.is_present("nextera")),
        small_rna: Some(matches.is_present("small_rna")),
        consider_already_trimmed: Some(
            matches.is_present("consider_already_trimmed"),
        ),
        max_length,
        stringency,
        error_rate,
        gzip: Some(matches.is_present("gzip")),
        dont_gzip: Some(matches.is_present("dont_gzip")),
        length,
        max_n,
        trim_n: Some(matches.is_present("trim_n")),
        clip_r1,
        clip_r2,
        three_prime_clip_r1,
        three_prime_clip_r2,
        nextseq,
        basename,
        hardtrim5,
        hardtrim3,
        rbbs: Some(matches.is_present("rbbs")),
        non_directional: Some(matches.is_present("non_directional")),
        keep: Some(matches.is_present("keep")),
        trim: Some(matches.is_present("trim")),
        retain_unpaired: Some(matches.is_present("retain_unpaired")),
        length_1,
        length_2,
    })
}

// --------------------------------------------------
pub fn run(config: Config) -> MyResult<()> {
    let files = find_files(&config.query)?;

    if files.is_empty() {
        let msg = format!("No input files from query \"{:?}\"", &config.query);
        return Err(From::from(msg));
    }

    let (pairs, singles) = classify(&files)?;

    println!(
        "Processing {} pair, {} single.",
        pairs.keys().len(),
        singles.len()
    );

    let jobs = make_jobs(&config, pairs, singles)?;

    run_jobs(
        &jobs,
        "Running trim_galore",
        config.num_concurrent_jobs.unwrap_or(8),
        config.num_halt.unwrap_or(1),
    )?;

    println!("Done, see output in \"{}\"", &config.out_dir.display());

    Ok(())
}

// --------------------------------------------------
fn make_jobs(
    config: &Config,
    pairs: ReadPairLookup,
    singles: SingleReads,
) -> Result<Vec<String>, Box<dyn Error>> {
    let mut args: Vec<String> = vec![];

    if config.phred_64.unwrap_or(false) {
        args.push("--phred64".to_string());
    }

    if let Some(quality) = config.quality {
        args.push(format!("--quality {}", quality));
    }

    if let Some(adapter) = &config.adapter {
        args.push(format!("--adapter {}", adapter));
    }

    if let Some(adapter2) = &config.adapter2 {
        args.push(format!("--adapter2 {}", adapter2));
    }

    if config.illumina.unwrap_or(false) {
        args.push("--illumina".to_string());
    }

    if config.nextera.unwrap_or(false) {
        args.push("--nextera".to_string());
    }

    if config.small_rna.unwrap_or(false) {
        args.push("--small_rna".to_string());
    }

    if config.consider_already_trimmed.unwrap_or(false) {
        args.push("--consider_already_trimmed".to_string());
    }

    if let Some(max_length) = config.max_length {
        args.push(format!("--max_length {}", max_length));
    }

    if let Some(stringency) = config.stringency {
        args.push(format!("--stringency {}", stringency));
    }

    if let Some(error_rate) = config.error_rate {
        args.push(format!("-e {}", error_rate));
    }

    if config.gzip.unwrap_or(false) {
        args.push("--gzip".to_string());
    }

    if config.dont_gzip.unwrap_or(false) {
        args.push("--dont_gzip".to_string());
    }

    if let Some(length) = config.length {
        args.push(format!("--length {}", length));
    }

    if let Some(max_n) = config.max_n {
        args.push(format!("--max_n {}", max_n));
    }

    if config.trim_n.unwrap_or(false) {
        args.push("--trim_n".to_string());
    }

    if let Some(clip_r1) = config.clip_r1 {
        args.push(format!("--clip_r1 {}", clip_r1));
    }

    if let Some(clip_r2) = config.clip_r2 {
        args.push(format!("--clip_r2 {}", clip_r2));
    }

    if let Some(three_prime_clip_r1) = config.three_prime_clip_r1 {
        args.push(format!("--three_prime_clip_r1 {}", three_prime_clip_r1));
    }

    if let Some(three_prime_clip_r2) = config.three_prime_clip_r2 {
        args.push(format!("--three_prime_clip_r2 {}", three_prime_clip_r2));
    }

    if let Some(nextseq) = config.nextseq {
        args.push(format!("--nextseq {}", nextseq));
    }

    if let Some(basename) = &config.basename {
        args.push(format!("--basename {}", basename));
    }

    if let Some(hardtrim5) = config.hardtrim5 {
        args.push(format!("--hardtrim5 {}", hardtrim5));
    }

    if let Some(hardtrim3) = config.hardtrim3 {
        args.push(format!("--hardtrim3 {}", hardtrim3));
    }

    if config.rbbs.unwrap_or(false) {
        args.push("--rbbs".to_string());
    }

    if config.non_directional.unwrap_or(false) {
        args.push("--non_directional".to_string());
    }

    if config.keep.unwrap_or(false) {
        args.push("--keep".to_string());
    }

    if config.trim.unwrap_or(false) {
        args.push("--trim".to_string());
    }

    if config.retain_unpaired.unwrap_or(false) {
        args.push("--retain_unpaired".to_string());
    }

    if let Some(length_1) = config.length_1 {
        args.push(format!("--length_1 {}", length_1));
    }

    if let Some(length_2) = config.length_2 {
        args.push(format!("--length_2 {}", length_2));
    }

    let mut jobs: Vec<String> = vec![];
    for (i, (sample, val)) in pairs.iter().enumerate() {
        println!("{:3}: {}", i + 1, sample);

        if let (Some(fwd), Some(rev)) = (
            val.get(&ReadDirection::Forward),
            val.get(&ReadDirection::Reverse),
        ) {
            let out_dir = &config.out_dir.join(sample);
            if !out_dir.is_dir() {
                DirBuilder::new().recursive(true).create(&out_dir)?;
            }

            jobs.push(format!(
                "trim_galore -o {} {} --paired --fastqc {} {}",
                out_dir.display(),
                args.join(" "),
                fwd,
                rev,
            ));
        }
    }

    for (i, file) in singles.iter().enumerate() {
        let path = Path::new(file);
        let basename = path.file_name().expect("basename");
        let SplitPath { stem, .. } =
            split_filename(basename.to_string_lossy().to_string());
        let out_dir = &config.out_dir.join(&stem);
        if !out_dir.is_dir() {
            DirBuilder::new().recursive(true).create(&out_dir)?;
        }

        println!("{:3}: Single {}", i + 1, basename);

        jobs.push(format!(
            "trim_galore --fastqc -o {} {} {}",
            out_dir.display(),
            args.join(" "),
            file,
        ));
    }

    Ok(jobs)
}

// --------------------------------------------------
fn find_files(paths: &[String]) -> Result<Vec<String>, Box<dyn Error>> {
    let mut files = vec![];
    for path in paths {
        let meta = fs::metadata(path)?;
        if meta.is_file() {
            files.push(path.to_owned());
        } else {
            for entry in fs::read_dir(path)? {
                let entry = entry?;
                let meta = entry.metadata()?;
                if meta.is_file() {
                    files.push(entry.path().display().to_string());
                }
            }
        };
    }

    if files.is_empty() {
        return Err(From::from("No input files"));
    }

    Ok(files)
}

// --------------------------------------------------
fn classify(
    paths: &[String],
) -> Result<(ReadPairLookup, SingleReads), Box<dyn Error>> {
    let paths = paths.iter().map(Path::new);
    let mut exts: Vec<String> =
        paths.clone().map(get_extension).filter_map(|x| x).collect();
    exts.dedup();

    let dots = Regex::new(r"\.").unwrap();
    let exts: Vec<String> = exts
        .into_iter()
        .map(|x| dots.replace(&x, r"\.").to_string())
        .collect();

    let pattern = format!(r"(.+)[_-][Rr]?([12])?\.(?:{})$", exts.join("|"));
    let re = Regex::new(&pattern).unwrap();
    let mut pairs: ReadPairLookup = HashMap::new();
    let mut singles: Vec<String> = vec![];

    for path in paths.map(Path::new) {
        let path_str = path.to_str().expect("Convert path");

        if let Some(file_name) = path.file_name() {
            let basename = file_name.to_string_lossy();
            if let Some(cap) = re.captures(&basename) {
                let sample_name = &cap[1];
                let direction = if &cap[2] == "1" {
                    ReadDirection::Forward
                } else {
                    ReadDirection::Reverse
                };

                if !pairs.contains_key(sample_name) {
                    let mut pair: ReadPair = HashMap::new();
                    pair.insert(direction, path_str.to_string());
                    pairs.insert(sample_name.to_string(), pair);
                } else if let Some(pair) = pairs.get_mut(sample_name) {
                    pair.insert(direction, path_str.to_string());
                }
            } else {
                singles.push(path_str.to_string());
            }
        }
    }

    let bad: Vec<String> = pairs
        .iter()
        .filter_map(|(k, v)| {
            if !v.contains_key(&ReadDirection::Forward)
                || !v.contains_key(&ReadDirection::Reverse)
            {
                Some(k.to_string())
            } else {
                None
            }
        })
        .collect();

    // Push unpaired samples to the singles
    for key in bad {
        pairs.remove(&key);
        singles.push(key);
    }

    Ok((pairs, singles))
}

// --------------------------------------------------
/// Returns the extension plus optional ".gz"
fn get_extension(path: &Path) -> Option<String> {
    let re = Regex::new(r"\.([^.]+(?:\.gz)?)$").unwrap();
    if let Some(basename) = path.file_name() {
        let basename = basename.to_string_lossy();
        if let Some(cap) = re.captures(&basename) {
            return Some(cap[1].to_string());
        }
    }
    None
}

// --------------------------------------------------
fn run_jobs(
    jobs: &[String],
    msg: &str,
    num_concurrent_jobs: u32,
    num_halt: u32,
) -> MyResult<()> {
    let num_jobs = jobs.len();

    if num_jobs > 0 {
        println!(
            "{} (# {} job{} @ {})",
            msg,
            num_jobs,
            if num_jobs == 1 { "" } else { "s" },
            num_concurrent_jobs,
        );

        let mut args: Vec<String> =
            vec!["-j".to_string(), num_concurrent_jobs.to_string()];

        if num_halt > 0 {
            args.push("--halt".to_string());
            args.push(format!("soon,fail={}", num_halt.to_string()));
        }

        let mut process = Command::new("parallel")
            .args(args)
            .stdin(Stdio::piped())
            .stdout(Stdio::null())
            .spawn()?;

        {
            let stdin = process.stdin.as_mut().expect("Failed to open stdin");
            stdin
                .write_all(jobs.join("\n").as_bytes())
                .expect("Failed to write to stdin");
        }

        let result = process.wait()?;
        if !result.success() {
            return Err(From::from("Failed to run jobs in parallel"));
        }
    }

    Ok(())
}

// --------------------------------------------------
fn split_filename(filename: String) -> SplitPath {
    let re = Regex::new(r"([^/]+?)\.([^.]+(?:\.gz)?)$").unwrap();
    if let Some(cap) = re.captures(&filename) {
        SplitPath {
            stem: cap[1].to_string(),
            ext: Some(cap[2].to_string()),
        }
    } else {
        SplitPath {
            stem: filename,
            ext: None,
        }
    }
}

// --------------------------------------------------
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_get_extension() {
        assert_eq!(
            get_extension(Path::new("foo.fna")),
            Some("fna".to_string())
        );

        assert_eq!(
            get_extension(Path::new("foo.fasta.gz")),
            Some("fasta.gz".to_string())
        );

        assert_eq!(
            get_extension(Path::new("foo.fa.gz")),
            Some("fa.gz".to_string())
        );

        assert_eq!(
            get_extension(Path::new("foo.fasta")),
            Some("fasta".to_string())
        );

        assert_eq!(get_extension(Path::new("foo.fq")), Some("fq".to_string()));

        assert_eq!(get_extension(Path::new("foo")), None);
    }

    #[test]
    fn test_classify() {
        let res = classify(&["ERR1711926.fastq.gz".to_string()]);
        assert!(res.is_ok());

        if let Ok((pairs, singles)) = res {
            assert_eq!(pairs.len(), 0);
            assert_eq!(singles.len(), 1);
        }

        let res = classify(&[
            "/foo/bar/ERR1711926_1.fastq.gz".to_string(),
            "/foo/bar/ERR1711926_2.fastq.gz".to_string(),
            "/foo/bar/ERR1711927-R1.fastq.gz".to_string(),
            "/foo/bar/ERR1711927_R2.fastq.gz".to_string(),
            "/foo/bar/ERR1711928.fastq.gz".to_string(),
            "/foo/bar/ERR1711929_1.fastq.gz".to_string(),
        ]);
        assert!(res.is_ok());

        if let Ok((pairs, singles)) = res {
            assert_eq!(pairs.len(), 2);
            assert_eq!(singles.len(), 2);

            assert!(pairs.contains_key("ERR1711926"));
            assert!(pairs.contains_key("ERR1711927"));

            //assert!(!singles.contains_key("ERR1711928"));
            //assert!(!singles.contains_key("ERR1711929"));

            if let Some(val) = pairs.get("ERR1711926") {
                assert!(val.contains_key(&ReadDirection::Forward));
                assert!(val.contains_key(&ReadDirection::Reverse));

                if let Some(fwd) = val.get(&ReadDirection::Forward) {
                    assert_eq!(fwd, &"/foo/bar/ERR1711926_1.fastq.gz");
                }
                if let Some(rev) = val.get(&ReadDirection::Reverse) {
                    assert_eq!(rev, &"/foo/bar/ERR1711926_2.fastq.gz");
                }
            }

            if let Some(val) = pairs.get("ERR1711927") {
                assert!(val.contains_key(&ReadDirection::Forward));
                assert!(val.contains_key(&ReadDirection::Reverse));

                if let Some(fwd) = val.get(&ReadDirection::Forward) {
                    assert_eq!(fwd, &"/foo/bar/ERR1711927-R1.fastq.gz");
                }
                if let Some(rev) = val.get(&ReadDirection::Reverse) {
                    assert_eq!(rev, &"/foo/bar/ERR1711927_R2.fastq.gz");
                }
            }
        }
    }

    #[test]
    fn test_split_filename() {
        let res = split_filename("foo.fa".to_string());
        assert_eq!(res.stem, "foo".to_string());
        assert_eq!(res.ext, Some("fa".to_string()));

        let res = split_filename("foo.fa.gz".to_string());
        assert_eq!(res.stem, "foo".to_string());
        assert_eq!(res.ext, Some("fa.gz".to_string()));

        let res = split_filename("/foo/bar/baz.fa.gz".to_string());
        assert_eq!(res.stem, "baz".to_string());
        assert_eq!(res.ext, Some("fa.gz".to_string()));
    }
}
