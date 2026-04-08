mod cli;
mod io;

use crate::cli::Cli;
use crate::io::alignment_start_times_from_path;
use crate::io::Fastx;
use crate::io::TimeExt;
use anyhow::{anyhow, Context, Result};
use clap::Parser;
use env_logger::Builder;
use itertools::Itertools;
use itertools::MinMaxResult::{MinMax, NoElements, OneElement};
use log::info;
use log::LevelFilter;
use ontime::{valid_selection, DurationExt};
use std::io::stdout;
use time::format_description::well_known::Rfc3339;
use time::format_description::FormatItem;
use time::macros::format_description;
use time::{Duration, PrimitiveDateTime};

const TIME_FMT: &[FormatItem<'_>] =
    format_description!("[year]-[month]-[day]T[hour]:[minute]:[second].[subsecond]Z");

#[derive(Debug, PartialEq, Eq, Clone, Copy)]
enum FileFormat {
    Alignment,
    Fastx,
}

enum TimeArg {
    Absent,
    Timestamp(PrimitiveDateTime),
    Relative,
}

fn classify_time_arg(value: &Option<String>) -> TimeArg {
    match value {
        None => TimeArg::Absent,
        Some(s) => PrimitiveDateTime::parse(s, &Rfc3339)
            .map(TimeArg::Timestamp)
            .unwrap_or(TimeArg::Relative),
    }
}

fn parse_relative_arg(value: &Option<String>) -> Result<Option<Duration>> {
    match value {
        None => Ok(None),
        Some(s) => match PrimitiveDateTime::parse(s, &Rfc3339) {
            Ok(_) => Ok(None),
            Err(_) => Ok(Some(Duration::from_str(s)?)),
        },
    }
}

fn main() -> Result<()> {
    let args = Cli::parse();
    // setup logging
    let mut log_builder = Builder::new();
    log_builder
        .filter(None, LevelFilter::Info)
        .format_module_path(false)
        .format_target(false)
        .init();

    let input_format = match &args.input.extension().and_then(|ext| ext.to_str()) {
        Some("sam" | "bam") => FileFormat::Alignment,
        Some("fastq" | "fq" | "fasta" | "fa") => FileFormat::Fastx,
        Some("gz") => {
            let p = args.input.with_extension("");
            let ext = p.extension().and_then(|ext| ext.to_str());
            match ext {
                Some("sam" | "bam") => FileFormat::Alignment,
                Some("fastq" | "fq" | "fasta" | "fa") => FileFormat::Fastx,
                _ => return Err(anyhow!("Unrecognized file extension for input file")),
            }
        }
        _ => return Err(anyhow!("Unrecognized file extension for input file")),
    };

    let output_type = match &args.output {
        None => input_format,
        Some(p) => match &p.extension().and_then(|ext| ext.to_str()) {
            Some("sam" | "bam") => FileFormat::Alignment,
            Some("fastq" | "fq" | "fasta" | "fa") => FileFormat::Fastx,
            Some("gz") => {
                let p = p.with_extension("");
                let ext = p.extension().and_then(|ext| ext.to_str());
                match ext {
                    Some("sam" | "bam") => FileFormat::Alignment,
                    Some("fastq" | "fq" | "fasta" | "fa") => FileFormat::Fastx,
                    _ => return Err(anyhow!("Unrecognized file extension for output file")),
                }
            }
            _ => return Err(anyhow!("Unrecognized file extension for output file")),
        },
    };
    if input_format != output_type {
        return Err(anyhow!("Input and output file formats do not match"));
    }

    let input_fastx = Fastx::from_path(&args.input);

    let absolute_earliest = classify_time_arg(&args.earliest);
    let absolute_latest = classify_time_arg(&args.latest);

    let relative_earliest = parse_relative_arg(&args.earliest)?;
    let relative_latest = parse_relative_arg(&args.latest)?;

    if args.assume_start_time_sorted
        && !args.show
        && !matches!(absolute_earliest, TimeArg::Timestamp(_))
        && !matches!(absolute_latest, TimeArg::Timestamp(_))
        && !relative_earliest.map_or(false, |dur| dur.is_negative())
        && !relative_latest.map_or(false, |dur| dur.is_negative())
        && (relative_earliest.is_some() || relative_latest.is_some())
    {
        info!("Using single-pass extraction for sorted relative bounds");

        let nb_reads_to_keep = match output_type {
            FileFormat::Fastx => {
                let mut output_handle = match &args.output {
                    None => match args.output_type {
                        None => Box::new(stdout()),
                        Some(fmt) => {
                            niffler::basic::get_writer(Box::new(stdout()), fmt, args.compress_level)?
                        }
                    },
                    Some(p) => {
                        let out_fastx = Fastx::from_path(p);
                        out_fastx
                            .create(args.compress_level, args.output_type)
                            .context("Failed to create the output file")?
                    }
                };

                let (nb_reads_seen, nb_reads_written) = input_fastx
                    .extract_reads_relative_to_first_into(
                        relative_earliest,
                        relative_latest,
                        &mut output_handle,
                    )
                    .context("Failed to extract start times")?;
                if nb_reads_seen == 0 {
                    return Err(anyhow!("Did not find any start times in the input"));
                }
                nb_reads_written
            }
            FileFormat::Alignment => {
                let mut writer = match &args.output {
                    None => noodles_util::alignment::io::writer::Builder::default()
                        .build_from_writer(Box::new(stdout()))?,
                    Some(p) => {
                        noodles_util::alignment::io::writer::Builder::default().build_from_path(p)?
                    }
                };

                let mut bam_reader = noodles_util::alignment::io::reader::Builder::default()
                    .build_from_path(&args.input)?;
                let (nb_reads_seen, nb_reads_written) = bam_reader
                    .extract_reads_relative_to_first_into(
                        relative_earliest,
                        relative_latest,
                        &mut writer,
                    )
                    .context("Failed to extract start times")?;
                if nb_reads_seen == 0 {
                    return Err(anyhow!("Did not find any start times in the input"));
                }
                nb_reads_written
            }
        };

        info!("Done! Kept {} reads", nb_reads_to_keep);
        return Ok(());
    }

    if !args.show
        && !matches!(absolute_earliest, TimeArg::Relative)
        && !matches!(absolute_latest, TimeArg::Relative)
        && (matches!(absolute_earliest, TimeArg::Timestamp(_))
            || matches!(absolute_latest, TimeArg::Timestamp(_)))
    {
        let earliest = match absolute_earliest {
            TimeArg::Timestamp(ts) => Some(ts),
            TimeArg::Absent => None,
            TimeArg::Relative => unreachable!(),
        };
        let latest = match absolute_latest {
            TimeArg::Timestamp(ts) => Some(ts),
            TimeArg::Absent => None,
            TimeArg::Relative => unreachable!(),
        };

        if let (Some(earliest), Some(latest)) = (earliest, latest) {
            if latest < earliest {
                return Err(anyhow!(
                    "The earliest timestamp is after the latest timestamp"
                ));
            }
        }

        info!("Using single-pass extraction for absolute timestamp bounds");
        info!(
            "Extracting reads with a start time between {:?} and {:?}...",
            earliest, latest
        );

        let nb_reads_to_keep = match output_type {
            FileFormat::Fastx => {
                let mut output_handle = match &args.output {
                    None => match args.output_type {
                        None => Box::new(stdout()),
                        Some(fmt) => {
                            niffler::basic::get_writer(Box::new(stdout()), fmt, args.compress_level)?
                        }
                    },
                    Some(p) => {
                        let out_fastx = Fastx::from_path(p);
                        out_fastx
                            .create(args.compress_level, args.output_type)
                            .context("Failed to create the output file")?
                    }
                };

                let (nb_reads_seen, nb_reads_written) = input_fastx
                    .extract_reads_between_into(earliest.as_ref(), latest.as_ref(), &mut output_handle)
                    .context("Failed to extract start times")?;
                if nb_reads_seen == 0 {
                    return Err(anyhow!("Did not find any start times in the input"));
                }
                nb_reads_written
            }
            FileFormat::Alignment => {
                let mut writer = match &args.output {
                    None => noodles_util::alignment::io::writer::Builder::default()
                        .build_from_writer(Box::new(stdout()))?,
                    Some(p) => {
                        noodles_util::alignment::io::writer::Builder::default().build_from_path(p)?
                    }
                };

                let mut bam_reader = noodles_util::alignment::io::reader::Builder::default()
                    .build_from_path(&args.input)?;
                let (nb_reads_seen, nb_reads_written) = bam_reader
                    .extract_reads_between_into(earliest.as_ref(), latest.as_ref(), &mut writer)
                    .context("Failed to extract start times")?;
                if nb_reads_seen == 0 {
                    return Err(anyhow!("Did not find any start times in the input"));
                }
                nb_reads_written
            }
        };

        info!("Done! Kept {} reads", nb_reads_to_keep);
        return Ok(());
    }

    info!("Extracting read start times...");

    let start_times = match input_format {
        FileFormat::Fastx => input_fastx.start_times(),
        FileFormat::Alignment => alignment_start_times_from_path(&args.input),
    }
    .context("Failed to extract start times")?;

    if start_times.is_empty() {
        return Err(anyhow!("Did not find any start times in the input"));
    }

    info!("Gathered start times for {} reads", start_times.len());

    // safe to unwrap as we know start times is not empty
    let (first_timestamp, last_timestamp) = match start_times.iter().minmax() {
        NoElements => return Err(anyhow!("No start times in input fastq")),
        OneElement(el) => (*el, *el),
        MinMax(x, y) => (*x, *y),
    };

    if args.show {
        println!("Earliest: {}", first_timestamp.format(TIME_FMT)?);
        println!("Latest  : {}", last_timestamp.format(TIME_FMT)?);
        return Ok(());
    }
    info!(
        "First and last timestamps in the input are {} and {}",
        first_timestamp.format(TIME_FMT)?,
        last_timestamp.format(TIME_FMT)?
    );

    let earliest = match args.earliest {
        None => first_timestamp.to_owned(),
        Some(s) => match PrimitiveDateTime::parse(&s, &Rfc3339) {
            Ok(t) => t,
            Err(_) => {
                let duration = Duration::from_str(&s)?;
                if duration.is_negative() {
                    last_timestamp
                        .checked_add(duration)
                        .context("Subtracting --from from the last timestamp caused an overflow")?
                } else {
                    first_timestamp
                        .checked_add(duration)
                        .context("Adding --from to the first timestamp caused an overflow")?
                }
            }
        },
    };

    let latest = match args.latest {
        None => last_timestamp.to_owned(),
        Some(s) => match PrimitiveDateTime::parse(&s, &Rfc3339) {
            Ok(t) => t,
            Err(_) => {
                let duration = Duration::from_str(&s)?;
                if duration.is_negative() {
                    last_timestamp
                        .checked_add(duration)
                        .context("Subtracting --to from the last timestamp caused an overflow")?
                } else {
                    first_timestamp
                        .checked_add(duration)
                        .context("Adding --to to the first timestamp caused an overflow")?
                }
            }
        },
    };

    if latest < earliest {
        return Err(anyhow!(
            "The earliest timestamp is after the latest timestamp"
        ));
    }

    info!(
        "Extracting reads with a start time between {} and {}...",
        earliest, latest
    );
    let reads_to_keep = valid_selection(&start_times, &earliest, &latest);
    let nb_reads_to_keep = reads_to_keep.keep_count();

    match output_type {
        FileFormat::Fastx => {
            let mut output_handle = match &args.output {
                None => match args.output_type {
                    None => Box::new(stdout()),
                    Some(fmt) => {
                        niffler::basic::get_writer(Box::new(stdout()), fmt, args.compress_level)?
                    }
                },
                Some(p) => {
                    let out_fastx = Fastx::from_path(p);
                    out_fastx
                        .create(args.compress_level, args.output_type)
                        .context("Failed to create the output file")?
                }
            };

            input_fastx.extract_reads_in_timeframe_into(
                &reads_to_keep,
                &mut output_handle,
            )?;
        }
        FileFormat::Alignment => {
            let mut writer = match &args.output {
                None => noodles_util::alignment::io::writer::Builder::default()
                    .build_from_writer(Box::new(stdout()))?,
                Some(p) => {
                    noodles_util::alignment::io::writer::Builder::default().build_from_path(p)?
                }
            };

            let mut bam_reader = noodles_util::alignment::io::reader::Builder::default()
                .build_from_path(&args.input)?;
            bam_reader.extract_reads_in_timeframe_into(
                &reads_to_keep,
                &mut writer,
            )?;
        }
    };

    info!("Done! Kept {} reads", nb_reads_to_keep);

    Ok(())
}
