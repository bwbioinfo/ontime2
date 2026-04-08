use crate::cli::CompressionExt;
use anyhow::anyhow;
use needletail::errors::ParseErrorKind::EmptyFile;
use needletail::parse_fastx_file;
use noodles_bam as bam;
use noodles_bgzf as bgzf;
use noodles_sam::alignment::record::data::field::Tag;
use noodles_util::alignment::io::Writer;
use ontime::{parse_rfc3339_bytes, FastxRecordExt, ReadSelection};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::num::NonZeroUsize;
use std::path::{Path, PathBuf};
use std::time::UNIX_EPOCH;
use thiserror::Error;
use time::format_description::FormatItem;
use time::macros::format_description;
use time::Duration;
use time::PrimitiveDateTime;

const SIDECAR_TIME_FMT: &[FormatItem<'_>] =
    format_description!("[year]-[month]-[day]T[hour]:[minute]:[second].[subsecond]Z");
const SIDECAR_MAGIC: &str = "# ontime-sidecar-v1";

/// A `Struct` used for seamlessly dealing with either compressed or uncompressed fasta/fastq files.
#[derive(Debug, PartialEq, Eq)]
pub struct Fastx {
    /// The path for the file.
    path: PathBuf,
}

/// A collection of custom errors relating to the working with files for this package.
#[derive(Error, Debug)]
pub enum IOError {
    /// Indicates that the specified input file could not be opened/read.
    #[error("Read error")]
    ReadError {
        source: needletail::errors::ParseError,
    },

    /// Indicates that a sequence record could not be parsed.
    #[error("Failed to parse record")]
    ParseError {
        source: needletail::errors::ParseError,
    },

    /// Indicates that the specified output file could not be created.
    #[error("Output file could not be created")]
    CreateError { source: std::io::Error },

    /// The fastq record is missing the start time
    #[error("Missing start_time in fastq record start at line {0}")]
    MissingTime(u64),

    /// Indicates and error trying to create the compressor
    #[error(transparent)]
    CompressOutputError(#[from] niffler::Error),

    /// Indicates that some indices we expected to find in the input file weren't found.
    #[error("Some expected indices were not in the input file")]
    IndicesNotFound,

    /// Indicates that writing to the output file failed.
    #[error("Could not write to output file")]
    WriteError { source: anyhow::Error },

    /// Indicates there was an error reading the header of the input file.
    #[error("Could not read the header of the input file")]
    ReadHeaderError { source: anyhow::Error },

    /// Indicates that the alignment file record could not be parsed.
    #[error("Failed to parse alignment record")]
    ParseAlignmentError { source: anyhow::Error },

    /// Indicates an issue reading or writing the timestamp sidecar.
    #[error("Failed to access timestamp sidecar")]
    SidecarError { source: std::io::Error },
}

fn sidecar_path(path: &Path) -> PathBuf {
    let mut sidecar = path.as_os_str().to_os_string();
    sidecar.push(".ontime-index");
    PathBuf::from(sidecar)
}

fn file_fingerprint(path: &Path) -> std::io::Result<(u64, u128)> {
    let metadata = std::fs::metadata(path)?;
    let modified = metadata
        .modified()?
        .duration_since(UNIX_EPOCH)
        .map_err(|err| std::io::Error::new(std::io::ErrorKind::Other, err))?;
    Ok((metadata.len(), modified.as_nanos()))
}

fn read_sidecar(path: &Path) -> Result<Option<Vec<PrimitiveDateTime>>, IOError> {
    let sidecar = sidecar_path(path);
    if !sidecar.exists() {
        return Ok(None);
    }

    let expected = file_fingerprint(path).map_err(|source| IOError::SidecarError { source })?;
    let reader = BufReader::new(
        File::open(sidecar).map_err(|source| IOError::SidecarError { source })?,
    );
    let mut lines = reader.lines();

    let Some(Ok(magic)) = lines.next() else {
        return Ok(None);
    };
    if magic != SIDECAR_MAGIC {
        return Ok(None);
    }

    let Some(Ok(size_line)) = lines.next() else {
        return Ok(None);
    };
    let Some(size) = size_line.strip_prefix("size=") else {
        return Ok(None);
    };
    let Ok(size) = size.parse::<u64>() else {
        return Ok(None);
    };

    let Some(Ok(mtime_line)) = lines.next() else {
        return Ok(None);
    };
    let Some(mtime) = mtime_line.strip_prefix("mtime_nanos=") else {
        return Ok(None);
    };
    let Ok(mtime) = mtime.parse::<u128>() else {
        return Ok(None);
    };

    if (size, mtime) != expected {
        return Ok(None);
    }

    let mut timestamps = Vec::new();
    for line in lines {
        let line = line.map_err(|source| IOError::SidecarError { source })?;
        let Ok(timestamp) = PrimitiveDateTime::parse(&line, SIDECAR_TIME_FMT) else {
            return Ok(None);
        };
        timestamps.push(timestamp);
    }

    Ok(Some(timestamps))
}

fn write_sidecar(path: &Path, timestamps: &[PrimitiveDateTime]) -> Result<(), IOError> {
    let sidecar = sidecar_path(path);
    let (size, mtime_nanos) =
        file_fingerprint(path).map_err(|source| IOError::SidecarError { source })?;
    let file = File::create(sidecar).map_err(|source| IOError::SidecarError { source })?;
    let mut writer = BufWriter::new(file);

    writeln!(writer, "{SIDECAR_MAGIC}").map_err(|source| IOError::SidecarError { source })?;
    writeln!(writer, "size={size}").map_err(|source| IOError::SidecarError { source })?;
    writeln!(writer, "mtime_nanos={mtime_nanos}")
        .map_err(|source| IOError::SidecarError { source })?;
    for timestamp in timestamps {
        let line = timestamp
            .format(SIDECAR_TIME_FMT)
            .map_err(|err| std::io::Error::new(std::io::ErrorKind::Other, err))
            .map_err(|source| IOError::SidecarError { source })?;
        writeln!(writer, "{line}").map_err(|source| IOError::SidecarError { source })?;
    }

    Ok(())
}

pub fn alignment_start_times_from_path(path: &Path) -> Result<Vec<PrimitiveDateTime>, IOError> {
    if let Some(timestamps) = read_sidecar(path)? {
        return Ok(timestamps);
    }

    let timestamps = match path.extension().and_then(|ext| ext.to_str()) {
        Some("bam") => bam_start_times_from_path_parallel(path)?,
        _ => {
            let mut reader = noodles_util::alignment::io::reader::Builder::default()
                .build_from_path(path)
                .map_err(|source| IOError::ReadHeaderError {
                    source: anyhow::Error::from(source),
                })?;
            reader.start_times()?
        }
    };
    let _ = write_sidecar(path, &timestamps);
    Ok(timestamps)
}

fn bam_start_times_from_path_parallel(path: &Path) -> Result<Vec<PrimitiveDateTime>, IOError> {
    let worker_count = std::thread::available_parallelism()
        .ok()
        .unwrap_or_else(|| NonZeroUsize::new(1).unwrap());
    let file = File::open(path).map_err(|source| IOError::CreateError { source })?;
    let bgzf_reader = bgzf::io::MultithreadedReader::with_worker_count(worker_count, file);
    let mut reader = bam::io::Reader::from(bgzf_reader);
    reader
        .read_header()
        .map_err(|source| IOError::ReadHeaderError {
            source: anyhow::Error::from(source),
        })?;

    let tag = Tag::new(b's', b't');
    let mut start_times = Vec::new();

    for (i, record) in reader.records().enumerate() {
        let record = record.map_err(|source| IOError::ParseAlignmentError {
            source: anyhow! { source.to_string() },
        })?;
        let data = record.data();
        let start_time = data
            .get(&tag)
            .ok_or(IOError::MissingTime(i as u64))?
            .map_err(|_| IOError::MissingTime(i as u64))?;
        let start_time = match start_time {
            noodles_sam::alignment::record::data::field::Value::String(s) => s,
            _ => return Err(IOError::MissingTime(i as u64)),
        };
        let start_time =
            parse_rfc3339_bytes(start_time).ok_or(IOError::MissingTime(i as u64))?;
        start_times.push(start_time);
    }

    Ok(start_times)
}

impl Fastx {
    /// Create a `Fastx` object from a `std::path::Path`.
    ///
    /// # Example
    ///
    /// ```rust
    /// let path = std::path::Path::new("input.fa.gz");
    /// let fastx = Fastx::from_path(path);
    /// ```
    pub fn from_path(path: &Path) -> Self {
        Fastx {
            path: path.to_path_buf(),
        }
    }
    /// Create the file associated with this `Fastx` object for writing.
    ///
    /// # Errors
    /// If the file cannot be created then an `Err` containing a variant of [`FastxError`](#fastxerror) is
    /// returned.
    ///
    /// # Example
    ///
    /// ```rust
    /// let path = std::path::Path::new("output.fa");
    /// let fastx = Fastx{ path };
    /// { // this scoping means the file handle is closed afterwards.
    ///     let file_handle = fastx.create(6, None)?;
    ///     write!(file_handle, ">read1\nACGT\n")?
    /// }
    /// ```
    pub fn create(
        &self,
        compression_lvl: niffler::compression::Level,
        compression_fmt: Option<niffler::compression::Format>,
    ) -> Result<Box<dyn Write>, IOError> {
        let file = File::create(&self.path).map_err(|source| IOError::CreateError { source })?;
        let file_handle = Box::new(BufWriter::new(file));
        let fmt = match compression_fmt {
            None => niffler::Format::from_path(&self.path),
            Some(f) => f,
        };
        niffler::get_writer(file_handle, fmt, compression_lvl).map_err(IOError::CompressOutputError)
    }
    /// Returns a vector containing the start time of each read.
    ///
    /// # Errors
    /// If the file cannot be opened or there is an issue parsing any records then an
    /// `Err` containing a variant of [`IOError`](#ioerror) is returned.
    pub fn start_times(&self) -> Result<Vec<PrimitiveDateTime>, IOError> {
        let mut start_times: Vec<PrimitiveDateTime> = vec![];
        let mut reader = match parse_fastx_file(&self.path) {
            Ok(rdr) => rdr,
            Err(e) if e.kind == EmptyFile => return Ok(start_times),
            Err(source) => return Err(IOError::ReadError { source }),
        };

        while let Some(record) = reader.next() {
            match record {
                Ok(rec) => {
                    let start_time = match rec.start_time() {
                        Some(t) => t,
                        None => return Err(IOError::MissingTime(rec.start_line_number())),
                    };
                    start_times.push(start_time)
                }
                Err(err) => return Err(IOError::ParseError { source: err }),
            }
        }
        Ok(start_times)
    }

    pub fn extract_reads_in_timeframe_into<T: Write>(
        &self,
        selection: &ReadSelection,
        write_to: &mut T,
    ) -> Result<(), IOError> {
        let mut reader =
            parse_fastx_file(&self.path).map_err(|source| IOError::ReadError { source })?;
        let mut read_idx: usize = 0;
        let mut nb_reads_written = 0;
        let nb_reads_keep = selection.keep_count();
        let mut next_sparse_idx = 0;

        while let Some(record) = reader.next() {
            match record {
                Err(source) => return Err(IOError::ParseError { source }),
                Ok(rec) if match selection {
                    ReadSelection::Dense(mask) => mask[read_idx],
                    ReadSelection::Sparse(indices) => {
                        if next_sparse_idx < indices.len() && indices[next_sparse_idx] == read_idx {
                            next_sparse_idx += 1;
                            true
                        } else {
                            false
                        }
                    }
                } => {
                    rec.write(write_to, None)
                        .map_err(|err| IOError::WriteError {
                            source: anyhow::Error::from(err),
                        })?;
                    nb_reads_written += 1;
                    if nb_reads_keep == nb_reads_written {
                        break;
                    }
                }
                Ok(_) => (),
            }

            read_idx += 1;
        }

        if nb_reads_written == nb_reads_keep {
            Ok(())
        } else {
            Err(IOError::IndicesNotFound)
        }
    }

    pub fn extract_reads_between_into<T: Write>(
        &self,
        earliest: Option<&PrimitiveDateTime>,
        latest: Option<&PrimitiveDateTime>,
        write_to: &mut T,
    ) -> Result<(usize, usize), IOError> {
        let mut reader =
            parse_fastx_file(&self.path).map_err(|source| IOError::ReadError { source })?;
        let mut nb_reads_seen = 0;
        let mut nb_reads_written = 0;

        while let Some(record) = reader.next() {
            match record {
                Err(source) => return Err(IOError::ParseError { source }),
                Ok(rec) => {
                    let start_time = rec
                        .start_time()
                        .ok_or(IOError::MissingTime(rec.start_line_number()))?;
                    nb_reads_seen += 1;

                    let keep = earliest.map_or(true, |min| &start_time >= min)
                        && latest.map_or(true, |max| &start_time <= max);

                    if keep {
                        rec.write(write_to, None)
                            .map_err(|err| IOError::WriteError {
                                source: anyhow::Error::from(err),
                            })?;
                        nb_reads_written += 1;
                    }
                }
            }
        }

        Ok((nb_reads_seen, nb_reads_written))
    }

    pub fn extract_reads_relative_to_first_into<T: Write>(
        &self,
        from: Option<Duration>,
        to: Option<Duration>,
        write_to: &mut T,
    ) -> Result<(usize, usize), IOError> {
        let mut reader =
            parse_fastx_file(&self.path).map_err(|source| IOError::ReadError { source })?;
        let mut nb_reads_seen = 0;
        let mut nb_reads_written = 0;
        let mut anchor: Option<PrimitiveDateTime> = None;
        let mut earliest_bound = None;
        let mut latest_bound = None;

        while let Some(record) = reader.next() {
            match record {
                Err(source) => return Err(IOError::ParseError { source }),
                Ok(rec) => {
                    let start_time = rec
                        .start_time()
                        .ok_or(IOError::MissingTime(rec.start_line_number()))?;
                    nb_reads_seen += 1;

                    if anchor.is_none() {
                        anchor = Some(start_time);
                        earliest_bound = from.and_then(|dur| start_time.checked_add(dur));
                        latest_bound = to.and_then(|dur| start_time.checked_add(dur));
                    }

                    if latest_bound.map_or(false, |max| start_time > max) {
                        break;
                    }

                    let keep = earliest_bound.map_or(true, |min| start_time >= min)
                        && latest_bound.map_or(true, |max| start_time <= max);

                    if keep {
                        rec.write(write_to, None)
                            .map_err(|err| IOError::WriteError {
                                source: anyhow::Error::from(err),
                            })?;
                        nb_reads_written += 1;
                    }
                }
            }
        }

        Ok((nb_reads_seen, nb_reads_written))
    }
}

pub trait TimeExt {
    fn start_times(&mut self) -> Result<Vec<PrimitiveDateTime>, IOError>;
    fn extract_reads_relative_to_first_into(
        &mut self,
        from: Option<Duration>,
        to: Option<Duration>,
        writer: &mut Writer,
    ) -> Result<(usize, usize), IOError>;
    fn extract_reads_between_into(
        &mut self,
        earliest: Option<&PrimitiveDateTime>,
        latest: Option<&PrimitiveDateTime>,
        writer: &mut Writer,
    ) -> Result<(usize, usize), IOError>;
    fn extract_reads_in_timeframe_into(
        &mut self,
        selection: &ReadSelection,
        writer: &mut Writer,
    ) -> Result<(), IOError>;
}

impl TimeExt for noodles_util::alignment::io::reader::Reader<Box<dyn BufRead>> {
    fn start_times(&mut self) -> Result<Vec<PrimitiveDateTime>, IOError> {
        let mut start_times: Vec<PrimitiveDateTime> = vec![];
        let header = self
            .read_header()
            .map_err(|source| IOError::ReadHeaderError {
                source: anyhow::Error::from(source),
            })?;
        let records = self.records(&header);
        let tag = Tag::new(b's', b't');

        for (i, record) in records.enumerate() {
            let record = record.map_err(|source| IOError::ParseAlignmentError {
                source: anyhow! { source.to_string() },
            })?;
            let data = record.data();
            let start_time = data
                .get(&tag)
                .ok_or(IOError::MissingTime(i as u64))?
                .map_err(|_| IOError::MissingTime(i as u64))?;
            let start_time = match start_time {
                noodles_sam::alignment::record::data::field::Value::String(s) => s,
                _ => return Err(IOError::MissingTime(i as u64)),
            };
            let start_time = parse_rfc3339_bytes(start_time).ok_or(IOError::MissingTime(i as u64))?;
            start_times.push(start_time);
        }
        Ok(start_times)
    }

    fn extract_reads_relative_to_first_into(
        &mut self,
        from: Option<Duration>,
        to: Option<Duration>,
        writer: &mut Writer,
    ) -> Result<(usize, usize), IOError> {
        let header = self
            .read_header()
            .map_err(|source| IOError::ReadHeaderError {
                source: anyhow::Error::from(source),
            })?;
        let records = self.records(&header);
        let tag = Tag::new(b's', b't');
        let mut nb_reads_seen = 0;
        let mut nb_reads_written = 0;
        let mut anchor: Option<PrimitiveDateTime> = None;
        let mut earliest_bound = None;
        let mut latest_bound = None;

        writer
            .write_header(&header)
            .map_err(|source| IOError::WriteError {
                source: anyhow::Error::from(source),
            })?;

        for (i, record) in records.enumerate() {
            let record = record.map_err(|source| IOError::ParseAlignmentError {
                source: anyhow! { source.to_string() },
            })?;
            let data = record.data();
            let start_time = data
                .get(&tag)
                .ok_or(IOError::MissingTime(i as u64))?
                .map_err(|_| IOError::MissingTime(i as u64))?;
            let start_time = match start_time {
                noodles_sam::alignment::record::data::field::Value::String(s) => s,
                _ => return Err(IOError::MissingTime(i as u64)),
            };
            let start_time =
                parse_rfc3339_bytes(start_time).ok_or(IOError::MissingTime(i as u64))?;
            nb_reads_seen += 1;

            if anchor.is_none() {
                anchor = Some(start_time);
                earliest_bound = from.and_then(|dur| start_time.checked_add(dur));
                latest_bound = to.and_then(|dur| start_time.checked_add(dur));
            }

            if latest_bound.map_or(false, |max| start_time > max) {
                break;
            }

            let keep = earliest_bound.map_or(true, |min| start_time >= min)
                && latest_bound.map_or(true, |max| start_time <= max);

            if keep {
                writer
                    .write_record(&header, &record)
                    .map_err(|source| IOError::WriteError {
                        source: anyhow::Error::from(source),
                    })?;
                nb_reads_written += 1;
            }
        }

        writer.finish(&header).map_err(|source| IOError::WriteError {
            source: anyhow::Error::from(source),
        })?;

        Ok((nb_reads_seen, nb_reads_written))
    }

    fn extract_reads_between_into(
        &mut self,
        earliest: Option<&PrimitiveDateTime>,
        latest: Option<&PrimitiveDateTime>,
        writer: &mut Writer,
    ) -> Result<(usize, usize), IOError> {
        let header = self
            .read_header()
            .map_err(|source| IOError::ReadHeaderError {
                source: anyhow::Error::from(source),
            })?;
        let records = self.records(&header);
        let tag = Tag::new(b's', b't');
        let mut nb_reads_seen = 0;
        let mut nb_reads_written = 0;

        writer
            .write_header(&header)
            .map_err(|source| IOError::WriteError {
                source: anyhow::Error::from(source),
            })?;

        for (i, record) in records.enumerate() {
            let record = record.map_err(|source| IOError::ParseAlignmentError {
                source: anyhow! { source.to_string() },
            })?;
            let data = record.data();
            let start_time = data
                .get(&tag)
                .ok_or(IOError::MissingTime(i as u64))?
                .map_err(|_| IOError::MissingTime(i as u64))?;
            let start_time = match start_time {
                noodles_sam::alignment::record::data::field::Value::String(s) => s,
                _ => return Err(IOError::MissingTime(i as u64)),
            };
            let start_time =
                parse_rfc3339_bytes(start_time).ok_or(IOError::MissingTime(i as u64))?;
            nb_reads_seen += 1;

            let keep = earliest.map_or(true, |min| &start_time >= min)
                && latest.map_or(true, |max| &start_time <= max);

            if keep {
                writer
                    .write_record(&header, &record)
                    .map_err(|source| IOError::WriteError {
                        source: anyhow::Error::from(source),
                    })?;
                nb_reads_written += 1;
            }
        }

        writer.finish(&header).map_err(|source| IOError::WriteError {
            source: anyhow::Error::from(source),
        })?;

        Ok((nb_reads_seen, nb_reads_written))
    }

    fn extract_reads_in_timeframe_into(
        &mut self,
        selection: &ReadSelection,
        writer: &mut Writer,
    ) -> Result<(), IOError> {
        let header = self
            .read_header()
            .map_err(|source| IOError::ReadHeaderError {
                source: anyhow::Error::from(source),
            })?;
        let records = self.records(&header);
        let mut nb_reads_written = 0;
        let nb_reads_keep = selection.keep_count();
        let mut next_sparse_idx = 0;

        writer
            .write_header(&header)
            .map_err(|source| IOError::WriteError {
                source: anyhow::Error::from(source),
            })?;

        for (i, record) in records.enumerate() {
            let record = record.map_err(|source| IOError::ParseAlignmentError {
                source: anyhow! { source.to_string() },
            })?;
            let keep = match selection {
                ReadSelection::Dense(mask) => mask[i],
                ReadSelection::Sparse(indices) => {
                    if next_sparse_idx < indices.len() && indices[next_sparse_idx] == i {
                        next_sparse_idx += 1;
                        true
                    } else {
                        false
                    }
                }
            };
            if keep {
                writer
                    .write_record(&header, &record)
                    .map_err(|source| IOError::WriteError {
                        source: anyhow::Error::from(source),
                    })?;
                nb_reads_written += 1;
            }
        }

        writer.finish(&header).map_err(|source| IOError::WriteError {
            source: anyhow::Error::from(source),
        })?;

        if nb_reads_written == nb_reads_keep {
            Ok(())
        } else {
            Err(IOError::IndicesNotFound)
        }
    }
}
