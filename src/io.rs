use crate::cli::CompressionExt;
use anyhow::anyhow;
use needletail::errors::ParseErrorKind::EmptyFile;
use needletail::parse_fastx_file;
use noodles_bam as bam;
use noodles_bgzf::{self as bgzf, VirtualPosition};
use noodles_sam::alignment::record::data::field::Tag;
use noodles_util::alignment::io::Writer;
use ontime::{parse_rfc3339_bytes, FastxRecordExt, ReadSelection};
use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Seek, SeekFrom, Write};
use std::num::NonZeroUsize;
use std::path::{Path, PathBuf};
use std::time::UNIX_EPOCH;
use thiserror::Error;
use time::Duration;
use time::{OffsetDateTime, PrimitiveDateTime};

const SIDECAR_MAGIC: [u8; 8] = *b"ONTIDX4\0";
const PREVIOUS_BINARY_SIDECAR_MAGIC: [u8; 8] = *b"ONTIDX3\0";
const LEGACY_SIDECAR_PREFIX: [u8; 8] = *b"# ontime";
const SIDECAR_SPAN_SIZE: usize = 4096;

#[derive(Debug, Clone, PartialEq, Eq)]
struct SidecarSpanSummary {
    start_idx: u64,
    end_idx: u64,
    min_time: PrimitiveDateTime,
    max_time: PrimitiveDateTime,
    start_virtual_offset: Option<u64>,
    end_virtual_offset: Option<u64>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct SidecarIndex {
    timestamps: Vec<PrimitiveDateTime>,
    spans: Vec<SidecarSpanSummary>,
}

#[derive(Debug)]
struct ActiveSpan {
    start_idx: u64,
    len: usize,
    min_time: PrimitiveDateTime,
    max_time: PrimitiveDateTime,
    start_virtual_offset: Option<u64>,
    end_virtual_offset: Option<u64>,
}

struct SidecarStreamWriter {
    writer: BufWriter<File>,
    temp_path: PathBuf,
    final_path: PathBuf,
    count: u64,
    spans: Vec<SidecarSpanSummary>,
    active_span: Option<ActiveSpan>,
}

pub struct AlignmentTimeIndex {
    timestamps: Vec<PrimitiveDateTime>,
    spans: Vec<SidecarSpanSummary>,
}

impl AlignmentTimeIndex {
    pub fn start_times(&self) -> &[PrimitiveDateTime] {
        &self.timestamps
    }

    pub fn first_last(&self) -> Option<(PrimitiveDateTime, PrimitiveDateTime)> {
        self.spans.iter().fold(None, |acc, span| {
            Some(match acc {
                None => (span.min_time, span.max_time),
                Some((first, last)) => (first.min(span.min_time), last.max(span.max_time)),
            })
        })
    }

    pub fn valid_selection(
        &self,
        earliest: &PrimitiveDateTime,
        latest: &PrimitiveDateTime,
    ) -> ReadSelection {
        let mut selected_indices = Vec::new();

        for span in &self.spans {
            if span.max_time < *earliest || span.min_time > *latest {
                continue;
            }

            let start = span.start_idx as usize;
            let end = span.end_idx as usize;

            if *earliest <= span.min_time && span.max_time <= *latest {
                selected_indices.extend(start..end);
                continue;
            }

            selected_indices.extend(
                self.timestamps[start..end]
                    .iter()
                    .enumerate()
                    .filter_map(|(offset, timestamp)| {
                        (earliest <= timestamp && timestamp <= latest).then_some(start + offset)
                    }),
            );
        }

        if selected_indices.len() * std::mem::size_of::<usize>() < self.timestamps.len() {
            ReadSelection::Sparse(selected_indices)
        } else {
            let mut to_keep = vec![false; self.timestamps.len()];
            for idx in selected_indices {
                to_keep[idx] = true;
            }
            ReadSelection::Dense(to_keep)
        }
    }

    pub fn bam_candidate_spans(
        &self,
        earliest: &PrimitiveDateTime,
        latest: &PrimitiveDateTime,
    ) -> Vec<BamCandidateSpan> {
        self.spans
            .iter()
            .filter_map(|span| {
                let (start, end) = match (span.start_virtual_offset, span.end_virtual_offset) {
                    (Some(start), Some(end)) => (start, end),
                    _ => return None,
                };

                if span.max_time < *earliest || span.min_time > *latest {
                    return None;
                }

                Some(BamCandidateSpan {
                    start: VirtualPosition::from(start),
                    end: VirtualPosition::from(end),
                    fully_selected: *earliest <= span.min_time && span.max_time <= *latest,
                })
            })
            .collect()
    }
}

pub struct BamCandidateSpan {
    start: VirtualPosition,
    end: VirtualPosition,
    fully_selected: bool,
}

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

fn sidecar_temp_path(path: &Path) -> PathBuf {
    let mut temp = sidecar_path(path).as_os_str().to_os_string();
    temp.push(".tmp");
    PathBuf::from(temp)
}

fn file_fingerprint(path: &Path) -> std::io::Result<(u64, u128)> {
    let metadata = std::fs::metadata(path)?;
    let modified = metadata
        .modified()?
        .duration_since(UNIX_EPOCH)
        .map_err(|err| std::io::Error::new(std::io::ErrorKind::Other, err))?;
    Ok((metadata.len(), modified.as_nanos()))
}

fn read_sidecar(path: &Path) -> Result<Option<SidecarIndex>, IOError> {
    let sidecar = sidecar_path(path);
    if !sidecar.exists() {
        return Ok(None);
    }

    let expected = file_fingerprint(path).map_err(|source| IOError::SidecarError { source })?;
    let mut reader = BufReader::new(
        File::open(sidecar).map_err(|source| IOError::SidecarError { source })?,
    );

    let mut magic = [0u8; SIDECAR_MAGIC.len()];
    match reader.read_exact(&mut magic) {
        Ok(()) => {}
        Err(err) if err.kind() == std::io::ErrorKind::UnexpectedEof => return Ok(None),
        Err(source) => return Err(IOError::SidecarError { source }),
    }

    if magic == LEGACY_SIDECAR_PREFIX || magic == PREVIOUS_BINARY_SIDECAR_MAGIC {
        return Ok(None);
    }

    if magic != SIDECAR_MAGIC {
        return Ok(None);
    }

    let Some(size) = read_u64(&mut reader)? else {
        return Ok(None);
    };
    let Some(mtime) = read_u128(&mut reader)? else {
        return Ok(None);
    };
    let Some(count) = read_u64(&mut reader)? else {
        return Ok(None);
    };
    let Some(span_count) = read_u64(&mut reader)? else {
        return Ok(None);
    };

    if (size, mtime) != expected {
        return Ok(None);
    }

    let mut timestamps = Vec::with_capacity(count as usize);
    for _ in 0..count {
        let Some(timestamp_nanos) = read_i128(&mut reader)? else {
            return Ok(None);
        };
        let Ok(timestamp) = OffsetDateTime::from_unix_timestamp_nanos(timestamp_nanos) else {
            return Ok(None);
        };
        timestamps.push(timestamp.date().with_time(timestamp.time()));
    }

    let mut spans = Vec::with_capacity(span_count as usize);
    for _ in 0..span_count {
        let Some(start_idx) = read_u64(&mut reader)? else {
            return Ok(None);
        };
        let Some(end_idx) = read_u64(&mut reader)? else {
            return Ok(None);
        };
        let Some(min_time) = read_timestamp(&mut reader)? else {
            return Ok(None);
        };
        let Some(max_time) = read_timestamp(&mut reader)? else {
            return Ok(None);
        };
        let Some(start_virtual_offset) = read_u64(&mut reader)? else {
            return Ok(None);
        };
        let Some(end_virtual_offset) = read_u64(&mut reader)? else {
            return Ok(None);
        };
        spans.push(SidecarSpanSummary {
            start_idx,
            end_idx,
            min_time,
            max_time,
            start_virtual_offset: (start_virtual_offset != u64::MAX).then_some(start_virtual_offset),
            end_virtual_offset: (end_virtual_offset != u64::MAX).then_some(end_virtual_offset),
        });
    }

    Ok(Some(SidecarIndex { timestamps, spans }))
}

#[cfg(test)]
fn write_sidecar(path: &Path, index: &SidecarIndex) -> Result<(), IOError> {
    let sidecar = sidecar_path(path);
    let (size, mtime_nanos) =
        file_fingerprint(path).map_err(|source| IOError::SidecarError { source })?;
    let file = File::create(sidecar).map_err(|source| IOError::SidecarError { source })?;
    let mut writer = BufWriter::new(file);

    writer
        .write_all(&SIDECAR_MAGIC)
        .map_err(|source| IOError::SidecarError { source })?;
    writer
        .write_all(&size.to_le_bytes())
        .map_err(|source| IOError::SidecarError { source })?;
    writer
        .write_all(&mtime_nanos.to_le_bytes())
        .map_err(|source| IOError::SidecarError { source })?;
    writer
        .write_all(&(index.timestamps.len() as u64).to_le_bytes())
        .map_err(|source| IOError::SidecarError { source })?;
    writer
        .write_all(&(index.spans.len() as u64).to_le_bytes())
        .map_err(|source| IOError::SidecarError { source })?;
    for timestamp in &index.timestamps {
        let timestamp_nanos = timestamp_to_nanos(*timestamp);
        writer
            .write_all(&timestamp_nanos.to_le_bytes())
            .map_err(|source| IOError::SidecarError { source })?;
    }
    for span in &index.spans {
        writer
            .write_all(&span.start_idx.to_le_bytes())
            .map_err(|source| IOError::SidecarError { source })?;
        writer
            .write_all(&span.end_idx.to_le_bytes())
            .map_err(|source| IOError::SidecarError { source })?;
        writer
            .write_all(&timestamp_to_nanos(span.min_time).to_le_bytes())
            .map_err(|source| IOError::SidecarError { source })?;
        writer
            .write_all(&timestamp_to_nanos(span.max_time).to_le_bytes())
            .map_err(|source| IOError::SidecarError { source })?;
        writer
            .write_all(&span.start_virtual_offset.unwrap_or(u64::MAX).to_le_bytes())
            .map_err(|source| IOError::SidecarError { source })?;
        writer
            .write_all(&span.end_virtual_offset.unwrap_or(u64::MAX).to_le_bytes())
            .map_err(|source| IOError::SidecarError { source })?;
    }

    Ok(())
}

impl SidecarStreamWriter {
    fn new(path: &Path) -> Result<Self, IOError> {
        let final_path = sidecar_path(path);
        let temp_path = sidecar_temp_path(path);
        let (size, mtime_nanos) =
            file_fingerprint(path).map_err(|source| IOError::SidecarError { source })?;
        let file = File::create(&temp_path).map_err(|source| IOError::SidecarError { source })?;
        let mut writer = BufWriter::new(file);

        writer
            .write_all(&SIDECAR_MAGIC)
            .map_err(|source| IOError::SidecarError { source })?;
        writer
            .write_all(&size.to_le_bytes())
            .map_err(|source| IOError::SidecarError { source })?;
        writer
            .write_all(&mtime_nanos.to_le_bytes())
            .map_err(|source| IOError::SidecarError { source })?;
        writer
            .write_all(&0u64.to_le_bytes())
            .map_err(|source| IOError::SidecarError { source })?;
        writer
            .write_all(&0u64.to_le_bytes())
            .map_err(|source| IOError::SidecarError { source })?;

        Ok(Self {
            writer,
            temp_path,
            final_path,
            count: 0,
            spans: Vec::new(),
            active_span: None,
        })
    }

    fn push_timestamp(
        &mut self,
        timestamp: PrimitiveDateTime,
        start_virtual_offset: Option<u64>,
        end_virtual_offset: Option<u64>,
    ) -> Result<(), IOError> {
        self.writer
            .write_all(&timestamp_to_nanos(timestamp).to_le_bytes())
            .map_err(|source| IOError::SidecarError { source })?;

        match &mut self.active_span {
            Some(span) if span.len < SIDECAR_SPAN_SIZE => {
                span.len += 1;
                span.min_time = span.min_time.min(timestamp);
                span.max_time = span.max_time.max(timestamp);
                span.end_virtual_offset = end_virtual_offset;
            }
            Some(_) => {
                self.flush_active_span();
                self.active_span = Some(ActiveSpan {
                    start_idx: self.count,
                    len: 1,
                    min_time: timestamp,
                    max_time: timestamp,
                    start_virtual_offset,
                    end_virtual_offset,
                });
            }
            None => {
                self.active_span = Some(ActiveSpan {
                    start_idx: self.count,
                    len: 1,
                    min_time: timestamp,
                    max_time: timestamp,
                    start_virtual_offset,
                    end_virtual_offset,
                });
            }
        }

        self.count += 1;
        Ok(())
    }

    fn finish(mut self) -> Result<Vec<SidecarSpanSummary>, IOError> {
        self.flush_active_span();
        self.writer.flush().map_err(|source| IOError::SidecarError { source })?;

        let mut file = self
            .writer
            .into_inner()
            .map_err(|err| IOError::SidecarError {
                source: err.into_error(),
            })?;

        for span in &self.spans {
            file.write_all(&span.start_idx.to_le_bytes())
                .map_err(|source| IOError::SidecarError { source })?;
            file.write_all(&span.end_idx.to_le_bytes())
                .map_err(|source| IOError::SidecarError { source })?;
            file.write_all(&timestamp_to_nanos(span.min_time).to_le_bytes())
                .map_err(|source| IOError::SidecarError { source })?;
            file.write_all(&timestamp_to_nanos(span.max_time).to_le_bytes())
                .map_err(|source| IOError::SidecarError { source })?;
            file.write_all(&span.start_virtual_offset.unwrap_or(u64::MAX).to_le_bytes())
                .map_err(|source| IOError::SidecarError { source })?;
            file.write_all(&span.end_virtual_offset.unwrap_or(u64::MAX).to_le_bytes())
                .map_err(|source| IOError::SidecarError { source })?;
        }

        let count_offset = (SIDECAR_MAGIC.len() + std::mem::size_of::<u64>() + std::mem::size_of::<u128>()) as u64;
        file.seek(SeekFrom::Start(count_offset))
            .map_err(|source| IOError::SidecarError { source })?;
        file.write_all(&self.count.to_le_bytes())
            .map_err(|source| IOError::SidecarError { source })?;
        file.write_all(&(self.spans.len() as u64).to_le_bytes())
            .map_err(|source| IOError::SidecarError { source })?;
        file.flush().map_err(|source| IOError::SidecarError { source })?;

        std::fs::rename(&self.temp_path, &self.final_path)
            .map_err(|source| IOError::SidecarError { source })?;

        Ok(self.spans)
    }

    fn flush_active_span(&mut self) {
        if let Some(span) = self.active_span.take() {
            self.spans.push(SidecarSpanSummary {
                start_idx: span.start_idx,
                end_idx: span.start_idx + span.len as u64,
                min_time: span.min_time,
                max_time: span.max_time,
                start_virtual_offset: span.start_virtual_offset,
                end_virtual_offset: span.end_virtual_offset,
            });
        }
    }
}

fn read_u64<R: Read>(reader: &mut R) -> Result<Option<u64>, IOError> {
    let mut buf = [0u8; std::mem::size_of::<u64>()];
    match reader.read_exact(&mut buf) {
        Ok(()) => Ok(Some(u64::from_le_bytes(buf))),
        Err(err) if err.kind() == std::io::ErrorKind::UnexpectedEof => Ok(None),
        Err(source) => Err(IOError::SidecarError { source }),
    }
}

fn read_u128<R: Read>(reader: &mut R) -> Result<Option<u128>, IOError> {
    let mut buf = [0u8; std::mem::size_of::<u128>()];
    match reader.read_exact(&mut buf) {
        Ok(()) => Ok(Some(u128::from_le_bytes(buf))),
        Err(err) if err.kind() == std::io::ErrorKind::UnexpectedEof => Ok(None),
        Err(source) => Err(IOError::SidecarError { source }),
    }
}

fn read_i128<R: Read>(reader: &mut R) -> Result<Option<i128>, IOError> {
    let mut buf = [0u8; std::mem::size_of::<i128>()];
    match reader.read_exact(&mut buf) {
        Ok(()) => Ok(Some(i128::from_le_bytes(buf))),
        Err(err) if err.kind() == std::io::ErrorKind::UnexpectedEof => Ok(None),
        Err(source) => Err(IOError::SidecarError { source }),
    }
}

fn read_timestamp<R: Read>(reader: &mut R) -> Result<Option<PrimitiveDateTime>, IOError> {
    let Some(timestamp_nanos) = read_i128(reader)? else {
        return Ok(None);
    };
    let Ok(timestamp) = OffsetDateTime::from_unix_timestamp_nanos(timestamp_nanos) else {
        return Ok(None);
    };
    Ok(Some(timestamp.date().with_time(timestamp.time())))
}

fn timestamp_to_nanos(timestamp: PrimitiveDateTime) -> i128 {
    timestamp.assume_utc().unix_timestamp_nanos()
}

fn build_sidecar_index(timestamps: Vec<PrimitiveDateTime>) -> SidecarIndex {
    let spans = timestamps
        .chunks(SIDECAR_SPAN_SIZE)
        .enumerate()
        .filter_map(|(chunk_idx, chunk)| {
            let (&first, rest) = chunk.split_first()?;
            let (min_time, max_time) = rest.iter().copied().fold((first, first), |acc, timestamp| {
                (acc.0.min(timestamp), acc.1.max(timestamp))
            });
            let start_idx = (chunk_idx * SIDECAR_SPAN_SIZE) as u64;
            let end_idx = start_idx + chunk.len() as u64;
            Some(SidecarSpanSummary {
                start_idx,
                end_idx,
                min_time,
                max_time,
                start_virtual_offset: None,
                end_virtual_offset: None,
            })
        })
        .collect();

    SidecarIndex { timestamps, spans }
}

fn alignment_sidecar_index_from_path(path: &Path) -> Result<SidecarIndex, IOError> {
    if let Some(index) = read_sidecar(path)? {
        return Ok(index);
    }

    match path.extension().and_then(|ext| ext.to_str()) {
        Some("bam") => bam_time_index_from_path_parallel(path),
        _ => sam_time_index_from_path(path),
    }
}

pub fn alignment_time_index_from_path(path: &Path) -> Result<AlignmentTimeIndex, IOError> {
    alignment_sidecar_index_from_path(path).map(|index| AlignmentTimeIndex {
        timestamps: index.timestamps,
        spans: index.spans,
    })
}

fn sam_time_index_from_path(path: &Path) -> Result<SidecarIndex, IOError> {
    let mut reader = noodles_util::alignment::io::reader::Builder::default()
        .build_from_path(path)
        .map_err(|source| IOError::ReadHeaderError {
            source: anyhow::Error::from(source),
        })?;
    let header = reader
        .read_header()
        .map_err(|source| IOError::ReadHeaderError {
            source: anyhow::Error::from(source),
        })?;
    let records = reader.records(&header);
    let tag = Tag::new(b's', b't');
    let mut timestamps = Vec::new();
    let mut sidecar_writer = SidecarStreamWriter::new(path).ok();

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
        if let Some(writer) = sidecar_writer.as_mut() {
            if writer.push_timestamp(start_time, None, None).is_err() {
                sidecar_writer = None;
            }
        }
        timestamps.push(start_time);
    }

    let spans = match sidecar_writer {
        Some(writer) => writer.finish().unwrap_or_else(|_| build_sidecar_index(timestamps.clone()).spans),
        None => build_sidecar_index(timestamps.clone()).spans,
    };

    Ok(SidecarIndex { timestamps, spans })
}

fn bam_time_index_from_path_parallel(path: &Path) -> Result<SidecarIndex, IOError> {
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
    let mut timestamps = Vec::new();
    let mut sidecar_writer = SidecarStreamWriter::new(path).ok();
    let mut record = bam::Record::default();
    let mut i = 0u64;

    loop {
        let start_virtual_offset = u64::from(reader.get_ref().virtual_position());
        let block_size = reader
            .read_record(&mut record)
            .map_err(|source| IOError::ParseAlignmentError {
                source: anyhow! { source.to_string() },
            })?;
        if block_size == 0 {
            break;
        }
        let end_virtual_offset = u64::from(reader.get_ref().virtual_position());
        let data = record.data();
        let start_time = data
            .get(&tag)
            .ok_or(IOError::MissingTime(i))?
            .map_err(|_| IOError::MissingTime(i))?;
        let start_time = match start_time {
            noodles_sam::alignment::record::data::field::Value::String(s) => s,
            _ => return Err(IOError::MissingTime(i)),
        };
        let start_time = parse_rfc3339_bytes(start_time).ok_or(IOError::MissingTime(i))?;
        if let Some(writer) = sidecar_writer.as_mut() {
            if writer
                .push_timestamp(
                    start_time,
                    Some(start_virtual_offset),
                    Some(end_virtual_offset),
                )
                .is_err()
            {
                sidecar_writer = None;
            }
        }
        timestamps.push(start_time);
        i += 1;
    }

    let spans = match sidecar_writer {
        Some(writer) => writer.finish().unwrap_or_else(|_| build_sidecar_index(timestamps.clone()).spans),
        None => build_sidecar_index(timestamps.clone()).spans,
    };

    Ok(SidecarIndex { timestamps, spans })
}

pub fn extract_bam_candidate_spans_into<W: Write>(
    path: &Path,
    spans: &[BamCandidateSpan],
    earliest: &PrimitiveDateTime,
    latest: &PrimitiveDateTime,
    writer: &mut Writer<W>,
) -> Result<usize, IOError> {
    let file = File::open(path).map_err(|source| IOError::ReadHeaderError {
        source: anyhow::Error::from(source),
    })?;
    let mut reader = bam::io::Reader::new(file);
    let header = reader
        .read_header()
        .map_err(|source| IOError::ReadHeaderError {
            source: anyhow::Error::from(source),
        })?;
    let tag = Tag::new(b's', b't');
    let mut record = bam::Record::default();
    let mut nb_reads_written = 0;

    writer
        .write_header(&header)
        .map_err(|source| IOError::WriteError {
            source: anyhow::Error::from(source),
        })?;

    for span in spans {
        bgzf::io::Seek::seek_to_virtual_position(reader.get_mut(), span.start).map_err(
            |source| IOError::ReadHeaderError {
                source: anyhow::Error::from(source),
            },
        )?;

        while reader.get_ref().virtual_position() < span.end {
            let block_size = reader
                .read_record(&mut record)
                .map_err(|source| IOError::ParseAlignmentError {
                    source: anyhow! { source.to_string() },
                })?;

            if block_size == 0 {
                break;
            }

            let keep = if span.fully_selected {
                true
            } else {
                let data = record.data();
                let start_time = data
                    .get(&tag)
                    .ok_or(IOError::MissingTime(0))?
                    .map_err(|_| IOError::MissingTime(0))?;
                let start_time = match start_time {
                    noodles_sam::alignment::record::data::field::Value::String(s) => s,
                    _ => return Err(IOError::MissingTime(0)),
                };
                let start_time = parse_rfc3339_bytes(start_time).ok_or(IOError::MissingTime(0))?;
                earliest <= &start_time && &start_time <= latest
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
    }

    writer.finish(&header).map_err(|source| IOError::WriteError {
        source: anyhow::Error::from(source),
    })?;

    Ok(nb_reads_written)
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
    fn extract_reads_relative_to_first_into<W: Write>(
        &mut self,
        from: Option<Duration>,
        to: Option<Duration>,
        writer: &mut Writer<W>,
    ) -> Result<(usize, usize), IOError>;
    fn extract_reads_between_into<W: Write>(
        &mut self,
        earliest: Option<&PrimitiveDateTime>,
        latest: Option<&PrimitiveDateTime>,
        writer: &mut Writer<W>,
    ) -> Result<(usize, usize), IOError>;
    fn extract_reads_in_timeframe_into<W: Write>(
        &mut self,
        selection: &ReadSelection,
        writer: &mut Writer<W>,
    ) -> Result<(), IOError>;
}

impl<R> TimeExt for noodles_util::alignment::io::reader::Reader<R>
where
    R: Read,
{
    fn extract_reads_relative_to_first_into<W: Write>(
        &mut self,
        from: Option<Duration>,
        to: Option<Duration>,
        writer: &mut Writer<W>,
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

    fn extract_reads_between_into<W: Write>(
        &mut self,
        earliest: Option<&PrimitiveDateTime>,
        latest: Option<&PrimitiveDateTime>,
        writer: &mut Writer<W>,
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

    fn extract_reads_in_timeframe_into<W: Write>(
        &mut self,
        selection: &ReadSelection,
        writer: &mut Writer<W>,
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

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::NamedTempFile;
    use time::macros::{date, time};

    #[test]
    fn build_sidecar_index_adds_span_summaries() {
        let timestamps: Vec<_> = (0..(SIDECAR_SPAN_SIZE + 2))
            .map(|seconds| {
                PrimitiveDateTime::new(
                    date!(2023 - 09 - 23),
                    time!(16:00:00) + Duration::seconds(seconds as i64),
                )
            })
            .collect();

        let index = build_sidecar_index(timestamps.clone());

        assert_eq!(index.timestamps, timestamps);
        assert_eq!(index.spans.len(), 2);
        assert_eq!(index.spans[0].start_idx, 0);
        assert_eq!(index.spans[0].end_idx, SIDECAR_SPAN_SIZE as u64);
        assert_eq!(index.spans[1].start_idx, SIDECAR_SPAN_SIZE as u64);
        assert_eq!(index.spans[1].end_idx, (SIDECAR_SPAN_SIZE + 2) as u64);
        assert_eq!(index.spans[0].min_time, timestamps[0]);
        assert_eq!(index.spans[0].max_time, timestamps[SIDECAR_SPAN_SIZE - 1]);
        assert_eq!(index.spans[1].min_time, timestamps[SIDECAR_SPAN_SIZE]);
        assert_eq!(index.spans[1].max_time, timestamps[SIDECAR_SPAN_SIZE + 1]);
    }

    #[test]
    fn binary_sidecar_roundtrips_with_spans() {
        let mut input = NamedTempFile::new().unwrap();
        input.write_all(b"alignment input").unwrap();

        let timestamps = vec![
            PrimitiveDateTime::new(date!(2023 - 09 - 23), time!(16:02:50.388)),
            PrimitiveDateTime::new(date!(2023 - 09 - 23), time!(18:10:34.160)),
            PrimitiveDateTime::new(date!(2023 - 09 - 23), time!(18:20:00)),
        ];
        let expected = build_sidecar_index(timestamps);

        write_sidecar(input.path(), &expected).unwrap();
        let actual = read_sidecar(input.path()).unwrap().unwrap();

        assert_eq!(actual, expected);
    }

    #[test]
    fn streaming_sidecar_writer_roundtrips_with_spans() {
        let mut input = NamedTempFile::new().unwrap();
        input.write_all(b"alignment input").unwrap();

        let timestamps = vec![
            PrimitiveDateTime::new(date!(2023 - 09 - 23), time!(16:02:50.388)),
            PrimitiveDateTime::new(date!(2023 - 09 - 23), time!(18:10:34.160)),
            PrimitiveDateTime::new(date!(2023 - 09 - 23), time!(18:20:00)),
        ];

        let mut writer = SidecarStreamWriter::new(input.path()).unwrap();
        for timestamp in &timestamps {
            writer.push_timestamp(*timestamp, None, None).unwrap();
        }
        let spans = writer.finish().unwrap();
        let actual = read_sidecar(input.path()).unwrap().unwrap();

        assert_eq!(actual.timestamps, timestamps);
        assert_eq!(actual.spans, spans);
    }

    #[test]
    fn alignment_time_index_valid_selection_prunes_by_span() {
        let timestamps: Vec<_> = (0..(SIDECAR_SPAN_SIZE * 3))
            .map(|seconds| {
                PrimitiveDateTime::new(
                    date!(2023 - 09 - 23),
                    time!(16:00:00) + Duration::seconds(seconds as i64),
                )
            })
            .collect();
        let index = AlignmentTimeIndex {
            timestamps: timestamps.clone(),
            spans: build_sidecar_index(timestamps).spans,
        };

        let earliest = PrimitiveDateTime::new(
            date!(2023 - 09 - 23),
            time!(16:00:00) + Duration::seconds(SIDECAR_SPAN_SIZE as i64),
        );
        let latest = PrimitiveDateTime::new(
            date!(2023 - 09 - 23),
            time!(16:00:00) + Duration::seconds((SIDECAR_SPAN_SIZE * 2 - 1) as i64),
        );

        let selection = index.valid_selection(&earliest, &latest);

        assert_eq!(selection.keep_count(), SIDECAR_SPAN_SIZE);
        match selection {
            ReadSelection::Sparse(indices) => {
                assert_eq!(indices.first().copied(), Some(SIDECAR_SPAN_SIZE));
                assert_eq!(indices.last().copied(), Some(SIDECAR_SPAN_SIZE * 2 - 1));
            }
            ReadSelection::Dense(mask) => {
                assert_eq!(mask.iter().filter(|keep| **keep).count(), SIDECAR_SPAN_SIZE);
                assert!(!mask[SIDECAR_SPAN_SIZE - 1]);
                assert!(mask[SIDECAR_SPAN_SIZE]);
                assert!(mask[SIDECAR_SPAN_SIZE * 2 - 1]);
                assert!(!mask[SIDECAR_SPAN_SIZE * 2]);
            }
        }
    }

    #[test]
    fn bam_candidate_spans_require_offsets_and_prune_outside_ranges() {
        let timestamps = vec![
            PrimitiveDateTime::new(date!(2023 - 09 - 23), time!(16:00:00)),
            PrimitiveDateTime::new(date!(2023 - 09 - 23), time!(16:00:01)),
            PrimitiveDateTime::new(date!(2023 - 09 - 23), time!(16:00:02)),
        ];
        let index = AlignmentTimeIndex {
            timestamps,
            spans: vec![
                SidecarSpanSummary {
                    start_idx: 0,
                    end_idx: 1,
                    min_time: PrimitiveDateTime::new(date!(2023 - 09 - 23), time!(15:59:00)),
                    max_time: PrimitiveDateTime::new(date!(2023 - 09 - 23), time!(15:59:59)),
                    start_virtual_offset: Some(10),
                    end_virtual_offset: Some(20),
                },
                SidecarSpanSummary {
                    start_idx: 1,
                    end_idx: 2,
                    min_time: PrimitiveDateTime::new(date!(2023 - 09 - 23), time!(16:00:00)),
                    max_time: PrimitiveDateTime::new(date!(2023 - 09 - 23), time!(16:00:10)),
                    start_virtual_offset: Some(20),
                    end_virtual_offset: Some(30),
                },
                SidecarSpanSummary {
                    start_idx: 2,
                    end_idx: 3,
                    min_time: PrimitiveDateTime::new(date!(2023 - 09 - 23), time!(16:00:05)),
                    max_time: PrimitiveDateTime::new(date!(2023 - 09 - 23), time!(16:00:06)),
                    start_virtual_offset: None,
                    end_virtual_offset: None,
                },
            ],
        };

        let spans = index.bam_candidate_spans(
            &PrimitiveDateTime::new(date!(2023 - 09 - 23), time!(16:00:00)),
            &PrimitiveDateTime::new(date!(2023 - 09 - 23), time!(16:00:10)),
        );

        assert_eq!(spans.len(), 1);
        assert_eq!(u64::from(spans[0].start), 20);
        assert_eq!(u64::from(spans[0].end), 30);
        assert!(spans[0].fully_selected);
    }
}
