use duration_str::DError;
use needletail::parser::SequenceRecord;
use time::format_description::well_known::Rfc3339;
use time::{Duration, PrimitiveDateTime};

pub trait FastxRecordExt {
    fn start_time(&self) -> Option<PrimitiveDateTime>;
}

impl FastxRecordExt for SequenceRecord<'_> {
    fn start_time(&self) -> Option<PrimitiveDateTime> {
        extract_embedded_timestamp(self.id())
    }
}

pub trait DurationExt {
    fn from_str(s: &str) -> Result<Self, DError>
    where
        Self: Sized;
}

impl DurationExt for Duration {
    fn from_str(s: &str) -> Result<Self, DError> {
        if let Some(pos_s) = s.strip_prefix('-') {
            let dur = duration_str::parse_time(pos_s)?;
            Ok(-1 * dur)
        } else {
            duration_str::parse_time(s)
        }
    }
}

pub enum ReadSelection {
    Dense(Vec<bool>),
    Sparse(Vec<usize>),
}

impl ReadSelection {
    pub fn keep_count(&self) -> usize {
        match self {
            Self::Dense(mask) => mask.iter().filter(|keep| **keep).count(),
            Self::Sparse(indices) => indices.len(),
        }
    }
}

pub fn parse_rfc3339_bytes(bytes: &[u8]) -> Option<PrimitiveDateTime> {
    let timestamp = std::str::from_utf8(bytes).ok()?;
    PrimitiveDateTime::parse(timestamp, &Rfc3339).ok()
}

pub fn extract_embedded_timestamp(bytes: &[u8]) -> Option<PrimitiveDateTime> {
    const FASTQ_PREFIX: &[u8] = b"start_time=";
    const BAM_PREFIX: &[u8] = b"st:Z:";

    let start = find_subslice(bytes, FASTQ_PREFIX)
        .map(|idx| idx + FASTQ_PREFIX.len())
        .or_else(|| find_subslice(bytes, BAM_PREFIX).map(|idx| idx + BAM_PREFIX.len()))?;

    let end = bytes[start..]
        .iter()
        .position(|b| b.is_ascii_whitespace())
        .map_or(bytes.len(), |offset| start + offset);

    parse_rfc3339_bytes(&bytes[start..end])
}

fn find_subslice(haystack: &[u8], needle: &[u8]) -> Option<usize> {
    haystack.windows(needle.len()).position(|window| window == needle)
}

pub fn valid_selection(
    timestamps: &[PrimitiveDateTime],
    earliest: &PrimitiveDateTime,
    latest: &PrimitiveDateTime,
) -> ReadSelection {
    let selected_indices: Vec<usize> = timestamps
        .iter()
        .enumerate()
        .filter_map(|(i, t)| (earliest <= t && t <= latest).then_some(i))
        .collect();

    if selected_indices.len() * std::mem::size_of::<usize>() < timestamps.len() {
        ReadSelection::Sparse(selected_indices)
    } else {
        let mut to_keep: Vec<bool> = vec![false; timestamps.len()];
        for idx in selected_indices {
            to_keep[idx] = true;
        }
        ReadSelection::Dense(to_keep)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use needletail::parse_fastx_file;
    use std::io::Write;
    use tempfile::Builder;
    use time::macros::{date, time};
    use time::Duration;

    #[test]
    fn test_no_start_time() {
        let text = "@read1\nA\n+\n1";
        let mut file = Builder::new().suffix(".fa").tempfile().unwrap();
        file.write_all(text.as_bytes()).unwrap();

        let mut reader = parse_fastx_file(file.path()).unwrap();
        let rec = reader.next().unwrap();
        let record = rec.unwrap();

        let actual = record.start_time();
        let expected = None;

        assert_eq!(actual, expected)
    }

    #[test]
    fn test_start_time_old_valid() {
        let text = "@read1 ch=352 start_time=2022-12-12T18:39:27Z model_version_id=2021\nA\n+\n1";
        let mut file = Builder::new().suffix(".fa").tempfile().unwrap();
        file.write_all(text.as_bytes()).unwrap();

        let mut reader = parse_fastx_file(file.path()).unwrap();
        let rec = reader.next().unwrap();
        let record = rec.unwrap();

        let actual = record.start_time().unwrap();
        let expected = PrimitiveDateTime::new(date!(2022 - 12 - 12), time!(18:39:27));

        assert_eq!(actual, expected)
    }

    #[test]
    fn test_start_time_offset_valid() {
        let text =
            "@read1 ch=352 start_time=2021-07-08T17:47:25+01:00 model_version_id=2021\nA\n+\n1";
        let mut file = Builder::new().suffix(".fa").tempfile().unwrap();
        file.write_all(text.as_bytes()).unwrap();

        let mut reader = parse_fastx_file(file.path()).unwrap();
        let rec = reader.next().unwrap();
        let record = rec.unwrap();

        let actual = record.start_time().unwrap();
        let expected = PrimitiveDateTime::new(date!(2021 - 07 - 08), time!(17:47:25));

        assert_eq!(actual, expected)
    }

    #[test]
    fn test_start_time_offset_with_micro_valid() {
        let text = "@read1 ch=352 start_time=2021-07-08T17:47:25.558027+01:00 model_version_id=2021\nA\n+\n1";
        let mut file = Builder::new().suffix(".fa").tempfile().unwrap();
        file.write_all(text.as_bytes()).unwrap();

        let mut reader = parse_fastx_file(file.path()).unwrap();
        let rec = reader.next().unwrap();
        let record = rec.unwrap();

        let actual = record.start_time().unwrap();
        let expected = PrimitiveDateTime::new(date!(2021 - 07 - 08), time!(17:47:25.558027));

        assert_eq!(actual, expected)
    }

    #[test]
    fn test_start_time_invalid_without_z() {
        let text = "@read1 ch=352 start_time=2022-12-12T18:39:27 model_version_id=2021\nA\n+\n1";
        let mut file = Builder::new().suffix(".fa").tempfile().unwrap();
        file.write_all(text.as_bytes()).unwrap();

        let mut reader = parse_fastx_file(file.path()).unwrap();
        let rec = reader.next().unwrap();
        let record = rec.unwrap();

        let actual = record.start_time();
        assert!(actual.is_none())
    }

    #[test]
    fn test_start_time_invalid() {
        let text = "@read1 ch=352 start_time=2022-12-12T18:39Z model_version_id=2021\nA\n+\n1";
        let mut file = Builder::new().suffix(".fa").tempfile().unwrap();
        file.write_all(text.as_bytes()).unwrap();

        let mut reader = parse_fastx_file(file.path()).unwrap();
        let rec = reader.next().unwrap();
        let record = rec.unwrap();

        let actual = record.start_time();
        let expected = None;

        assert_eq!(actual, expected)
    }

    #[test]
    fn test_bam_tag_start_time_is_valid() {
        let text = "@read1 st:Z:2023-08-07T13:14:42.356+00:00\nA\n+\n1";
        let mut file = Builder::new().suffix(".fa").tempfile().unwrap();
        file.write_all(text.as_bytes()).unwrap();

        let mut reader = parse_fastx_file(file.path()).unwrap();
        let rec = reader.next().unwrap();
        let record = rec.unwrap();

        let actual = record.start_time().unwrap();
        let expected = PrimitiveDateTime::new(date!(2023 - 08 - 07), time!(13:14:42.356));

        assert_eq!(actual, expected)
    }

    #[test]
    fn test_duration_from_str_negative() {
        let s = "-1h";
        let actual = Duration::from_str(s).unwrap();
        let expected = Duration::hours(-1);

        assert_eq!(actual, expected)
    }

    #[test]
    fn test_duration_from_str_negative_invalid() {
        let s = "1d-1h";
        let actual = Duration::from_str(s);
        assert!(actual.is_err())
    }

    #[test]
    fn test_duration_from_str() {
        let s = "11h30min";
        let actual = Duration::from_str(s).unwrap();
        let expected = Duration::seconds(41_400);

        assert_eq!(actual, expected)
    }

    #[test]
    fn test_duration_from_str_invalid() {
        let s = "11h30min12foo";
        let actual = Duration::from_str(s);
        assert!(actual.is_err())
    }

    #[test]
    fn test_valid_selection_prefers_sparse_for_sparse_windows() {
        let timestamps = vec![
            PrimitiveDateTime::new(date!(2022 - 12 - 12), time!(12:00:00)),
            PrimitiveDateTime::new(date!(2022 - 12 - 12), time!(13:00:00)),
            PrimitiveDateTime::new(date!(2022 - 12 - 12), time!(14:00:00)),
            PrimitiveDateTime::new(date!(2022 - 12 - 12), time!(15:00:00)),
            PrimitiveDateTime::new(date!(2022 - 12 - 12), time!(16:00:00)),
            PrimitiveDateTime::new(date!(2022 - 12 - 12), time!(17:00:00)),
            PrimitiveDateTime::new(date!(2022 - 12 - 12), time!(18:00:00)),
            PrimitiveDateTime::new(date!(2022 - 12 - 12), time!(19:00:00)),
            PrimitiveDateTime::new(date!(2022 - 12 - 12), time!(20:00:00)),
        ];

        let selection = valid_selection(
            &timestamps,
            &PrimitiveDateTime::new(date!(2022 - 12 - 12), time!(12:00:00)),
            &PrimitiveDateTime::new(date!(2022 - 12 - 12), time!(12:00:00)),
        );

        match selection {
            ReadSelection::Sparse(indices) => assert_eq!(indices, vec![0]),
            ReadSelection::Dense(_) => panic!("expected sparse selection"),
        }
    }

    #[test]
    fn test_valid_selection_prefers_dense_for_dense_windows() {
        let timestamps = vec![
            PrimitiveDateTime::new(date!(2022 - 12 - 12), time!(12:00:00)),
            PrimitiveDateTime::new(date!(2022 - 12 - 12), time!(13:00:00)),
            PrimitiveDateTime::new(date!(2022 - 12 - 12), time!(14:00:00)),
            PrimitiveDateTime::new(date!(2022 - 12 - 12), time!(15:00:00)),
        ];

        let selection = valid_selection(
            &timestamps,
            &PrimitiveDateTime::new(date!(2022 - 12 - 12), time!(12:00:00)),
            &PrimitiveDateTime::new(date!(2022 - 12 - 12), time!(15:00:00)),
        );

        match selection {
            ReadSelection::Dense(mask) => assert_eq!(mask, vec![true, true, true, true]),
            ReadSelection::Sparse(_) => panic!("expected dense selection"),
        }
    }
}
