#!/usr/bin/env bash

set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BIN="${ROOT_DIR}/target/release/ontime"
WORK_DIR="${TMPDIR:-/tmp}/ontime-bench"
SOURCE_SAM="${ROOT_DIR}/tests/cases/test.sam"
RECORD_GROUPS="${RECORD_GROUPS:-20}"
READS_PER_GROUP="${READS_PER_GROUP:-1000}"

require_cmd() {
  if ! command -v "$1" >/dev/null 2>&1; then
    echo "missing required command: $1" >&2
    exit 1
  fi
}

time_run() {
  local label="$1"
  shift
  local start_ns end_ns elapsed_ms
  start_ns="$(date +%s%N)"
  "$@" >/dev/null
  end_ns="$(date +%s%N)"
  elapsed_ms="$(( (end_ns - start_ns) / 1000000 ))"
  printf '%s\t%s ms\n' "$label" "$elapsed_ms"
}

make_clustered_sam() {
  local output="$1"
  python3 - "$SOURCE_SAM" "$output" "$RECORD_GROUPS" "$READS_PER_GROUP" <<'PY'
from datetime import datetime, timedelta, timezone
from pathlib import Path
import sys

source = Path(sys.argv[1])
output = Path(sys.argv[2])
groups = int(sys.argv[3])
reads_per_group = int(sys.argv[4])

lines = source.read_text().splitlines()
header = [line for line in lines if line.startswith("@")]
template = next(line for line in lines if not line.startswith("@")).split("\t")
base = datetime(2023, 9, 22, 0, 0, 0, tzinfo=timezone.utc)

with output.open("w") as handle:
    for line in header:
        handle.write(line + "\n")

    for group in range(groups):
        ts = (base + timedelta(hours=group)).isoformat().replace("+00:00", "Z")
        for i in range(reads_per_group):
            fields = template.copy()
            fields[0] = f"read_{group}_{i}"
            for idx, field in enumerate(fields):
                if field.startswith("st:Z:"):
                    fields[idx] = f"st:Z:{ts}"
                    break
            handle.write("\t".join(fields) + "\n")
PY
}

make_mixed_sam() {
  local output="$1"
  python3 - "$SOURCE_SAM" "$output" "$RECORD_GROUPS" "$READS_PER_GROUP" <<'PY'
from pathlib import Path
import sys

source = Path(sys.argv[1])
output = Path(sys.argv[2])
groups = int(sys.argv[3])
reads_per_group = int(sys.argv[4])

lines = source.read_text().splitlines()
header = [line for line in lines if line.startswith("@")]
body = [line for line in lines if not line.startswith("@")]

with output.open("w") as handle:
    for line in header:
        handle.write(line + "\n")

    for group in range(groups):
        for line in body:
            for i in range(reads_per_group // len(body)):
                fields = line.split("\t")
                fields[0] = f"{fields[0]}_{group}_{i}"
                handle.write("\t".join(fields) + "\n")
PY
}

benchmark_relative_window() {
  local label="$1"
  local bam="$2"
  local from_arg="$3"
  local to_arg="$4"

  rm -f "${bam}.ontime-index"
  printf '\n[%s]\n' "$label"
  time_run "first" "$BIN" -f "$from_arg" -t "$to_arg" "$bam"
  time_run "repeat1" "$BIN" -f "$from_arg" -t "$to_arg" "$bam"
  time_run "repeat2" "$BIN" -f "$from_arg" -t "$to_arg" "$bam"
}

main() {
  require_cmd cargo
  require_cmd python3
  require_cmd samtools

  mkdir -p "$WORK_DIR"

  cargo build --release --manifest-path "${ROOT_DIR}/Cargo.toml" >/dev/null

  local clustered_sam="${WORK_DIR}/clustered.sam"
  local clustered_bam="${WORK_DIR}/clustered.bam"
  local mixed_sam="${WORK_DIR}/mixed.sam"
  local mixed_bam="${WORK_DIR}/mixed.bam"

  make_clustered_sam "$clustered_sam"
  samtools view -bS -o "$clustered_bam" "$clustered_sam"

  make_mixed_sam "$mixed_sam"
  samtools view -bS -o "$mixed_bam" "$mixed_sam"

  printf 'groups=%s reads_per_group=%s\n' "$RECORD_GROUPS" "$READS_PER_GROUP"
  benchmark_relative_window "clustered" "$clustered_bam" "5h" "5h"
  benchmark_relative_window "mixed" "$mixed_bam" "0h" "-4h"
}

main "$@"
