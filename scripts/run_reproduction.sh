#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "$0")/.." && pwd)"
export FRACTIONAL_LASER_PROJECT_ROOT="$repo_root"
export FRACTIONAL_LASER_DATA_DIR="${FRACTIONAL_LASER_DATA_DIR:-$repo_root/data}"

required=(
  GSE168760_series_matrix.txt.gz
  GSE206495_series_matrix.txt.gz
  GPL13667_platform.soft.tsv
  GPL15207_platform.soft.tsv
)

for filename in "${required[@]}"; do
  if [[ ! -f "$FRACTIONAL_LASER_DATA_DIR/$filename" ]]; then
    echo "Missing required input: $FRACTIONAL_LASER_DATA_DIR/$filename" >&2
    exit 2
  fi
done

mkdir -p "$repo_root/reproduction_output"
Rscript "$repo_root/scripts/phase1b_temporal_validation.R"
Rscript "$repo_root/scripts/phase2b_formal_dual_cohort.R"
python3 "$repo_root/scripts/phase3b_generate_all.py"

echo "Reproduction scripts completed. Review stopping checks before interpretation."

