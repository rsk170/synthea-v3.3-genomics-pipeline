#!/usr/bin/env bash

set -u

if [[ $# -lt 1 || $# -gt 2 ]]; then
  echo "Usage: $0 SOURCE_RUN_DIR [RUN_LABEL]" >&2
  echo "Example: $0 output_runs/benchmark_bc_generation_final/bc_p50_seed12345" >&2
  exit 2
fi

SOURCE_RUN="$1"
RUN_LABEL="${2:-$(basename "$SOURCE_RUN")}"

cd /home/rkisleva/SDG_tools/synthea_developer3.3.0/synthea || exit 1

BENCH_ROOT="${BENCH_ROOT:-output_runs/benchmark_bc_postprocessing}"
WORK_RUN="$BENCH_ROOT/$RUN_LABEL"
LOGDIR="$BENCH_ROOT/logs/$RUN_LABEL"
CSV_DIR="$WORK_RUN/csv"
PASSENGER_POOL="$WORK_RUN/passenger_pool/breast_cancer_passenger_only_from_maf.maf"

if [[ ! -d "$SOURCE_RUN/csv" ]]; then
  echo "Source run is missing a csv/ directory: $SOURCE_RUN" >&2
  exit 3
fi

if [[ -e "$WORK_RUN" || -e "$LOGDIR" ]]; then
  echo "Refusing to overwrite existing benchmark output:" >&2
  echo "  $WORK_RUN" >&2
  echo "  $LOGDIR" >&2
  echo "Choose a different RUN_LABEL or BENCH_ROOT, or move the existing output first." >&2
  exit 4
fi

mkdir -p "$WORK_RUN" "$LOGDIR" "$(dirname "$PASSENGER_POOL")"

lscpu > "$BENCH_ROOT/machine_lscpu.txt"
free -h > "$BENCH_ROOT/machine_memory.txt"
python3 --version > "$BENCH_ROOT/python_version.txt" 2>&1

{
  echo "run_label=$RUN_LABEL"
  echo "source_run=$SOURCE_RUN"
  echo "work_run=$WORK_RUN"
  echo "start_time=$(date -Is)"
  echo "copy_start_time=$(date -Is)"
} > "$LOGDIR/run_status.txt"

cp -a "$SOURCE_RUN/csv" "$WORK_RUN/"

{
  echo "copy_end_time=$(date -Is)"
  echo "copied_csv_bytes=$(du -sb "$CSV_DIR" | awk '{print $1}')"
} >> "$LOGDIR/run_status.txt"

write_storage_snapshot() {
  local step="$1"
  du -sb "$WORK_RUN" > "$LOGDIR/${step}.storage_total_bytes.txt"
  du -sh "$WORK_RUN" > "$LOGDIR/${step}.storage_total_human.txt"
  if [[ -d "$CSV_DIR" ]]; then
    du -sb "$CSV_DIR" > "$LOGDIR/${step}.storage_csv_bytes.txt"
    du -sh "$CSV_DIR" > "$LOGDIR/${step}.storage_csv_human.txt"
  fi
  if [[ -d "$WORK_RUN/maf_files" ]]; then
    du -sb "$WORK_RUN/maf_files" > "$LOGDIR/${step}.storage_maf_files_bytes.txt"
    du -sh "$WORK_RUN/maf_files" > "$LOGDIR/${step}.storage_maf_files_human.txt"
  fi
}

run_step() {
  local step="$1"
  shift
  local -a cmd=("$@")

  {
    echo "step_start_time=$step $(date -Is)"
    printf "step_command_%s=" "$step"
    printf "%q " "${cmd[@]}"
    echo
  } >> "$LOGDIR/run_status.txt"

  /usr/bin/time -v -o "$LOGDIR/${step}.time.log" \
    "${cmd[@]}" \
    > "$LOGDIR/${step}.stdout.log" \
    2> "$LOGDIR/${step}.stderr.log"
  local status=$?

  {
    echo "step_exit_status=$step $status"
    echo "step_end_time=$step $(date -Is)"
  } >> "$LOGDIR/run_status.txt"

  write_storage_snapshot "$step"

  if [[ "$status" -ne 0 ]]; then
    echo "Step failed: $step (exit status $status)" >&2
    echo "See logs: $LOGDIR" >&2
    exit "$status"
  fi
}

write_storage_snapshot "initial"

run_step "01_clone_groups" \
  python3 scripts/build_breast_cancer_clones.py \
    --observations "$CSV_DIR/observations.csv" \
    --output "$CSV_DIR/breast_cancer_clone_groups.csv"

run_step "02_clone_proportions" \
  python3 scripts/build_breast_cancer_clone_proportions.py \
    --clone-groups "$CSV_DIR/breast_cancer_clone_groups.csv" \
    --output "$CSV_DIR/breast_cancer_clone_proportions.csv"

run_step "03_pruned_observations" \
  python3 scripts/build_breast_cancer_pruned_observations.py \
    --clone-proportions "$CSV_DIR/breast_cancer_clone_proportions.csv" \
    --observations "$CSV_DIR/observations.csv" \
    --medications "$CSV_DIR/medications.csv" \
    --civic-driver-variants scripts/civic_breast_cancer_driver_variants_from_maf.csv \
    --driver-variants scripts/breast_cancer_driver_variants_from_maf.csv \
    --non-disruptive-variants scripts/breast_cancer_non_disruptive_variants_from_maf.csv \
    --output "$CSV_DIR/observations_pruned_by_clone_vaf.csv"

run_step "04_passenger_mutations" \
  python3 scripts/build_breast_cancer_passenger_mutations.py \
    --maf scripts/final_consensus_passonly.snv_mnv_indel.icgc.controlled.maf \
    --driver-tsv scripts/TableS3_panorama_driver_mutations_ICGC_samples.controlled.tsv \
    --patients "$CSV_DIR/patients.csv" \
    --observations "$CSV_DIR/observations.csv" \
    --passenger-maf "$PASSENGER_POOL" \
    --output "$CSV_DIR/breast_cancer_assigned_passenger_mutations.tsv"

run_step "05_complete_maf_files" \
  python3 scripts/build_breast_cancer_complete_maf_files.py \
    --pruned-observations "$CSV_DIR/observations_pruned_by_clone_vaf.csv" \
    --assigned-passengers "$CSV_DIR/breast_cancer_assigned_passenger_mutations.tsv" \
    --clone-groups "$CSV_DIR/breast_cancer_clone_groups.csv" \
    --clone-proportions "$CSV_DIR/breast_cancer_clone_proportions.csv" \
    --driver-variants scripts/breast_cancer_driver_variants_from_maf.csv \
    --civic-driver-variants scripts/civic_breast_cancer_driver_variants_from_maf.csv \
    --non-disruptive-variants scripts/breast_cancer_non_disruptive_variants_from_maf.csv \
    --maf scripts/final_consensus_passonly.snv_mnv_indel.icgc.controlled.maf \
    --output-dir "$WORK_RUN/maf_files"

count_rows() {
  local file="$1"
  if [[ -f "$file" ]]; then
    tail -n +2 "$file" | wc -l
  else
    echo ""
  fi
}

{
  echo "end_time=$(date -Is)"
  echo "patients_csv_rows=$(count_rows "$CSV_DIR/patients.csv")"
  echo "clone_groups_rows=$(count_rows "$CSV_DIR/breast_cancer_clone_groups.csv")"
  echo "clone_proportions_rows=$(count_rows "$CSV_DIR/breast_cancer_clone_proportions.csv")"
  echo "pruned_observations_rows=$(count_rows "$CSV_DIR/observations_pruned_by_clone_vaf.csv")"
  echo "assigned_passenger_rows=$(count_rows "$CSV_DIR/breast_cancer_assigned_passenger_mutations.tsv")"
  if [[ -d "$WORK_RUN/maf_files" ]]; then
    echo "maf_files=$(find "$WORK_RUN/maf_files" -type f -name '*.maf' | wc -l)"
  else
    echo "maf_files=0"
  fi
  if [[ -f "$PASSENGER_POOL" ]]; then
    echo "passenger_pool_bytes=$(stat -c '%s' "$PASSENGER_POOL")"
  fi
} >> "$LOGDIR/run_status.txt"

echo "Post-processing benchmark finished successfully"
echo "Logs: $LOGDIR"
echo "Working output: $WORK_RUN"
