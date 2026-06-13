#!/usr/bin/env bash

set -u

if [[ $# -lt 2 || $# -gt 3 ]]; then
  echo "Usage: $0 COHORT_SIZE SEED [TIME_LIMIT]" >&2
  echo "Example: $0 50 12345 4h" >&2
  exit 2
fi

COHORT_SIZE="$1"
SEED="$2"
TIME_LIMIT="${3:-4h}"

cd /home/rkisleva/SDG_tools/synthea_developer3.3.0/synthea || exit 1

BENCH_ROOT="${BENCH_ROOT:-output_runs/benchmark_bc_generation_single}"
AGE_RANGE="${AGE_RANGE:-45-90}"
CONFIG_FILE="${CONFIG_FILE:-src/main/resources/synthea.properties}"
KEEP_FILE="${KEEP_FILE:-src/main/resources/keep_modules/keep_breast_cancer.json}"
MODULES_DIR="${MODULES_DIR:-src/main/resources/modules}"
MAX_ATTEMPTS="${MAX_ATTEMPTS:-1000}"
REFERENCE_DATE="${REFERENCE_DATE:-20260507}"
END_DATE="${END_DATE:-20260507}"
JAVA_MAX_HEAP="${JAVA_MAX_HEAP:-2g}"

RUN="bc_p${COHORT_SIZE}_seed${SEED}"
OUTDIR="$BENCH_ROOT/$RUN"
LOGDIR="$BENCH_ROOT/logs/$RUN"

if [[ -e "$OUTDIR" || -e "$LOGDIR" ]]; then
  echo "Refusing to overwrite existing benchmark output:" >&2
  echo "  $OUTDIR" >&2
  echo "  $LOGDIR" >&2
  echo "Choose a different SEED, BENCH_ROOT, or move the existing output first." >&2
  exit 3
fi

mkdir -p "$OUTDIR" "$LOGDIR"

lscpu > "$BENCH_ROOT/machine_lscpu.txt"
free -h > "$BENCH_ROOT/machine_memory.txt"
java -version > "$BENCH_ROOT/java_version.txt" 2>&1
python3 --version > "$BENCH_ROOT/python_version.txt" 2>&1

SYNTHEA_JAR="${SYNTHEA_JAR:-build/libs/synthea-with-dependencies.jar}"

if [[ ! -f "$SYNTHEA_JAR" ]]; then
  echo "Synthea jar not found: $SYNTHEA_JAR" >&2
  echo "Build it first with: ./gradlew build" >&2
  exit 4
fi

SYNTHEA_ARGS=(
  -c "$CONFIG_FILE"
  -s "$SEED"
  -cs "$SEED"
  -r "$REFERENCE_DATE"
  -e "$END_DATE"
  -p "$COHORT_SIZE"
  -a "$AGE_RANGE"
  -d "$MODULES_DIR"
  -k "$KEEP_FILE"
  --exporter.csv.export=true
  --exporter.baseDirectory="./$OUTDIR/"
  --generate.max_attempts_to_keep_patient="$MAX_ATTEMPTS"
)

CMD=(
  java
  "-Xmx$JAVA_MAX_HEAP"
  -jar
  "$SYNTHEA_JAR"
  "${SYNTHEA_ARGS[@]}"
)

{
  echo "run=$RUN"
  echo "cohort_size=$COHORT_SIZE"
  echo "seed=$SEED"
  echo "time_limit=$TIME_LIMIT"
  echo "age_range=$AGE_RANGE"
  echo "config_file=$CONFIG_FILE"
  echo "keep_file=$KEEP_FILE"
  echo "modules_dir=$MODULES_DIR"
  echo "synthea_jar=$SYNTHEA_JAR"
  echo "java_max_heap=$JAVA_MAX_HEAP"
  echo "max_attempts=$MAX_ATTEMPTS"
  echo "reference_date=$REFERENCE_DATE"
  echo "end_date=$END_DATE"
  echo "start_time=$(date -Is)"
  printf "synthea_args="
  printf "%q " "${SYNTHEA_ARGS[@]}"
  echo
  printf "command="
  printf "%q " "${CMD[@]}"
  echo
} > "$LOGDIR/run_status.txt"

/usr/bin/time -v -o "$LOGDIR/synthea.time.log" \
  timeout --kill-after=60s "$TIME_LIMIT" "${CMD[@]}" \
  > "$LOGDIR/synthea.stdout.log" \
  2> "$LOGDIR/synthea.stderr.log"
STATUS=$?

{
  echo "exit_status=$STATUS"
  if [[ "$STATUS" -eq 124 || "$STATUS" -eq 137 ]]; then
    echo "timed_out=true"
  else
    echo "timed_out=false"
  fi
  echo "end_time=$(date -Is)"
} >> "$LOGDIR/run_status.txt"

if [[ -d "$OUTDIR" ]]; then
  du -sh "$OUTDIR" > "$LOGDIR/storage_total.txt"
  du -sb "$OUTDIR" > "$LOGDIR/storage_total_bytes.txt"
fi

if [[ -d "$OUTDIR/csv" ]]; then
  du -sh "$OUTDIR/csv" > "$LOGDIR/storage_csv.txt"
  du -sb "$OUTDIR/csv" > "$LOGDIR/storage_csv_bytes.txt"
fi

if [[ -f "$OUTDIR/csv/patients.csv" ]]; then
  PATIENT_ROWS=$(tail -n +2 "$OUTDIR/csv/patients.csv" | wc -l)
  echo "patients_csv_rows=$PATIENT_ROWS" >> "$LOGDIR/run_status.txt"
fi

{
  grep -En "Population:|Seed:|Provider Seed:|Reference Time:|Records:|RNG=|Clinician RNG=" \
    "$LOGDIR/synthea.stdout.log" || true
  grep -En "Elapsed|Maximum resident|User time|System time|Exit status" \
    "$LOGDIR/synthea.time.log" || true
} > "$LOGDIR/summary_metrics.txt"

echo "Benchmark run finished with exit status $STATUS"
echo "Logs: $LOGDIR"
echo "Output: $OUTDIR"

exit "$STATUS"
