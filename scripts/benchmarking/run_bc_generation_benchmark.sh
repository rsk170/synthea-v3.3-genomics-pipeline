#!/usr/bin/env bash

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
cd "$REPO_ROOT" || exit 1

BENCH_ROOT="${BENCH_ROOT:-output_runs/benchmark_bc_generation}"
AGE_RANGE="${AGE_RANGE:-45-90}"
KEEP_FILE="${KEEP_FILE:-src/main/resources/keep_modules/keep_breast_cancer.json}"
MAX_ATTEMPTS="${MAX_ATTEMPTS:-1000}"
mkdir -p "$BENCH_ROOT/logs"

lscpu > "$BENCH_ROOT/machine_lscpu.txt"
free -h > "$BENCH_ROOT/machine_memory.txt"
java -version > "$BENCH_ROOT/java_version.txt" 2>&1
python3 --version > "$BENCH_ROOT/python_version.txt" 2>&1

for SIZE in 10 50 100 500 1000; do
  for SEED in 12345 23456 34567; do
    RUN="bc_p${SIZE}_seed${SEED}"
    OUTDIR="$BENCH_ROOT/$RUN"
    LOGDIR="$BENCH_ROOT/logs/$RUN"

    mkdir -p "$OUTDIR" "$LOGDIR"

    /usr/bin/time -v -o "$LOGDIR/synthea.time.log" \
    ./run_synthea \
      -s "$SEED" \
      -p "$SIZE" \
      -a "$AGE_RANGE" \
      -k "$KEEP_FILE" \
      --exporter.csv.export=true \
      --exporter.baseDirectory="./$OUTDIR/" \
      --generate.max_attempts_to_keep_patient="$MAX_ATTEMPTS" \
      > "$LOGDIR/synthea.stdout.log" \
      2> "$LOGDIR/synthea.stderr.log"

    du -sh "$OUTDIR" > "$LOGDIR/storage_total.txt"
    du -sh "$OUTDIR/csv" > "$LOGDIR/storage_csv.txt"

    grep -E "Records:|Elapsed|Maximum resident|User time|System time|Exit status" \
      "$LOGDIR/synthea.stdout.log" "$LOGDIR/synthea.time.log" \
      > "$LOGDIR/summary_metrics.txt"
  done
done
