#!/usr/bin/env bash
set -euo pipefail

root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
scratch=$(mktemp -d "${TMPDIR:-/tmp}/scpcqtl-resume-reports.XXXXXX")
trap 'rm -rf "$scratch"' EXIT

cat > "$scratch/main.nf" <<'EOF'
nextflow.enable.dsl=2

process PROBE {
    output:
    path 'probe.txt'
    script:
    "printf 'ok\\n' > probe.txt"
}

workflow { PROBE() }
EOF

run=(nextflow -c "$root/nextflow.config" run "$scratch/main.nf"
  -profile test -ansi-log false -work-dir "$scratch/work"
  --outdir "$scratch/output")
"${run[@]}" > "$scratch/first.log" 2>&1
"${run[@]}" -resume > "$scratch/resume.log" 2>&1

for report in execution_report.html execution_timeline.html \
  execution_trace.txt pipeline_dag.html; do
  [[ -s "$scratch/output/pipeline_info/$report" ]]
done
grep -Fq 'CACHED' "$scratch/output/pipeline_info/execution_trace.txt"
if grep -Eiq 'already exists|null.?writer|nullpointer' "$scratch/resume.log"; then
  printf 'Resume emitted a telemetry overwrite/observer error.\n' >&2
  cat "$scratch/resume.log" >&2
  exit 1
fi
printf 'Resume telemetry test passed.\n'
