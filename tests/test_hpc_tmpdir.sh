#!/usr/bin/env bash
set -euo pipefail

root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
scratch=$(mktemp -d "${TMPDIR:-/tmp}/scpcqtl-hpc-tmpdir.XXXXXX")
trap 'rm -rf "$scratch"' EXIT

for runtime in apptainer singularity; do
  config=$(TMPDIR="$scratch" nextflow config -profile "${runtime},example" -flat "$root")
  grep -Fq -- "--env TMPDIR=${scratch}" <<< "$config"
  grep -Fq -- "--bind ${scratch}:${scratch}" <<< "$config"
  grep -Fq -- '--env R_LIBS_USER=/opt/scpcqtl-user-library-disabled' <<< "$config"
done

sge_config=$(TMPDIR="$scratch" nextflow config -profile apptainer,sge,example -flat "$root")
grep -Fq "process.executor = 'sge'" <<< "$sge_config"
grep -Fq 'executor.queueSize = 100' <<< "$sge_config"

printf 'not a directory\n' > "$scratch/not-a-directory"
if TMPDIR="$scratch/not-a-directory" \
  nextflow config -profile apptainer,example -flat "$root" >/dev/null 2>&1; then
  printf 'Apptainer profile accepted a non-directory TMPDIR.\n' >&2
  exit 1
fi
printf 'HPC TMPDIR configuration tests passed.\n'
