#!/usr/bin/env bash
set -euo pipefail

root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
runtime=$(mktemp -d "${TMPDIR:-/tmp}/scpcqtl-qtl-test.XXXXXX")
trap 'rm -rf "${runtime}"' EXIT
profile=${SCPCQTL_TEST_PROFILE:-test}

run_rscript() {
  if [[ -n "${SCPCQTL_TEST_CONTAINER:-}" ]]; then
    docker run --rm -v "${root}:${root}" -v "${runtime}:${runtime}" -w "${root}" \
      "${SCPCQTL_TEST_CONTAINER}" Rscript "$@"
  else
    Rscript "$@"
  fi
}

for chromosome in $(seq 1 22); do
  for extension in bed bim fam; do
    printf 'mock\n' > "${runtime}/genotype_chr${chromosome}.${extension}"
  done
done
for extension in bed bim fam; do
  printf 'mock\n' > "${runtime}/vr.${extension}"
done

export PATH="${root}/tests/mocks:${PATH}"
nextflow run "${root}" -profile "${profile}" -ansi-log false \
  -work-dir "${runtime}/work" --outdir "${runtime}/output" \
  --run_qtl true --genotype_prefix "${runtime}/genotype_chr{chr}" \
  --variance_ratio_prefix "${runtime}/vr" "$@"
run_rscript "${root}/tests/assert_qtl_outputs.R" "${runtime}/output"
