#!/usr/bin/env bash
set -euo pipefail

output=''
gene='PC1'
for argument in "$@"; do
  case "${argument}" in
    --genePval_outputFile=*) output=${argument#*=} ;;
    --geneName=*) gene=${argument#*=} ;;
  esac
done
[[ -n "${output}" ]] || { echo 'Missing --genePval_outputFile' >&2; exit 2; }
printf 'gene\tACAT_p\ttop_MarkerID\ttop_pval\n%s\t0.01\t1:101000_A_G\t0.001\n' "${gene}" > "${output}"
