#!/usr/bin/env bash
set -euo pipefail

output=''
for argument in "$@"; do
  case "${argument}" in
    --SAIGEOutputFile=*) output=${argument#*=} ;;
  esac
done
[[ -n "${output}" ]] || { echo 'Missing --SAIGEOutputFile' >&2; exit 2; }
printf 'CHR\tPOS\tMarkerID\tAllele1\tAllele2\tAC_Allele2\tAF_Allele2\tMissingRate\tBETA\tSE\tTstat\tvar\tp.value\tN\n' > "${output}"
printf '1\t101000\t1:101000_A_G\tA\tG\t6\t0.30\t0\t0.40\t0.10\t4\t1\t0.001\t10\n' >> "${output}"
printf '1\t102000\t1:102000_C_T\tC\tT\t8\t0.40\t0\t0.02\t0.10\t0.2\t1\t0.5\t10\n' >> "${output}"
