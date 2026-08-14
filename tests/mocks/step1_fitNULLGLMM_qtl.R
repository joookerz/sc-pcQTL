#!/usr/bin/env bash
set -euo pipefail

prefix=''
for argument in "$@"; do
  case "${argument}" in
    --outputPrefix=*) prefix=${argument#*=} ;;
  esac
done
[[ -n "${prefix}" ]] || { echo 'Missing --outputPrefix' >&2; exit 2; }
printf 'mock model\n' > "${prefix}.rda"
printf 'mock variance ratio\n' > "${prefix}.varianceRatio.txt"
