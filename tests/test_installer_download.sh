#!/usr/bin/env bash
set -euo pipefail

root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
scratch=$(mktemp -d "${TMPDIR:-/tmp}/scpcqtl-installer-download.XXXXXX")
installer_bash=${SCPCQTL_TEST_BASH:-bash}
trap 'rm -rf "$scratch"' EXIT

SCPCQTL_INSTALL_ROOT="$scratch/install" \
SCPCQTL_BIN_DIR="$scratch/bin" \
SCPCQTL_LAUNCHER_SOURCE="$root/bin/sc-pcqtl" \
"$installer_bash" "$root/install.sh" --no-java-download >/dev/null

version_output=$("$scratch/bin/sc-pcqtl" nextflow -version 2>&1)
grep -Fq 'version 25.10.7' <<< "$version_output"
[[ -s "$scratch/install/nextflow-home/framework/25.10.7/nextflow-25.10.7-one.jar" ]]
printf 'Real Nextflow bootstrap test passed.\n'
