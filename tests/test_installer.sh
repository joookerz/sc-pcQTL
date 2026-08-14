#!/usr/bin/env bash
set -euo pipefail

root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
scratch=$(mktemp -d "${TMPDIR:-/tmp}/scpcqtl-installer-test.XXXXXX")
trap 'rm -rf "$scratch"' EXIT
mkdir -p "$scratch/source"

cat > "$scratch/source/nextflow" <<'EOF'
#!/usr/bin/env bash
if [[ "${1:-}" == -version ]]; then
  printf '%s\n' 'nextflow version 25.10.7'
fi
EOF
chmod 0755 "$scratch/source/nextflow"

SCPCQTL_INSTALL_ROOT="$scratch/install" \
SCPCQTL_BIN_DIR="$scratch/bin" \
SCPCQTL_NEXTFLOW_SOURCE="$scratch/source/nextflow" \
SCPCQTL_LAUNCHER_SOURCE="$root/bin/sc-pcqtl" \
bash "$root/install.sh" --no-java-download >/dev/null

[[ -x "$scratch/bin/sc-pcqtl" ]]
[[ -x "$scratch/install/nextflow/25.10.7/nextflow" ]]
"$scratch/bin/sc-pcqtl" version | grep -Fq 'sc-pcqtl 0.1.0-dev'
printf 'Installer tests passed.\n'
