#!/usr/bin/env bash
set -euo pipefail

root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
scratch=$(mktemp -d "${TMPDIR:-/tmp}/scpcqtl-clean-install.XXXXXX")
trap 'rm -rf "$scratch"' EXIT
mkdir -p "$scratch/java8/bin"

cat > "$scratch/java8/bin/java" <<'EOF'
#!/usr/bin/env bash
printf 'java version "1.8.0_452"\n' >&2
EOF
chmod 0755 "$scratch/java8/bin/java"

run_installer() {
  PATH="$scratch/java8/bin:/usr/bin:/bin" \
  JAVA_HOME="$scratch/java8" \
  JAVA_CMD="$scratch/java8/bin/java" \
  SCPCQTL_INSTALL_ROOT="$scratch/install" \
  SCPCQTL_BIN_DIR="$scratch/bin" \
  SCPCQTL_LAUNCHER_SOURCE="$root/bin/sc-pcqtl" \
  bash "$root/install.sh" >/dev/null
}

run_installer
java_home=$(head -n 1 "$scratch/install/java-home")
[[ -x "$java_home/bin/java" ]]
"$java_home/bin/java" -version 2>&1 | grep -Eq 'version "17([.]|\")'
"$scratch/bin/sc-pcqtl" nextflow -version 2>&1 | grep -Fq 'version 25.10.7'

run_installer
[[ "$(head -n 1 "$scratch/install/java-home")" == "$java_home" ]]
printf 'Clean-host Java and Nextflow installation tests passed.\n'
