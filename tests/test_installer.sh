#!/usr/bin/env bash
set -euo pipefail

root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
scratch=$(mktemp -d "${TMPDIR:-/tmp}/scpcqtl-installer-test.XXXXXX")
installer_bash=${SCPCQTL_TEST_BASH:-bash}
trap 'rm -rf "$scratch"' EXIT
mkdir -p "$scratch/source" "$scratch/java17/bin" "$scratch/java8/bin"

cat > "$scratch/java17/bin/java" <<'EOF'
#!/usr/bin/env bash
home=$(cd "$(dirname "$0")/.." && pwd)
if [[ "$*" == *-XshowSettings:properties* ]]; then
  printf '    java.home = %s\n' "$home" >&2
fi
printf 'openjdk version "17.0.15"\n' >&2
EOF
cat > "$scratch/java8/bin/java" <<'EOF'
#!/usr/bin/env bash
printf 'java version "1.8.0_452"\n' >&2
EOF
chmod 0755 "$scratch/java17/bin/java" "$scratch/java8/bin/java"

cat > "$scratch/source/nextflow" <<'EOF'
#!/usr/bin/env bash
if [[ "${1:-}" == -version ]]; then
  [[ "${JAVA_CMD:-}" == "${EXPECTED_JAVA_CMD:-}" ]] || {
    printf 'Nextflow received the wrong JAVA_CMD: %s\n' "${JAVA_CMD:-unset}" >&2
    exit 42
  }
  printf '%s\n' 'nextflow version 25.10.7'
fi
EOF
chmod 0755 "$scratch/source/nextflow"

run_installer() {
  PATH="$scratch/java17/bin:/usr/bin:/bin" \
  JAVA_HOME="$scratch/java8" \
  JAVA_CMD="$scratch/java8/bin/java" \
  EXPECTED_JAVA_CMD="$scratch/java17/bin/java" \
  SCPCQTL_INSTALL_ROOT="$1" \
  SCPCQTL_BIN_DIR="$2" \
  SCPCQTL_NEXTFLOW_SOURCE="$scratch/source/nextflow" \
  SCPCQTL_LAUNCHER_SOURCE="$root/bin/sc-pcqtl" \
  "$installer_bash" "$root/install.sh" --no-java-download >/dev/null
}

run_installer "$scratch/install" "$scratch/bin"
run_installer "$scratch/install" "$scratch/bin"

[[ -x "$scratch/bin/sc-pcqtl" ]]
[[ -x "$scratch/install/nextflow/25.10.7/nextflow" ]]
[[ -x "$scratch/install/libexec/sc-pcqtl.real" ]]
EXPECTED_JAVA_CMD="$scratch/java17/bin/java" \
  "$scratch/bin/sc-pcqtl" version | grep -Fq 'sc-pcqtl 0.1.0-dev'
EXPECTED_JAVA_CMD="$scratch/java17/bin/java" \
  "$scratch/bin/sc-pcqtl" nextflow -version | grep -Fq 'nextflow version 25.10.7'

run_installer "$scratch/collision" "$scratch/collision/bin"
[[ -x "$scratch/collision/bin/sc-pcqtl" ]]
[[ -x "$scratch/collision/libexec/sc-pcqtl.real" ]]
EXPECTED_JAVA_CMD="$scratch/java17/bin/java" \
  "$scratch/collision/bin/sc-pcqtl" version | grep -Fq 'sc-pcqtl 0.1.0-dev'

mkdir -p "$scratch/pointer/jvm/temurin-17-test/bin"
cp "$scratch/java17/bin/java" "$scratch/pointer/jvm/temurin-17-test/bin/java"
printf '%s\n' "$scratch/pointer/jvm/temurin-17-test" > "$scratch/pointer/java-home"
PATH="$scratch/java8/bin:/usr/bin:/bin" \
JAVA_HOME="$scratch/java8" JAVA_CMD="$scratch/java8/bin/java" \
EXPECTED_JAVA_CMD="$scratch/pointer/jvm/temurin-17-test/bin/java" \
SCPCQTL_INSTALL_ROOT="$scratch/pointer" \
SCPCQTL_BIN_DIR="$scratch/pointer-bin" \
SCPCQTL_NEXTFLOW_SOURCE="$scratch/source/nextflow" \
SCPCQTL_LAUNCHER_SOURCE="$root/bin/sc-pcqtl" \
"$installer_bash" "$root/install.sh" --no-java-download >/dev/null
EXPECTED_JAVA_CMD="$scratch/pointer/jvm/temurin-17-test/bin/java" \
  "$scratch/pointer-bin/sc-pcqtl" nextflow -version | \
  grep -Fq 'nextflow version 25.10.7'

if PATH="$scratch/java8/bin:/usr/bin:/bin" \
  JAVA_HOME="$scratch/java8" JAVA_CMD="$scratch/java8/bin/java" \
  SCPCQTL_INSTALL_ROOT="$scratch/incompatible" \
  SCPCQTL_BIN_DIR="$scratch/incompatible-bin" \
  SCPCQTL_NEXTFLOW_SOURCE="$scratch/source/nextflow" \
  SCPCQTL_LAUNCHER_SOURCE="$root/bin/sc-pcqtl" \
  "$installer_bash" "$root/install.sh" --no-java-download >/dev/null 2>&1; then
  printf 'Installer accepted Java 8 with --no-java-download.\n' >&2
  exit 1
fi
printf 'Installer tests passed.\n'
