#!/usr/bin/env bash
set -euo pipefail

nextflow_version=${NEXTFLOW_VERSION:-25.10.7}
install_root=${SCPCQTL_INSTALL_ROOT:-"$HOME/.local/share/sc-pcqtl"}
bin_dir=${SCPCQTL_BIN_DIR:-"$HOME/.local/bin"}
launcher_url=${SCPCQTL_LAUNCHER_URL:-https://raw.githubusercontent.com/joookerz/sc-pcQTL/main/bin/sc-pcqtl}

usage() {
  cat <<'EOF'
Install the sc-pcQTL launcher and a pinned Nextflow into the user account.

Usage: bash install.sh [options]

Options:
  --prefix DIR             Installation root (default: ~/.local/share/sc-pcqtl)
  --bin-dir DIR            Launcher directory (default: ~/.local/bin)
  --nextflow-version VER   Nextflow version (default: 25.10.7)
  --no-java-download       Require an existing compatible Java installation
  -h, --help               Show this help

Environment variables with the same names as the defaults are also supported.
The script does not install a container runtime because Docker Desktop and
Linux/HPC runtimes require platform- or site-specific administration.
EOF
}

allow_java_download=true
while [[ $# -gt 0 ]]; do
  case "$1" in
    --prefix)
      [[ $# -ge 2 ]] || { printf '%s requires a value\n' "$1" >&2; exit 2; }
      install_root=$2
      shift 2
      ;;
    --bin-dir)
      [[ $# -ge 2 ]] || { printf '%s requires a value\n' "$1" >&2; exit 2; }
      bin_dir=$2
      shift 2
      ;;
    --nextflow-version)
      [[ $# -ge 2 ]] || { printf '%s requires a value\n' "$1" >&2; exit 2; }
      nextflow_version=$2
      shift 2
      ;;
    --no-java-download) allow_java_download=false; shift ;;
    -h|--help) usage; exit 0 ;;
    *) printf 'Unknown option: %s\n' "$1" >&2; usage >&2; exit 2 ;;
  esac
done

for tool in bash curl tar sed install mktemp uname; do
  command -v "$tool" >/dev/null 2>&1 || {
    printf 'Required command is missing: %s\n' "$tool" >&2
    exit 1
  }
done

java_major() {
  local executable=$1 line version major
  line=$("$executable" -version 2>&1 | head -n 1) || return 1
  version=$(printf '%s\n' "$line" | sed -E 's/.*version "([^"]+)".*/\1/')
  major=${version%%.*}
  if [[ "$major" == '1' ]]; then
    version=${version#1.}
    major=${version%%.*}
  fi
  [[ "$major" =~ ^[0-9]+$ ]] || return 1
  printf '%s\n' "$major"
}

compatible_java() {
  local executable=$1 major
  [[ -x "$executable" ]] || return 1
  major=$(java_major "$executable") || return 1
  (( major >= 17 ))
}

bundled_java_home() {
  if [[ -x "$install_root/java17/bin/java" ]]; then
    printf '%s\n' "$install_root/java17"
  elif [[ -x "$install_root/java17/Contents/Home/bin/java" ]]; then
    printf '%s\n' "$install_root/java17/Contents/Home"
  else
    return 1
  fi
}

mkdir -p "$install_root" "$bin_dir"
temporary_dir=$(mktemp -d "${TMPDIR:-/tmp}/scpcqtl-install.XXXXXX")
trap 'rm -rf "$temporary_dir"' EXIT

java_home=''
java_executable=''
if [[ -n "${JAVA_HOME:-}" ]] && compatible_java "$JAVA_HOME/bin/java"; then
  java_home=$JAVA_HOME
  java_executable=$JAVA_HOME/bin/java
elif command -v java >/dev/null 2>&1 && compatible_java "$(command -v java)"; then
  java_executable=$(command -v java)
  java_home=$("$java_executable" -XshowSettings:properties -version 2>&1 | \
    sed -n 's/^[[:space:]]*java.home = //p' | head -n 1)
  [[ -d "$java_home" ]] || java_home=''
elif java_home=$(bundled_java_home 2>/dev/null); then
  java_executable=$java_home/bin/java
fi

if [[ -z "$java_executable" ]]; then
  [[ "$allow_java_download" == true ]] || {
    printf 'Java 17 or newer was not found.\n' >&2
    exit 1
  }
  case "$(uname -s)" in
    Linux) adoptium_os=linux ;;
    Darwin) adoptium_os=mac ;;
    *) printf 'Automatic Java installation supports Linux and macOS only.\n' >&2; exit 1 ;;
  esac
  case "$(uname -m)" in
    x86_64|amd64) adoptium_arch=x64 ;;
    arm64|aarch64) adoptium_arch=aarch64 ;;
    *) printf 'Unsupported CPU architecture: %s\n' "$(uname -m)" >&2; exit 1 ;;
  esac
  java_url="https://api.adoptium.net/v3/binary/latest/17/ga/${adoptium_os}/${adoptium_arch}/jre/hotspot/normal/eclipse?project=jdk"
  printf 'Installing Temurin Java 17 under %s\n' "$install_root/java17"
  curl -fL --retry 3 "$java_url" -o "$temporary_dir/java.tar.gz"
  mkdir -p "$temporary_dir/java17"
  tar -xzf "$temporary_dir/java.tar.gz" -C "$temporary_dir/java17" --strip-components=1
  if [[ ! -e "$install_root/java17" ]]; then
    mv "$temporary_dir/java17" "$install_root/java17"
  fi
  java_home=$(bundled_java_home)
  java_executable=$java_home/bin/java
fi

nextflow_dir="$install_root/nextflow/$nextflow_version"
nextflow_executable="$nextflow_dir/nextflow"
if [[ ! -x "$nextflow_executable" ]]; then
  mkdir -p "$nextflow_dir" "$temporary_dir/nextflow"
  if [[ -n "${SCPCQTL_NEXTFLOW_SOURCE:-}" ]]; then
    install -m 0755 "$SCPCQTL_NEXTFLOW_SOURCE" "$nextflow_executable"
  else
    printf 'Installing Nextflow %s under %s\n' "$nextflow_version" "$nextflow_dir"
    curl -fsSL https://get.nextflow.io -o "$temporary_dir/get-nextflow.sh"
    (
      cd "$temporary_dir/nextflow"
      PATH="$(dirname "$java_executable"):$PATH" NXF_VER="$nextflow_version" \
        bash "$temporary_dir/get-nextflow.sh" >/dev/null
    )
    install -m 0755 "$temporary_dir/nextflow/nextflow" "$nextflow_executable"
  fi
fi

launcher_source=''
script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" 2>/dev/null && pwd || true)
if [[ -n "${SCPCQTL_LAUNCHER_SOURCE:-}" ]]; then
  launcher_source=$SCPCQTL_LAUNCHER_SOURCE
elif [[ -n "$script_dir" && -f "$script_dir/bin/sc-pcqtl" ]]; then
  launcher_source=$script_dir/bin/sc-pcqtl
else
  launcher_source=$temporary_dir/sc-pcqtl
  curl -fsSL "$launcher_url" -o "$launcher_source"
fi
mkdir -p "$install_root/bin"
install -m 0755 "$launcher_source" "$install_root/bin/sc-pcqtl"

quoted_nextflow=$(printf '%q' "$nextflow_executable")
quoted_launcher=$(printf '%q' "$install_root/bin/sc-pcqtl")
{
  printf '%s\n' '#!/usr/bin/env bash'
  if [[ -n "$java_home" ]]; then
    quoted_java_home=$(printf '%q' "$java_home")
    printf 'export JAVA_HOME=%s\n' "$quoted_java_home"
    printf '%s\n' 'export PATH="$JAVA_HOME/bin:$PATH"'
  fi
  cat <<EOF
export SCPCQTL_NEXTFLOW=$quoted_nextflow
exec $quoted_launcher "\$@"
EOF
} > "$bin_dir/sc-pcqtl"
chmod 0755 "$bin_dir/sc-pcqtl"

printf '\nInstalled sc-pcqtl: %s\n' "$bin_dir/sc-pcqtl"
case ":$PATH:" in
  *":$bin_dir:"*) ;;
  *) printf 'Add this directory to PATH: export PATH="%s:$PATH"\n' "$bin_dir" ;;
esac
printf 'A Docker/Podman runtime is required on macOS and workstations; Apptainer is also supported on Linux/HPC.\n'
printf 'After starting the runtime, run: sc-pcqtl doctor\n'
