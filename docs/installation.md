# Installation

sc-pcQTL separates the lightweight workflow launcher from the analysis
software. The host needs Bash, a container runtime, Java 17 or newer, and
Nextflow 25.10.0 or newer. R, fasthurdle, PLINK, and SAIGE-QTL are supplied by
pinned containers.

The supported host combinations are Linux x86_64 with Docker, Podman,
Apptainer, or Singularity; macOS with Docker Desktop or Podman; and WSL2 with
Docker Desktop. ARM64 hosts require Docker/Podman amd64 emulation because the
current analysis images are `linux/amd64`.

## One-command user installation

Install and start a container runtime first, then run:

```bash
curl -fsSL https://raw.githubusercontent.com/joookerz/sc-pcQTL/main/install.sh | bash
export PATH="$HOME/.local/bin:$PATH"
sc-pcqtl doctor
sc-pcqtl doctor --deep
```

The installer works without administrator privileges on Linux and macOS. It:

1. reuses a compatible Java already on `PATH` or in `JAVA_HOME`;
2. otherwise installs a private Temurin Java 17 runtime under
   `~/.local/share/sc-pcqtl`;
3. installs Nextflow 25.10.7 under the same directory; and
4. installs `sc-pcqtl` into `~/.local/bin`.

It does not install Docker, Podman, Apptainer, or Singularity because those
runtimes require platform- or site-specific setup.

While the repository is private, collaborators should install from an
authenticated clone instead of the public `curl` command:

```bash
git clone git@github.com:joookerz/sc-pcQTL.git
cd sc-pcQTL
bash install.sh
```

They must also authenticate their selected container runtime to GHCR if the
analysis images remain private. GitHub package visibility is configured
separately from repository visibility; public package visibility permits
anonymous pulls.

After installation, persist the launcher path in the shell configuration:

```bash
printf '\nexport PATH="$HOME/.local/bin:$PATH"\n' >> "$HOME/.bashrc"
```

For the default macOS zsh, write the same line to `~/.zshrc` instead.

## Platform setup

### macOS

Install [Docker Desktop](https://docs.docker.com/desktop/setup/install/mac-install/)
or [Podman](https://podman.io/docs/installation) and start its Linux virtual machine.
Both Intel and Apple Silicon Macs are supported. The analysis images are
pinned to `linux/amd64`; Apple Silicon therefore uses transparent emulation
and may be slower for the SAIGE-QTL stage.

Keep inputs under a directory shared with the container runtime. Docker
Desktop shares the user home directory by default; external volumes may need
to be added in Docker Desktop settings.

```bash
sc-pcqtl doctor
sc-pcqtl example --outdir results/example_component_union
```

### Linux workstation

Install either [Docker Engine](https://docs.docker.com/engine/install/) or
[Podman](https://podman.io/docs/installation) using the operating-system package
manager. Ensure the current user can run the selected runtime, then install
the launcher and verify it:

```bash
sc-pcqtl doctor
sc-pcqtl doctor --deep
```

For rootless Podman, the operating system must provide subordinate UID/GID
ranges and cgroup delegation for the current user. A successful `podman info`
does not prove that image pulls, bind mounts, and Nextflow resource flags work;
`doctor --deep` tests those operations with a small public image. If the deep
check fails, configure rootless Podman with the system administrator or use
Docker/Apptainer instead.

The launcher prefers Apptainer, Singularity, Docker, then Podman on Linux. To
override detection:

```bash
export SCPCQTL_RUNTIME=docker
```

### Linux HPC

[Apptainer](https://apptainer.org/docs/admin/main/installation.html) is
recommended; Singularity remains supported. The runtime must be
available on compute nodes as well as the login node. Install the launcher in
the user account and add the Slurm profile with a site-specific configuration:

```bash
SCPCQTL_RUNTIME=apptainer SCPCQTL_EXTRA_PROFILES=slurm \
sc-pcqtl run \
  -c institutional.config \
  [parameters]
```

Set a shared cache to prevent repeated image downloads:

```bash
export NXF_APPTAINER_CACHEDIR=/shared/path/sc-pcqtl-containers
```

### Windows

Use Windows Subsystem for Linux 2 (WSL2), install Docker Desktop with WSL
integration, and run the Linux installer from the WSL shell. Native Windows
PowerShell execution is not supported by Nextflow.

## Advanced installer options

```bash
bash install.sh \
  --prefix "$HOME/apps/sc-pcqtl" \
  --bin-dir "$HOME/bin" \
  --nextflow-version 25.10.7
```

Use `--no-java-download` to require an existing Java 17+ installation. The
same settings can be supplied with `SCPCQTL_INSTALL_ROOT`,
`SCPCQTL_BIN_DIR`, and `NEXTFLOW_VERSION`.

## Direct Nextflow installation

Users who do not want the launcher can follow the standard Nextflow setup:

```bash
java -version
export NXF_VER=25.10.7
curl -fsSL https://get.nextflow.io | bash
install -m 0755 nextflow "$HOME/.local/bin/nextflow"
nextflow -version
```

See the [official Nextflow installation guide](https://docs.seqera.io/nextflow/install)
for supported Java and shell requirements.

The default `edge` images are built from the repository's `main` branch:

```text
ghcr.io/joookerz/sc-pcqtl-core:edge
ghcr.io/joookerz/sc-pcqtl-saigeqtl:edge
```

For a released analysis, run a tagged workflow revision and its matching
container tags. Container names can be overridden with `--core_container` and
`--saige_container`.

sc-pcQTL is CPU-only. It does not use or require a GPU; compute planning should
focus on CPUs, RAM, storage throughput, and scheduler concurrency.

## Conda launcher environment

`environment-launcher.yml` installs only Java and Nextflow through Conda. A
container runtime remains a host or cluster prerequisite.

```bash
conda env create -f environment-launcher.yml
conda activate sc-pcqtl-launcher
```

## Clone for development

```bash
git clone https://github.com/joookerz/sc-pcQTL.git
cd sc-pcQTL
nextflow run . --help
```
