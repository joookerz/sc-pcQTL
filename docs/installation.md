# Installation

## Recommended container installation

Install Java 17 and Nextflow 25.10.0 or newer, then provide one supported container
runtime. Analysis R packages do not need to be installed on the host.

```bash
java -version
export NXF_VER=25.10.7
curl -fsSL https://get.nextflow.io | bash
install -m 0755 nextflow "$HOME/.local/bin/nextflow"
nextflow -version
```

For local execution install Docker. On HPC, Apptainer is recommended; legacy
Singularity remains supported. The corresponding executable must be visible
on compute nodes, not only the login node.

The default `edge` images are built from the repository's `main` branch:

```text
ghcr.io/joookerz/sc-pcqtl-core:edge
ghcr.io/joookerz/sc-pcqtl-saigeqtl:edge
```

For a released analysis, run a tagged workflow revision and its matching
container tags. Container names can be overridden with `--core_container` and
`--saige_container`.

## Optional launcher environment

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
