# sc-pcQTL

[![CI](https://github.com/joookerz/sc-pcQTL/actions/workflows/ci.yml/badge.svg)](https://github.com/joookerz/sc-pcQTL/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Nextflow](https://img.shields.io/badge/Nextflow-25.10.0%2B-23aa62.svg)](https://www.nextflow.io/)

sc-pcQTL is a containerized Nextflow workflow for detecting local
co-expression clusters from single-cell RNA-seq data and mapping their
principal-component phenotypes as cis-pcQTLs. It covers the reusable method
from cell-type-specific expression tables through SAIGE-QTL summary
statistics.

## Workflow

1. Exclude cell types below a configurable cell-count threshold.
2. Filter genes by their nonzero-cell fraction.
3. Test local gene-gene associations with Poisson-binomial hurdle models.
4. Call non-overlapping local clusters with a greedy sliding window.
5. Build centered, unscaled cluster-PC phenotypes.
6. Test cluster-PC cis associations with SAIGE-QTL and report within-phenotype
   Benjamini-Hochberg results.

Fine-mapping, colocalization, SMR, S-LDSC, simulations, enrichment analyses,
and manuscript figures are intentionally outside this software repository.

## Requirements

The host needs:

- Linux x86_64;
- Java 17 or newer;
- Nextflow 25.10.0 or newer;
- Docker, Apptainer, or Singularity for analysis dependencies;
- Slurm only when using the `slurm` execution profile.

Java and Nextflow alone are not sufficient because the workflow executes
pinned R/fasthurdle and SAIGE-QTL containers. See
[installation instructions](docs/installation.md).

## Run the example

The bundled semi-synthetic example uses public 1000 Genomes EUR genotypes and
simulated single-cell expression. It runs from input validation through
SAIGE-QTL:

```bash
nextflow run joookerz/sc-pcQTL \
  -profile apptainer,example \
  --outdir results/example_component_union
```

Run the same data with the joint score pair test:

```bash
nextflow run joookerz/sc-pcQTL \
  -profile apptainer,example \
  --pair_test joint_score \
  --outdir results/example_joint_score
```

See [the example documentation](examples/README.md) for its design and
expected outputs. The dataset demonstrates execution and is not a performance
benchmark.

## Run your data

```bash
nextflow run joookerz/sc-pcQTL \
  -profile apptainer,slurm \
  --input samplesheet.csv \
  --gene_annotation genes.tsv \
  --genotype_prefix '/data/genotype_chr{chr}' \
  --variance_ratio_prefix /data/pruned_markers \
  --outdir results \
  -resume
```

Use `-profile docker` on a Docker workstation. Run
`nextflow run joookerz/sc-pcQTL --help` for a concise command summary.

## Statistical Modes

The manuscript-compatible defaults are `--pair_scope fast` and
`--pair_test component_union`. `fast` evaluates non-overlapping genomic blocks
of at most 50 genes. `complete` guarantees coverage of all pairs that could
occur together in a cluster of the configured maximum size, including pairs
across fast-block boundaries, and can evaluate additional nearby pairs.

`component_union` retains the smaller directional p-value separately for the
count and detection components and calls a pair when either component passes
`0.05/M`. `joint_score` uses the two-degree-of-freedom sum of count and
detection score statistics for each ordered direction and uses `0.05/(2M)`.
Both modes use the same global denominator `M`, the number of all unordered
within-chromosome autosomal pairs among filtered genes.

## Documentation

- [Installation](docs/installation.md)
- [Inputs and execution](docs/usage.md)
- [Parameters](docs/parameters.md)
- [Method definitions](docs/methods.md)
- [SAIGE-QTL customization](docs/saigeqtl.md)
- [Outputs](docs/outputs.md)
- [Troubleshooting](docs/troubleshooting.md)

The redistributable semi-synthetic example is stored under `examples/`.
Smaller fixtures under `tests/fixtures/` are used only for software tests.
