# Semi-synthetic example

This directory contains a complete, redistributable sc-pcQTL input dataset.
It combines public phased genotypes from the 1000 Genomes Project 30x GRCh38
release with simulated covariates and single-cell expression. The example is
designed to exercise every workflow stage through SAIGE-QTL; it is not a
method-performance benchmark.

## Design

- One simulated immune cell type, `sim_immune`.
- 503 unrelated European-ancestry 1000 Genomes donors.
- 12,000 simulated cells and 40 simulated genes on chromosome 22.
- Five local co-expression modules separated by more than 5 Mb: three
  two-gene modules, one four-gene module, and one five-gene module.
- Four modules carry an implanted cis genetic effect; the fifth is a
  co-expressed negative control without a genetic effect.
- Detection is generated separately from zero-truncated Poisson positive
  counts, with library size, donor covariates, donor effects, cell-level
  factors, and module-specific genotype effects.

The genotype source is the 1000 Genomes Project high-coverage 30x phased
GRCh38 call set described by Byrska-Bishop et al. (Cell, 2022). Chromosomes
1-21 contain compact LD-pruned background panels. Chromosome 22 additionally
contains dense local panels around the five simulated modules. The original
sample identifiers and alleles are retained.

## Files

| Path | Content |
|---|---|
| `samplesheet.csv` | sc-pcQTL samplesheet for `sim_immune` |
| `counts.tsv.gz` | Cell-by-gene count table with donor and cell covariates |
| `gene_annotation.tsv` | GRCh38 coordinates for 40 simulated genes |
| `genotype_chr{1..22}.{bed,bim,fam}` | Chromosome-specific PLINK1 genotypes |
| `variance_ratio.{bed,bim,fam}` | Precomputed LD-pruned SAIGE-QTL marker set |
| `genotype_pcs.tsv` | Genotype principal components used during simulation |
| `EXPECTED_RESULTS.md` | Expected workflow-level checks for both pair tests |

`sex` is encoded as `0/1`. The workflow treats it as categorical in the
hurdle screen and passes its numeric encoding to SAIGE-QTL. All FAM files use
the same 503 donors and ordering.

## Run

From a cloned repository, run the manuscript-default component-union mode:

```bash
nextflow run . -profile apptainer,example \
  --outdir results/example_component_union
```

Run the same input with the joint score test:

```bash
nextflow run . -profile apptainer,example \
  --pair_test joint_score \
  --outdir results/example_joint_score
```

Use `-profile docker,example` on a Docker workstation. Add `-resume` to reuse
completed tasks after interruption.

The maintainer-only data generator is isolated under `dev/example-data/` and
is not called during a normal workflow run.
