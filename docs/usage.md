# Inputs and execution

## Samplesheet

`--input` is a comma-separated file with one unique row per cell type:

```text
celltype,counts
cd8_nc,/absolute/path/cd8_nc_counts.tsv.gz
nk,/absolute/path/nk_counts.tsv.gz
```

Cell-type identifiers may contain letters, numbers, dots, underscores, and
hyphens. Each `counts` file is a wide, tab-separated table with one cell per
row. Relative count paths are resolved against the samplesheet directory. The
manuscript defaults expect these columns:

```text
individual barcode age sex pc1 pc2 pc3 pc4 pc5 pc6 pf1 pf2
total_read_counts log_total_read_counts GENE1 GENE2 ...
```

`individual` is the donor identifier and `barcode` is a unique cell
identifier. Gene columns must be finite and non-negative. The manuscript
analysis used SCTransform-corrected OneK1K counts; applying the workflow to a
different expression scale requires independent calibration.

Configured covariates are donor-level values and must be constant across all
cells from the same donor within a cell type. Library-size columns are
cell-level and are checked separately.

Column names and the covariate list are configurable. Categorical covariates
must also be listed in `--categorical_covariates`.

### Library size

Supplying `total_read_counts` and `log_total_read_counts` is recommended. The
count component uses log library size as an offset, while the detection
component uses it as an ordinary covariate. If either column is absent, the
workflow computes the missing value from expression columns matched to the
gene annotation and records the source in `celltype_qc.tsv`. If the annotation
does not cover the complete transcriptome, provide precomputed library sizes
instead.

## Gene annotation

`--gene_annotation` is a TSV with one unique row per gene and default columns:

```text
gene_name chromosome start end
GIMAP1    7          150417177 150423199
```

Chromosomes may be `1` or `chr1`. Version 0.1 analyzes autosomes only. The
annotation and genotype files must use the same genome build; sc-pcQTL does
not perform liftOver.

The complete autosomal annotation defines genomic adjacency for cluster
calling. Pairwise models are fitted only for genes passing the configured
expression filter. Cluster PCA reads any retained cluster gene directly from
the count table, including a sparsely expressed gene retained within an
otherwise dense window.

## Genotypes

`--genotype_prefix` is a PLINK1 prefix template containing `{chr}`. It must
resolve to BED, BIM, and FAM files for chromosomes 1 through 22, for example:

```text
/data/genotype_chr1.bed
/data/genotype_chr1.bim
/data/genotype_chr1.fam
...
```

Donor IDs in the phenotype table must map to genotype sample IDs. Extra
genotype samples are allowed. Use `--run_qtl false` to stop after cluster-PC
phenotype construction; genotype input is then optional.

## Variance-ratio markers

For exact reuse, provide a pruned PLINK prefix with
`--variance_ratio_prefix`. If omitted, the workflow merges the chromosome
files, filters and LD-prunes variants, and deterministically selects at most
10,000 markers. The generated PLINK files and SHA256 checksums are published
under `pipeline_info/variance_ratio/`.

## Execution profiles

Bundled example with Apptainer:

```bash
nextflow run joookerz/sc-pcQTL -profile apptainer,example \
  --outdir results/example_component_union
```

Docker workstation:

```bash
nextflow run joookerz/sc-pcQTL -profile docker [parameters]
```

Slurm with Apptainer:

```bash
nextflow run joookerz/sc-pcQTL -profile apptainer,slurm [parameters]
```

Add `-resume` after an interruption. Site-specific resources, queues, and
container cache paths belong in a separate Nextflow config:

```bash
nextflow run joookerz/sc-pcQTL -profile apptainer,slurm \
  -c institutional.config [parameters]
```

Do not edit the workflow's committed profiles for a single installation.
