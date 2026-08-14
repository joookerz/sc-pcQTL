# SAIGE-QTL configuration

The pinned SAIGE-QTL image is based on version 0.3.4. Built-in defaults are in
`assets/saigeqtl_defaults.tsv` and reproduce the manuscript analysis.

## Override table

Pass optional settings with `--saige_params overrides.tsv`:

```text
step	parameter	value
step1	maxiter	30
step2	markers_per_chunk	20000
```

`step` must be `step1`, `step2`, or `step3`. The workflow validates duplicate
and unknown parameters. Set `--allow_unknown_saige_params true` only when a
custom compatible image adds a parameter not recognized by this release.

## Workflow-owned arguments

File paths, output prefixes, phenotype/sample identifiers, chromosome,
region, quantitative trait type, covariate lists, and minimum MAF cannot be
overridden in the table. Use the corresponding sc-pcQTL parameter instead,
for example `--covariates` or `--qtl_maf`. This prevents a parameter table from
silently disconnecting SAIGE-QTL from workflow-managed inputs.

The fully resolved table is always written to
`pipeline_info/resolved_saigeqtl_params.tsv`.
