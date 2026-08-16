# Method definitions

## Hurdle pair screen

For response gene `Y`, predictor gene `X`, covariates `Z`, and log library
size `L`, the default positive-count model is

```text
log E(Y | Y > 0) = beta0 + beta1 X + Z beta + offset(L),
```

and the detection model is

```text
logit P(Y > 0) = gamma0 + gamma1 X + Z gamma + gamma_lib L.
```

Both ordered directions are fitted. In `component_union`, the smaller
directional p-value is retained separately for each component, and a pair is
significant if either is below `pair_alpha/M`. In `joint_score`, count and
detection score statistics are summed within each ordered direction and
evaluated against chi-squared with two degrees of freedom. An ordered test is
significant below `pair_alpha/(2M)`, and an unordered pair is retained if
either direction passes.

## Pair scopes

`fast` partitions filtered genes on each chromosome into non-overlapping
blocks of no more than `max_cluster_genes` and tests all pairs within each
block. This is the default block-scheduling strategy.

`complete` tests every pair separated by at most
`max_cluster_genes - 1` ranks among filtered genes. This guarantees coverage
of every filtered-gene pair that can coexist in a candidate sliding window and
can include additional pairs when filtered-out annotated genes lie between
them. It is more complete but more expensive.

## Cluster calling

All autosomal genes in the supplied annotation are ordered by genomic start,
matching the production workflow; significant edges are available only for
genes that passed the expression filter. Windows are evaluated from
`max_cluster_genes` down to `min_cluster_genes`. A window is retained when the
fraction of significant unordered pairwise edges reaches `cluster_density`.
Its genes are marked assigned, yielding non-overlapping, largest-first local
clusters.

## Cluster PCs and cis-QTL mapping

PCA is performed per cluster on centered, unscaled expression values. The
minimum number of leading PCs reaching `pca_variance` cumulative variance is
tested. SAIGE-QTL uses quantitative-trait mode, inverse normalization,
donor-level covariates repeated across cells, cluster-span plus/minus the cis
window, and the configured MAF threshold. Variant p-values are BH-adjusted
within each cluster-PC phenotype.
