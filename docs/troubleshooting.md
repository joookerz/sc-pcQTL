# Troubleshooting

## Container is unavailable

The `edge` image is built from `main`. Check GHCR visibility and test a pull
with the selected runtime. For restricted clusters, pre-pull Apptainer images
into a shared cache and set `NXF_APPTAINER_CACHEDIR`.

## Pair jobs wait in the scheduler

Each pair task requests one CPU. Adjust memory/time or the queue by supplying
a site config. Lower `pair_responses_per_task` to create shorter `complete`
tasks; increase it only after measuring memory and wall time.

## Library sizes were computed unexpectedly

Check the exact names configured by `total_library_col` and
`log_library_col`. The QC table records whether values came from input or were
computed.

## SAIGE-QTL reports missing samples

Confirm that phenotype donor IDs and PLINK FAM sample IDs use the same
identifier, and that all chromosome FAM files contain the same samples in the
same order.

## Resume after failure

Correct the input or site configuration and rerun the identical command with
`-resume`. Do not delete the Nextflow work directory before resuming.
