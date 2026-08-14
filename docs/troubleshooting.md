# Troubleshooting

## Environment diagnosis

Run `sc-pcqtl doctor` first. It checks Java, Nextflow, the selected container
runtime, and whether the Docker/Podman service is running. Select a runtime
explicitly when more than one is installed:

```bash
SCPCQTL_RUNTIME=apptainer sc-pcqtl doctor
```

## Docker or Podman is installed but unavailable

Start Docker Desktop or the Podman machine on macOS. On Linux, start the
service and verify that the current user has permission to run the runtime.
The launcher tests the daemon with `docker info` or `podman info`.

## Apple Silicon image warning

The workflow intentionally runs the pinned `linux/amd64` analysis images on
Apple Silicon. An architecture warning is expected. A failure usually means
Docker/Podman emulation is disabled or the virtual machine lacks sufficient
memory; allocate at least 8 GB for the bundled example and more for production
data.

## Container is unavailable

The `edge` image is built from `main`. Check GHCR visibility and test a pull
with the selected runtime. For restricted clusters, pre-pull Apptainer images
into a shared cache and set `NXF_APPTAINER_CACHEDIR`.

If GHCR authentication is required, log in before running the workflow:

```bash
printf '%s' "$GITHUB_TOKEN" | docker login ghcr.io -u USERNAME --password-stdin
```

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
