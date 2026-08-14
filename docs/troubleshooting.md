# Troubleshooting

## Environment diagnosis

Run `sc-pcqtl doctor` first. It checks Java, Nextflow, the selected container
runtime, and whether the Docker/Podman service is running. Select a runtime
explicitly when more than one is installed:

```bash
SCPCQTL_RUNTIME=apptainer sc-pcqtl doctor
```

If the basic check succeeds but a workflow task still cannot start, run:

```bash
sc-pcqtl doctor --deep
```

The deep check pulls a small public image and verifies amd64 execution, a
writable host bind mount, and the resource-limit arguments used with
Docker/Podman. It does not pull the much larger sc-pcQTL analysis images.
On a restricted network, set `SCPCQTL_SMOKE_IMAGE` to an accessible small
Linux image mirrored by the institution.

## Docker or Podman is installed but unavailable

Start Docker Desktop or the Podman machine on macOS. On Linux, start the
service and verify that the current user has permission to run the runtime.
The launcher tests the daemon with `docker info` or `podman info`.

Rootless Podman additionally requires subordinate UID/GID mappings and cgroup
delegation. Ask the system administrator to configure these when
`sc-pcqtl doctor --deep` fails despite a successful `podman info`; otherwise,
use Docker or Apptainer. Do not remove Nextflow resource controls merely to
work around an incomplete rootless Podman setup.

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

Repository visibility and GHCR package visibility are separate settings. A
public repository can still reference a private package, so confirm each
analysis image is public or authenticate the runtime explicitly.

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
`-resume`. Do not delete the Nextflow work directory before resuming. Workflow
reports are overwritten at their stable paths under `pipeline_info/`; cached
task outputs remain controlled by the Nextflow work directory.
