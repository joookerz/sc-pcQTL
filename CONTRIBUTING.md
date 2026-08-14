# Contributing

Open an issue before changing statistical defaults or output schemas. Pull
requests should include a focused rationale, tests for changed behavior, and
documentation updates.

Run the local checks from the repository root:

```bash
nextflow lint main.nf
Rscript tests/testthat.R
tests/run_core_integration.sh
tests/run_qtl_mock_integration.sh
```

Manual edits should not mix scientific-default changes with refactoring in a
single commit. Container dependencies must be pinned to a release or commit.
