#!/usr/bin/env Rscript

Sys.setenv(SCPCQTL_ROOT = normalizePath("."))
testthat::test_dir("tests/testthat", reporter = "summary")
