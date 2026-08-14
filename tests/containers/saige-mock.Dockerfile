ARG BASE_IMAGE=sc-pcqtl-core:test
FROM ${BASE_IMAGE}

COPY tests/mocks/step1_fitNULLGLMM_qtl.R /usr/local/bin/
COPY tests/mocks/step2_tests_qtl.R /usr/local/bin/
COPY tests/mocks/step3_gene_pvalue_qtl.R /usr/local/bin/

RUN chmod 0755 /usr/local/bin/step*_qtl.R
