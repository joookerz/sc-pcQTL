#!/usr/bin/env bash
set -euo pipefail

usage() {
  echo "Usage: run_saige_qtl.sh <celltype> <cluster> <phenotype> <chr> <pheno> <region> <bed> <bim> <fam> <vr-prefix> <params> <outdir>" >&2
  exit 2
}
[[ $# -eq 12 ]] || usage

celltype=$1
cluster=$2
phenotype=$3
chromosome=$4
pheno=$(realpath "$5")
region=$(realpath "$6")
bed=$(realpath "$7")
bim=$(realpath "$8")
fam=$(realpath "$9")
vr_prefix=$(realpath -m "${10}")
params=$(realpath "${11}")
mkdir -p "${12}"
outdir=$(realpath "${12}")
task_dir=$(dirname "${outdir}")

for path in "${pheno}" "${region}" "${bed}" "${bim}" "${fam}" \
            "${vr_prefix}.bed" "${vr_prefix}.bim" "${vr_prefix}.fam" "${params}"; do
  test -s "${path}" || { echo "Missing or empty required input: ${path}" >&2; exit 2; }
done
cmp -s "${fam}" "${vr_prefix}.fam" || {
  echo "Chromosome and variance-ratio PLINK FAM files differ" >&2
  exit 2
}

args_for_step() {
  local step=$1
  awk -F'\t' -v step="${step}" 'NR>1 && $1==step {printf "--%s=%s\n", $2, $3}' "${params}"
}
mapfile -t step1_args < <(args_for_step step1)
mapfile -t step2_args < <(args_for_step step2)
mapfile -t step3_args < <(args_for_step step3)

# Null-model artifacts are task-local intermediates. Keeping them outside the
# published result directory avoids duplicating large R objects per phenotype.
prefix="${task_dir}/saige_null_model"
association="${outdir}/association.tsv"
acat="${outdir}/acat.tsv"

step1_fitNULLGLMM_qtl.R \
  "${step1_args[@]}" \
  --phenoFile="${pheno}" \
  --phenoCol="${phenotype}" \
  --sampleIDColinphenoFile=individual \
  --outputPrefix="${prefix}" \
  --plinkFile="${vr_prefix}"

test -s "${prefix}.rda"
test -s "${prefix}.varianceRatio.txt"

step2_tests_qtl.R \
  "${step2_args[@]}" \
  --bedFile="${bed}" \
  --bimFile="${bim}" \
  --famFile="${fam}" \
  --SAIGEOutputFile="${association}" \
  --chrom="${chromosome}" \
  --GMMATmodelFile="${prefix}.rda" \
  --varianceRatioFile="${prefix}.varianceRatio.txt" \
  --rangestoIncludeFile="${region}"

test -s "${association}"
step3_gene_pvalue_qtl.R \
  "${step3_args[@]}" \
  --assocFile="${association}" \
  --geneName="${phenotype}" \
  --genePval_outputFile="${acat}"

printf 'celltype\tcluster_id\tphenotype_id\tchromosome\n%s\t%s\t%s\t%s\n' \
  "${celltype}" "${cluster}" "${phenotype}" "${chromosome}" > "${outdir}/metadata.tsv"
printf 'OK\n' > "${outdir}/COMPLETE"
