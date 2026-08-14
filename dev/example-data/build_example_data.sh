#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROJECT_DIR=$(cd "${SCRIPT_DIR}/../.." && pwd)
WORK_DIR=${WORK_DIR:-"${SCRIPT_DIR}/work"}
OUT_DIR=${OUT_DIR:-"${PROJECT_DIR}/examples"}
SEED=${SEED:-20260814}
N_CELLS=${N_CELLS:-12000}
SOURCE_BASE=${SOURCE_BASE:-"https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV"}
PANEL_URL=${PANEL_URL:-"https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel"}
BCFTOOLS_IMAGE=${BCFTOOLS_IMAGE:-"docker://quay.io/biocontainers/bcftools:1.20--h8b25389_0"}

find_executable() {
  local name=$1
  shift
  if command -v "${name}" >/dev/null 2>&1; then
    command -v "${name}"
    return
  fi
  local candidate
  for candidate in "$@"; do
    if [[ -x "${candidate}" ]]; then
      printf '%s\n' "${candidate}"
      return
    fi
  done
  return 1
}

PLINK=$(find_executable plink "/PHShome/jz1136/apps/plink/plink") || {
  echo "PLINK 1.9 was not found; set PLINK or add plink to PATH" >&2
  exit 2
}
RSCRIPT=$(find_executable Rscript "/apps/software/R/4.5.1_rbase/bin/Rscript") || {
  echo "Rscript was not found" >&2
  exit 2
}
CURL=$(find_executable curl "/usr/bin/curl") || {
  echo "curl was not found" >&2
  exit 2
}

if command -v bcftools >/dev/null 2>&1; then
  BCFTOOLS=("$(command -v bcftools)")
else
  APPTAINER=$(find_executable apptainer "/apps/software/Apptainer/1.4.2-1.el9/bin/apptainer") || {
    echo "Neither bcftools nor Apptainer was found" >&2
    exit 2
  }
  BCFTOOLS=("${APPTAINER}" exec "${BCFTOOLS_IMAGE}" bcftools)
fi

mkdir -p "${WORK_DIR}/vcf" "${WORK_DIR}/plink" "${WORK_DIR}/selection" \
  "${WORK_DIR}/final" "${OUT_DIR}"
"${CURL}" --fail --location --retry 4 --output "${WORK_DIR}/integrated_call_samples.panel" "${PANEL_URL}"
awk -F'\t' 'NR > 1 && $3 == "EUR" {print $1}' "${WORK_DIR}/integrated_call_samples.panel" > "${WORK_DIR}/eur_samples.txt"
if [[ $(wc -l < "${WORK_DIR}/eur_samples.txt") -ne 503 ]]; then
  echo "The 1000 Genomes panel did not yield exactly 503 EUR samples" >&2
  exit 2
fi

declare -A CHR_LENGTH=(
  [1]=248956422 [2]=242193529 [3]=198295559 [4]=190214555 [5]=181538259
  [6]=170805979 [7]=159345973 [8]=145138636 [9]=138394717 [10]=133797422
  [11]=135086622 [12]=133275309 [13]=114364328 [14]=107043718 [15]=101991189
  [16]=90338345 [17]=83257441 [18]=80373285 [19]=58617616 [20]=64444167
  [21]=46709983 [22]=50818468
)

cat > "${WORK_DIR}/signal_regions.tsv" <<'EOF'
module	center	start	end	genetic
C1_count	18000000	17450000	18550000	1
C2_detection	25000000	24450000	25550000	1
C3_joint	32000000	31450000	32550000	1
C4_mixed	39000000	38450000	39550000	1
C5_null	46000000	45450000	46550000	0
EOF

for chr in $(seq 1 22); do
  length=${CHR_LENGTH[${chr}]}
  if [[ ${chr} -eq 22 ]]; then
    background_center=21500000
  else
    background_center=$((length * 65 / 100))
  fi
  background_start=$((background_center - 2500000))
  background_end=$((background_center + 2500000))
  if [[ ${chr} -eq 22 ]]; then
    regions="chr22:${background_start}-${background_end},chr22:17450000-18550000,chr22:24450000-25550000,chr22:31450000-32550000,chr22:38450000-39550000,chr22:45450000-46550000"
  else
    regions="chr${chr}:${background_start}-${background_end}"
  fi

  vcf="${WORK_DIR}/vcf/chr${chr}.vcf.gz"
  source_vcf="${SOURCE_BASE}/1kGP_high_coverage_Illumina.chr${chr}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
  if [[ ! -s "${vcf}" ]]; then
    echo "Extracting 1000 Genomes chromosome ${chr}"
    (
      cd "${WORK_DIR}/vcf"
      "${BCFTOOLS[@]}" view \
        --samples-file "${WORK_DIR}/eur_samples.txt" \
        --regions "${regions}" \
        --min-alleles 2 --max-alleles 2 --types snps \
        --include 'INFO/MAF_EUR_unrel>=0.05 && INFO/MAF_EUR_unrel<=0.5' \
        --output-type z --output "${vcf}" "${source_vcf}"
    )
  fi

  full_prefix="${WORK_DIR}/plink/chr${chr}_full"
  "${PLINK}" --vcf "${vcf}" --double-id --keep-allele-order --allow-no-sex \
    --geno 0.01 --maf 0.05 --make-bed --out "${full_prefix}" \
    > "${WORK_DIR}/plink/chr${chr}_convert.log" 2>&1
  if [[ $(wc -l < "${full_prefix}.fam") -ne 503 ]]; then
    echo "Chromosome ${chr} does not contain 503 donors" >&2
    exit 2
  fi

  awk -v start="${background_start}" -v end="${background_end}" \
    '$4 >= start && $4 <= end {print $2}' "${full_prefix}.bim" > "${WORK_DIR}/selection/chr${chr}_background_candidates.txt"
  "${PLINK}" --bfile "${full_prefix}" \
    --extract "${WORK_DIR}/selection/chr${chr}_background_candidates.txt" \
    --indep-pairwise 200 50 0.2 --seed "${SEED}" --allow-no-sex \
    --out "${WORK_DIR}/selection/chr${chr}_background_prune" \
    > "${WORK_DIR}/selection/chr${chr}_prune.log" 2>&1
  "${RSCRIPT}" "${SCRIPT_DIR}/select_example_variants.R" \
    --mode background \
    --bim "${full_prefix}.bim" \
    --prune "${WORK_DIR}/selection/chr${chr}_background_prune.prune.in" \
    --n 160 \
    --output "${WORK_DIR}/selection/chr${chr}_background.txt"

  if [[ ${chr} -eq 22 ]]; then
    "${PLINK}" --bfile "${full_prefix}" --freq --missing --allow-no-sex \
      --out "${WORK_DIR}/selection/chr22_signal_qc" \
      > "${WORK_DIR}/selection/chr22_signal_qc.log" 2>&1
    "${RSCRIPT}" "${SCRIPT_DIR}/select_example_variants.R" \
      --mode signals \
      --bim "${full_prefix}.bim" \
      --freq "${WORK_DIR}/selection/chr22_signal_qc.frq" \
      --missing "${WORK_DIR}/selection/chr22_signal_qc.lmiss" \
      --regions "${WORK_DIR}/signal_regions.tsv" \
      --n-per-locus 600 --near-n 200 \
      --output "${WORK_DIR}/selection/chr22_signals.txt" \
      --causal-output "${WORK_DIR}/causal_variants.tsv"
    awk '!seen[$0]++' \
      "${WORK_DIR}/selection/chr22_background.txt" \
      "${WORK_DIR}/selection/chr22_signals.txt" \
      > "${WORK_DIR}/selection/chr22_selected.txt"
  else
    cp "${WORK_DIR}/selection/chr${chr}_background.txt" "${WORK_DIR}/selection/chr${chr}_selected.txt"
  fi

  output_prefix="${WORK_DIR}/final/genotype_chr${chr}"
  "${PLINK}" --bfile "${full_prefix}" \
    --extract "${WORK_DIR}/selection/chr${chr}_selected.txt" \
    --keep-allele-order --allow-no-sex --make-bed --out "${output_prefix}" \
    > "${WORK_DIR}/plink/chr${chr}_final.log" 2>&1
  cp "${output_prefix}.bed" "${OUT_DIR}/genotype_chr${chr}.bed"
  cp "${output_prefix}.bim" "${OUT_DIR}/genotype_chr${chr}.bim"
  cp "${output_prefix}.fam" "${OUT_DIR}/genotype_chr${chr}.fam"
done

reference_fam="${OUT_DIR}/genotype_chr1.fam"
for chr in $(seq 2 22); do
  cmp -s "${reference_fam}" "${OUT_DIR}/genotype_chr${chr}.fam" || {
    echo "Final PLINK FAM files differ between chromosomes 1 and ${chr}" >&2
    exit 2
  }
done

prefix_list="${WORK_DIR}/genotype_prefixes.txt"
for chr in $(seq 1 22); do
  printf '%s\n' "${OUT_DIR}/genotype_chr${chr}"
done > "${prefix_list}"

export PATH="$(dirname "${PLINK}"):${PATH}"
bash "${PROJECT_DIR}/bin/build_variance_ratio_markers.sh" \
  "${prefix_list}" "${WORK_DIR}/variance_ratio" 0.05 0.01 200 50 0.2 3000 "${SEED}" \
  > "${WORK_DIR}/variance_ratio_build.log" 2>&1
cp "${WORK_DIR}/variance_ratio.bed" "${OUT_DIR}/variance_ratio.bed"
cp "${WORK_DIR}/variance_ratio.bim" "${OUT_DIR}/variance_ratio.bim"
cp "${WORK_DIR}/variance_ratio.fam" "${OUT_DIR}/variance_ratio.fam"

"${PLINK}" --bfile "${WORK_DIR}/variance_ratio" --pca 6 header tabs --allow-no-sex \
  --out "${WORK_DIR}/genotype_pcs" > "${WORK_DIR}/genotype_pcs.log" 2>&1
cp "${WORK_DIR}/genotype_pcs.eigenvec" "${OUT_DIR}/genotype_pcs.tsv"

awk 'NR > 1 {print $2}' "${WORK_DIR}/causal_variants.tsv" > "${WORK_DIR}/causal_variants.txt"
"${PLINK}" --bfile "${OUT_DIR}/genotype_chr22" \
  --extract "${WORK_DIR}/causal_variants.txt" --recode A --keep-allele-order --allow-no-sex \
  --out "${WORK_DIR}/causal_dosage" > "${WORK_DIR}/causal_dosage.log" 2>&1

"${RSCRIPT}" "${SCRIPT_DIR}/simulate_expression.R" \
  --fam "${OUT_DIR}/genotype_chr22.fam" \
  --eigenvec "${WORK_DIR}/genotype_pcs.eigenvec" \
  --dosage "${WORK_DIR}/causal_dosage.raw" \
  --causal "${WORK_DIR}/causal_variants.tsv" \
  --outdir "${OUT_DIR}" \
  --n-cells "${N_CELLS}" \
  --seed "${SEED}"

echo "Example dataset generated in ${OUT_DIR}"
