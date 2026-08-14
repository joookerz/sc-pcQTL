#!/usr/bin/env bash
set -euo pipefail

[[ $# -eq 9 ]] || {
  echo "Usage: build_variance_ratio_markers.sh <prefix-list> <out-prefix> <maf> <geno> <window> <step> <r2> <n-markers> <seed>" >&2
  exit 2
}
prefix_list=$1
out_prefix=$2
maf=$3
geno=$4
window=$5
step=$6
r2=$7
n_markers=$8
seed=$9

mkdir -p "$(dirname "${out_prefix}")"
first_prefix=$(head -n1 "${prefix_list}")
[[ -n "${first_prefix}" ]] || { echo "No genotype prefixes were provided" >&2; exit 2; }
while IFS= read -r prefix; do
  [[ -n "${prefix}" ]] || continue
  cmp -s "${first_prefix}.fam" "${prefix}.fam" || {
    echo "PLINK FAM files differ across chromosomes: ${first_prefix}.fam vs ${prefix}.fam" >&2
    exit 2
  }
done < "${prefix_list}"
tail -n +2 "${prefix_list}" | awk '{print $0".bed "$0".bim "$0".fam"}' > "${out_prefix}.merge_list"
plink --bfile "${first_prefix}" --merge-list "${out_prefix}.merge_list" \
  --allow-no-sex --make-bed --out "${out_prefix}.merged"
plink --bfile "${out_prefix}.merged" --autosome --maf "${maf}" --geno "${geno}" \
  --indep-pairwise "${window}" "${step}" "${r2}" --seed "${seed}" \
  --out "${out_prefix}.prune"
test -s "${out_prefix}.prune.prune.in"
awk -v seed="${seed}" 'BEGIN{srand(seed)} {print rand(),$0}' "${out_prefix}.prune.prune.in" \
  | sort -k1,1n | sed -n "1,${n_markers}p" | cut -d' ' -f2- > "${out_prefix}.extract"
plink --bfile "${out_prefix}.merged" --extract "${out_prefix}.extract" \
  --make-bed --allow-no-sex --out "${out_prefix}"
sha256sum "${out_prefix}.bed" "${out_prefix}.bim" "${out_prefix}.fam" > "${out_prefix}.sha256"
printf 'parameter\tvalue\nmaf\t%s\ngeno\t%s\nprune_window\t%s\nprune_step\t%s\nprune_r2\t%s\nmax_markers\t%s\nseed\t%s\n' \
  "${maf}" "${geno}" "${window}" "${step}" "${r2}" "${n_markers}" "${seed}" \
  > "${out_prefix}.settings.tsv"
