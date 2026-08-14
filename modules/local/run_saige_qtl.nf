process RUN_SAIGE_QTL {
    tag "$celltype:$cluster_id:$phenotype_id"
    label 'process_qtl'
    publishDir "${params.outdir}/qtl/tasks", mode: 'copy', overwrite: true

    input:
    tuple val(celltype), val(task_id), val(cluster_id), val(phenotype_id), val(chromosome),
          path(phenotype_file), path(region_file), path(bed), path(bim), path(fam),
          path(vr_bed), path(vr_bim), path(vr_fam), path(saige_params)
    path workflow_bin

    output:
    tuple val(celltype), path("${celltype}__${cluster_id}__${phenotype_id}"), emit: result

    script:
    def vrPrefix = vr_bed.baseName
    """
    bash ${workflow_bin}/run_saige_qtl.sh \
      '${celltype}' '${cluster_id}' '${phenotype_id}' '${chromosome}' \
      '${phenotype_file}' '${region_file}' '${bed}' '${bim}' '${fam}' \
      '${vrPrefix}' '${saige_params}' '${celltype}__${cluster_id}__${phenotype_id}'
    """
}
