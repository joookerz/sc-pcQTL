process RUN_PAIR_TASK {
    tag "$celltype:chr$chromosome:$task_id"
    label 'process_pair'

    input:
    tuple val(celltype), val(task_id), val(chromosome), val(response_start),
          val(response_end), val(block_start), val(block_end), path(stage)
    path workflow_bin

    output:
    tuple val(celltype), val(task_id), path("pair_task_${task_id}.tsv.gz"), emit: result

    script:
    """
    Rscript ${workflow_bin}/run_pair_task.R \
      --stage '${stage}' \
      --out 'pair_task_${task_id}.tsv.gz' \
      --celltype '${celltype}' \
      --chromosome '${chromosome}' \
      --response_start '${response_start}' \
      --response_end '${response_end}' \
      --block_start '${block_start}' \
      --block_end '${block_end}' \
      --pair_scope '${params.pair_scope}' \
      --pair_test '${params.pair_test}' \
      --count_family '${params.count_family}' \
      --max_cluster_genes '${params.max_cluster_genes}' \
      --covariates '${params.covariates}'
    """
}
