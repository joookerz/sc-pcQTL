process PLAN_PAIR_TASKS {
    tag "$celltype"
    label 'process_low'

    input:
    tuple val(celltype), path(stage)
    path workflow_bin

    output:
    tuple val(celltype), path("${celltype}_pair_tasks.tsv"), emit: tasks
    tuple val(celltype), path("${celltype}_pair_tasks_summary.tsv"), emit: summary

    script:
    """
    Rscript ${workflow_bin}/plan_pair_tasks.R \
      --stage '${stage}' \
      --out '${celltype}_pair_tasks.tsv' \
      --pair_scope '${params.pair_scope}' \
      --max_cluster_genes '${params.max_cluster_genes}' \
      --responses_per_task '${params.pair_responses_per_task}'
    """
}
