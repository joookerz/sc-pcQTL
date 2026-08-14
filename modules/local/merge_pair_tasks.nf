process MERGE_PAIR_TASKS {
    tag "$celltype"
    label 'process_medium'
    publishDir "${params.outdir}/pairs", mode: 'copy', overwrite: true,
        saveAs: { name -> name == "${celltype}_pairs" ? celltype : name }

    input:
    tuple val(celltype), path(pair_files), path(stage)
    path workflow_bin

    output:
    tuple val(celltype), path("${celltype}_pairs"), emit: pairs

    script:
    def inputLines = pair_files.collect { pairFile -> pairFile.toString() }.join('\n')
    """
    printf '%s\n' '${inputLines}' > pair_inputs.txt
    Rscript ${workflow_bin}/merge_pair_tasks.R \
      --stage '${stage}' \
      --input_list pair_inputs.txt \
      --outdir '${celltype}_pairs' \
      --pair_scope '${params.pair_scope}' \
      --pair_test '${params.pair_test}' \
      --count_family '${params.count_family}' \
      --pair_alpha '${params.pair_alpha}' \
      --save_directional_tests '${params.save_directional_tests}'
    """
}
