process COLLECT_QTL_RESULTS {
    tag 'QTL summary'
    label 'process_medium'
    publishDir "${params.outdir}/summary", mode: 'copy', overwrite: true,
        saveAs: { name -> name.startsWith('qtl_summary/') ? name.substring(12) : name }

    input:
    path result_dirs
    path workflow_bin

    output:
    path 'qtl_summary/*', emit: summary

    script:
    def inputLines = result_dirs.collect { resultDir -> resultDir.toString() }.join('\n')
    """
    printf '%s\n' '${inputLines}' > qtl_inputs.txt
    Rscript ${workflow_bin}/collect_qtl_results.R \
      --input_list qtl_inputs.txt \
      --outdir qtl_summary \
      --qtl_fdr '${params.qtl_fdr}'
    """
}
