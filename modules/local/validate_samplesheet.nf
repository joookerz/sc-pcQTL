process VALIDATE_SAMPLESHEET {
    tag 'samplesheet'
    label 'process_low'
    publishDir "${params.outdir}/pipeline_info", mode: 'copy', overwrite: true

    input:
    path samplesheet
    val base_dir
    path workflow_bin

    output:
    path 'validated_samplesheet.csv', emit: samplesheet

    script:
    """
    Rscript ${workflow_bin}/validate_samplesheet.R \
      --input '${samplesheet}' \
      --base_dir '${base_dir}' \
      --out validated_samplesheet.csv
    """
}
