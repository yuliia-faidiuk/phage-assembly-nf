process MULTIQC {
    label 'process_low'
    container params.containers.multiqc
    publishDir path: { "${params.outdir}/multiqc" }, mode: 'copy'

    input:
    path files, stageAs: 'inputs/*'

    output:
    path "multiqc_report.html"
    path "multiqc_data"

    script:
    """
    export HOME=\$PWD
    multiqc -f -o . inputs > multiqc.log 2>&1
    """
}

process SUMMARY {
    // One table with every number the homework asks for, per sample, plus a readable Markdown version.
    label 'process_low'
    container params.containers.multiqc   // has python3
    publishDir path: { "${params.outdir}" }, mode: 'copy'

    input:
    path summaries, stageAs: 'summaries/*'
    path fastp_json, stageAs: 'fastp/*'
    path samplesheet

    output:
    path "summary.tsv"
    path "summary.md"

    script:
    """
    summarize.py --summaries summaries --fastp fastp --samplesheet $samplesheet --tsv summary.tsv --md summary.md ${params.skip_phageterm ? "--skip-phageterm" : ""}
    """
}
