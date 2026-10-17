#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { PREPARE_READS; FASTQC as FASTQC_RAW; FASTQC as FASTQC_TRIMMED; FASTP; HOST_REMOVE; CAP_READS; QUICK_ASSEMBLY; SUBSAMPLE } from './modules/reads.nf'
include { UNICYCLER; SPADES; BANDAGE_IMAGE; CONTIG_TABLE } from './modules/assembly.nf'
include { MAP_READS; VARIANTS; PILON; PHAGETERM; ROTATE } from './modules/finishing.nf'
include { MULTIQC; SUMMARY } from './modules/report.nf'

def helpMessage() {
    log.info """
    phage-assembly-nf  —  Illumina paired-end phage genome assembly (mirrors phage_genome_homework.html)

    Usage (simplest):   bash phage-assembly-nf/run.sh <folder_with_fastq.gz>
    Usage (nextflow):   nextflow run phage-assembly-nf/ --input <folder> [--genome_size 40000] [--outdir results] -resume

      --input        folder with paired FASTQ files: <sample>_R1*.fastq.gz / <sample>_R2*.fastq.gz,
                     Illumina facility names (phageX_S1_L001_R1_001.fastq.gz) and _1/_2 names work; lanes are merged
      --samplesheet  instead of --input: CSV with header  sample,fastq_1,fastq_2
      --outdir       output folder (default: results)
      --genome_size  expected genome size in bp; if omitted it is estimated from a quick assembly
      --coverage     subsample target (default 200); --fallback_coverage (100) is used if the assembly fails
      --cap_pairs    maximum read pairs ever given to the assembler (default 500000)
      --unicycler_modes  comma list of normal,conservative,bold (default: all three)
      --run_spades   also run SPAdes --isolate (default true)
      --primary      assembly used downstream (default unicycler_normal)
      --host_fasta   optional host genome FASTA; matching read pairs are removed before assembly
      --skip_phageterm   skip PhageTerm
      -profile test  run on the two public phage runs in test_data/
    """.stripIndent()
}

// Same rule as bin/detect_samples.py: strip _S1 / _L001 / _001 decorations, read number is R1/R2 or 1/2.
def parseReadName(String fname) {
    def base = fname.replaceAll(/\.(fastq|fq)\.gz$/, '')
    def m = base =~ /^(.+?)(?:_S\d+)?(?:_L\d{3})?[._]R?([12])(?:_001)?$/
    return m ? [id: m[0][1], read: m[0][2]] : null
}

workflow {
    main:
    if (params.help) { helpMessage(); exit 0 }
    if (!params.input && !params.samplesheet) { helpMessage(); error "Give --input <folder> or --samplesheet <csv>" }

    // ------------------------------------------------------------------ samples
    def ch_samples
    if (params.samplesheet) {
        ch_samples = Channel.fromPath(params.samplesheet, checkIfExists: true)
            .splitCsv(header: true)
            .map { row -> [row.sample.trim(), file(row.fastq_1.trim(), checkIfExists: true), file(row.fastq_2.trim(), checkIfExists: true)] }
            .groupTuple()
            .map { id, r1s, r2s -> [[id: id], r1s.sort { it.name }, r2s.sort { it.name }] }
    } else {
        def indir = file(params.input)
        if (!indir.isDirectory()) error "--input ${params.input} is not a folder"
        ch_samples = Channel.fromPath("${indir}/*.{fastq,fq}.gz")
            .map { f -> def p = parseReadName(f.name); if (!p) log.warn "Ignoring ${f.name}: cannot tell sample name / read number"; p ? [p.id, p.read, f] : null }
            .filter { it != null }
            .groupTuple()
            .map { id, reads, files ->
                def r1 = [files, reads].transpose().findAll { it[1] == '1' }.collect { it[0] }.sort { it.name }
                def r2 = [files, reads].transpose().findAll { it[1] == '2' }.collect { it[0] }.sort { it.name }
                if (r1.size() == 0 || r2.size() == 0) error "Sample '${id}' has R1 files ${r1*.name} but R2 files ${r2*.name}. This pipeline needs paired Illumina reads."
                if (r1.size() != r2.size()) error "Sample '${id}' has ${r1.size()} R1 files but ${r2.size()} R2 files."
                [[id: id], r1, r2]
            }
    }
    ch_samples = ch_samples.ifEmpty { error "No FASTQ files found. Expected *_R1*.fastq.gz / *_R2*.fastq.gz in ${params.input}" }
    ch_samples.map { meta, r1, r2 -> "  ${meta.id}: ${r1*.name.join(',')}  +  ${r2*.name.join(',')}" }
              .collect().view { "\nSamples detected:\n" + it.join('\n') + "\n" }

    ch_sheet = ch_samples.map { meta, r1, r2 -> [r1, r2].transpose().collect { a, b -> "${meta.id},${a},${b}" }.join('\n') }
                         .collectFile(name: 'samplesheet.csv', newLine: true, seed: 'sample,fastq_1,fastq_2', storeDir: "${params.outdir}/pipeline_info")

    // ------------------------------------------------------------------ reads
    PREPARE_READS(ch_samples)
    FASTQC_RAW(PREPARE_READS.out.reads, '01_qc_raw')
    FASTP(PREPARE_READS.out.reads)
    FASTQC_TRIMMED(FASTP.out.reads, '03_qc_trimmed')

    ch_clean = FASTP.out.reads
    ch_host_summary = Channel.empty()
    if (params.host_fasta) {
        HOST_REMOVE(FASTP.out.reads, file(params.host_fasta, checkIfExists: true))
        ch_clean = HOST_REMOVE.out.reads
        ch_host_summary = HOST_REMOVE.out.summary
    }

    CAP_READS(ch_clean)

    ch_gsize_summary = Channel.empty()
    if (params.genome_size) {
        ch_gsize = CAP_READS.out.reads.map { meta, r1, r2, info -> [meta, params.genome_size as Long] }
    } else {
        QUICK_ASSEMBLY(CAP_READS.out.quick)
        ch_gsize = QUICK_ASSEMBLY.out.size.map { meta, f -> [meta, f.text.trim() as Long] }
        ch_gsize_summary = QUICK_ASSEMBLY.out.summary
    }
    SUBSAMPLE(CAP_READS.out.reads.join(ch_gsize))

    // ------------------------------------------------------------------ assemblies
    modes = params.unicycler_modes.toString().split(',').collect { it.trim() }.findAll { it }
    modes.each { if (!(it in ['normal', 'conservative', 'bold'])) error "Unknown Unicycler mode '${it}' (use normal, conservative, bold)" }
    UNICYCLER(SUBSAMPLE.out.reads.combine(Channel.fromList(modes)))
    ch_assemblies = UNICYCLER.out.assembly
    ch_asm_summary = UNICYCLER.out.summary
    if (params.run_spades) {
        SPADES(SUBSAMPLE.out.reads)
        ch_assemblies = ch_assemblies.mix(SPADES.out.assembly)
        ch_asm_summary = ch_asm_summary.mix(SPADES.out.summary)
    }
    BANDAGE_IMAGE(ch_assemblies)
    CONTIG_TABLE(ch_assemblies)

    asm_names = modes.collect { "unicycler_${it}" } + (params.run_spades ? ['spades'] : [])
    if (!(params.primary in asm_names)) error "--primary '${params.primary}' is not among the assemblies that will run: ${asm_names.join(', ')}"
    ch_primary = ch_assemblies.filter { meta, name, dir -> name == params.primary }

    // ------------------------------------------------------------------ finishing
    MAP_READS(ch_clean.join(ch_primary))
    VARIANTS(MAP_READS.out.bam)
    PILON(MAP_READS.out.bam)

    ch_pt_summary = Channel.empty()
    if (!params.skip_phageterm) {
        PHAGETERM(ch_clean.join(PILON.out.fasta))
        ch_pt_summary = PHAGETERM.out.summary
        ch_for_rotate = PILON.out.fasta.join(MAP_READS.out.termini).join(PHAGETERM.out.seq, remainder: true).map { meta, pol, term, pt -> [meta, pol, term, pt ?: []] }
    } else {
        ch_for_rotate = PILON.out.fasta.join(MAP_READS.out.termini).map { meta, pol, term -> [meta, pol, term, []] }
    }
    ROTATE(ch_for_rotate)

    // ------------------------------------------------------------------ reports
    MULTIQC(FASTQC_RAW.out.zip.mix(FASTQC_TRIMMED.out.zip, FASTP.out.json, MAP_READS.out.flagstat, MAP_READS.out.stats).collect())
    ch_all_summaries = PREPARE_READS.out.summary.mix(ch_host_summary, CAP_READS.out.summary, ch_gsize_summary, SUBSAMPLE.out.summary,
                           ch_asm_summary, MAP_READS.out.summary, VARIANTS.out.summary, PILON.out.summary, ch_pt_summary, ROTATE.out.summary)
    SUMMARY(ch_all_summaries.collect(), FASTP.out.json.collect(), ch_sheet)

    onComplete:
        def trace = file("${params.outdir}/pipeline_info/trace.txt")
        def failed = []
        if (trace.exists()) {
            trace.readLines().drop(1).each { line ->
                def c = line.split('\t')
                if (c.size() >= 9 && c[3] == 'FAILED') failed << "  ${c[2]}  (exit ${c[4]})  ->  ${c[8]}"
            }
        }
        def msg = workflow.success ? "Pipeline finished OK." : "Pipeline FAILED."
        log.info """
        ${msg}
          results:   ${params.outdir}/<sample>/01_qc_raw ... 10_final
          summary:   ${params.outdir}/summary.md   (and summary.tsv)
          multiqc:   ${params.outdir}/multiqc/multiqc_report.html
          reports:   ${params.outdir}/pipeline_info/report.html  timeline.html  trace.txt
        """.stripIndent()
        if (failed) {
            log.info "Failed steps (open <workdir>/.command.log for the tool's own error message):\n" + failed.join('\n') +
                     "\n\nSee TROUBLESHOOTING.md in the pipeline folder. Fix the cause, then rerun the same command: finished steps are reused (-resume).\n"
        } else if (!workflow.success) {
            log.info "Error: ${workflow.errorMessage ?: workflow.errorReport}\nSee TROUBLESHOOTING.md in the pipeline folder\n"
        }
}
