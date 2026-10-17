// Step 0-4 of the homework: get the data, QC, trim, QC again, subsample.

process PREPARE_READS {
    // Concatenate lanes (if several files per read) and count raw pairs.
    tag "$meta.id"
    label 'process_low'
    container params.containers.seqkit
    publishDir path: { "${params.outdir}/${meta.id}/00_raw_reads" }, mode: 'copy', pattern: '*.tsv'

    input:
    tuple val(meta), path(r1s, stageAs: 'r1/*'), path(r2s, stageAs: 'r2/*')

    output:
    tuple val(meta), path("${meta.id}_R1.fastq.gz"), path("${meta.id}_R2.fastq.gz"), emit: reads
    path "${meta.id}.raw_stats.tsv"
    path "${meta.id}.00_raw.summary.tsv", emit: summary

    script:
    """
    set -euo pipefail
    # gzip members can simply be concatenated; sorting keeps lane order identical for R1 and R2
    cat \$(ls r1/* | sort) > ${meta.id}_R1.fastq.gz
    cat \$(ls r2/* | sort) > ${meta.id}_R2.fastq.gz
    seqkit stats -T -j ${task.cpus} ${meta.id}_R1.fastq.gz ${meta.id}_R2.fastq.gz > ${meta.id}.raw_stats.tsv
    n1=\$(awk 'NR==2{print \$4}' ${meta.id}.raw_stats.tsv)
    n2=\$(awk 'NR==3{print \$4}' ${meta.id}.raw_stats.tsv)
    if [ "\$n1" != "\$n2" ]; then
        echo "ERROR: R1 has \$n1 reads but R2 has \$n2 reads. Paired files must have the same number of reads in the same order." >&2
        echo "       Check that you did not mix files from different samples or lanes." >&2
        exit 1
    fi
    printf 'raw_pairs\\t%s\\nraw_files\\t%s\\nraw_mean_len\\t%s\\n' "\$n1" "$r1s" "\$(awk 'NR==2{print \$7}' ${meta.id}.raw_stats.tsv)" > ${meta.id}.00_raw.summary.tsv
    """
}

process FASTQC {
    tag "$meta.id"
    label 'process_low'
    container params.containers.fastqc
    publishDir path: { "${params.outdir}/${meta.id}/${stepdir}" }, mode: 'copy'

    input:
    tuple val(meta), path(r1), path(r2)
    val stepdir

    output:
    tuple val(meta), path("*.html"), emit: html
    path "*.zip", emit: zip

    script:
    """
    fastqc -o . -t ${task.cpus} --quiet $r1 $r2
    """
}

process FASTP {
    // Same options as the homework: adapter auto-detection, overlap correction, poly-G, sliding window.
    tag "$meta.id"
    label 'process_medium'
    container params.containers.fastp
    publishDir path: { "${params.outdir}/${meta.id}/02_trimmed" }, mode: 'copy'

    input:
    tuple val(meta), path(r1), path(r2)

    output:
    tuple val(meta), path("${meta.id}_R1.trim.fastq.gz"), path("${meta.id}_R2.trim.fastq.gz"), emit: reads
    path "${meta.id}.fastp.json", emit: json
    path "${meta.id}.fastp.html"

    script:
    """
    set -euo pipefail
    fastp -i $r1 -I $r2 \\
          -o ${meta.id}_R1.trim.fastq.gz -O ${meta.id}_R2.trim.fastq.gz \\
          --detect_adapter_for_pe -c --trim_poly_g \\
          --cut_right --cut_window_size 4 --cut_mean_quality ${params.cut_mean_quality} \\
          -l ${params.min_len} -w ${task.cpus} \\
          -h ${meta.id}.fastp.html -j ${meta.id}.fastp.json
    """
}

process HOST_REMOVE {
    // Optional: drop read pairs where either mate maps to the host genome.
    tag "$meta.id"
    label 'process_medium'
    container params.containers.bwasam
    publishDir path: { "${params.outdir}/${meta.id}/02_trimmed/host_removed" }, mode: 'copy', pattern: '*.{txt,tsv}'

    input:
    tuple val(meta), path(r1), path(r2)
    path host

    output:
    tuple val(meta), path("${meta.id}_R1.nohost.fastq.gz"), path("${meta.id}_R2.nohost.fastq.gz"), emit: reads
    path "${meta.id}.host_flagstat.txt"
    path "${meta.id}.02_host.summary.tsv", emit: summary

    script:
    """
    set -euo pipefail
    bwa index -p host $host 2> bwa_index.log
    bwa mem -t ${task.cpus} host $r1 $r2 2> bwa_mem.log \\
      | samtools view -b -o host.bam -
    samtools flagstat host.bam > ${meta.id}.host_flagstat.txt
    # -f 12: both mates unmapped -> these are the non-host pairs we keep
    samtools fastq -f 12 -1 ${meta.id}_R1.nohost.fastq.gz -2 ${meta.id}_R2.nohost.fastq.gz -0 /dev/null -s /dev/null -n host.bam
    total=\$(awk 'NR==1{print \$1}' ${meta.id}.host_flagstat.txt)
    kept=\$(( \$(zcat ${meta.id}_R1.nohost.fastq.gz | wc -l) / 4 ))
    printf 'host_pairs_removed\\t%s\\nhost_pct_removed\\t%.1f\\n' "\$(( total/2 - kept ))" "\$(awk -v t=\$total -v k=\$kept 'BEGIN{print (t/2-k)*100/(t/2)}')" > ${meta.id}.02_host.summary.tsv
    """
}

process CAP_READS {
    // Never feed the assembler more than --cap_pairs pairs; also make a small set for the quick genome-size assembly.
    tag "$meta.id"
    label 'process_low'
    container params.containers.seqtk
    publishDir path: { "${params.outdir}/${meta.id}/04_subsampling" }, mode: 'copy', pattern: '*.tsv'

    input:
    tuple val(meta), path(r1), path(r2)

    output:
    tuple val(meta), path("${meta.id}.cap_R1.fastq.gz"), path("${meta.id}.cap_R2.fastq.gz"), path("${meta.id}.cap_info.tsv"), emit: reads
    tuple val(meta), path("${meta.id}.quick_R1.fastq.gz"), path("${meta.id}.quick_R2.fastq.gz"), emit: quick
    path "${meta.id}.04_cap.summary.tsv", emit: summary

    script:
    """
    set -euo pipefail
    n=\$(zcat $r1 | awk 'END{print int(NR/4)}')
    # (zcat || true): head closes the pipe early, which would otherwise trip pipefail with exit 141
    len=\$( (zcat $r1 || true) | head -n 40000 | awk 'NR%4==2{s+=length(\$0);c++} END{printf "%d", (c?s/c:0)}')
    if [ "\$n" -gt ${params.cap_pairs} ]; then
        seqtk sample -s ${params.seed} $r1 ${params.cap_pairs} | gzip -1 > ${meta.id}.cap_R1.fastq.gz
        seqtk sample -s ${params.seed} $r2 ${params.cap_pairs} | gzip -1 > ${meta.id}.cap_R2.fastq.gz
        out=${params.cap_pairs}; capped=yes
    else
        ln -s \$(readlink -f $r1) ${meta.id}.cap_R1.fastq.gz
        ln -s \$(readlink -f $r2) ${meta.id}.cap_R2.fastq.gz
        out=\$n; capped=no
    fi
    if [ "\$out" -gt ${params.quick_pairs} ]; then
        seqtk sample -s ${params.seed} ${meta.id}.cap_R1.fastq.gz ${params.quick_pairs} | gzip -1 > ${meta.id}.quick_R1.fastq.gz
        seqtk sample -s ${params.seed} ${meta.id}.cap_R2.fastq.gz ${params.quick_pairs} | gzip -1 > ${meta.id}.quick_R2.fastq.gz
    else
        ln -s \$(readlink -f ${meta.id}.cap_R1.fastq.gz) ${meta.id}.quick_R1.fastq.gz
        ln -s \$(readlink -f ${meta.id}.cap_R2.fastq.gz) ${meta.id}.quick_R2.fastq.gz
    fi
    printf 'pairs_in\\t%s\\npairs_out\\t%s\\nmean_len\\t%s\\ncapped\\t%s\\n' "\$n" "\$out" "\$len" "\$capped" > ${meta.id}.cap_info.tsv
    printf 'pairs_after_cap\\t%s\\nread_cap_applied\\t%s\\n' "\$out" "\$capped" > ${meta.id}.04_cap.summary.tsv
    """
}

process QUICK_ASSEMBLY {
    // Fast SPAdes run on a small read set, only to estimate the genome size (longest contig).
    tag "$meta.id"
    label 'process_high'
    container params.containers.spades
    publishDir path: { "${params.outdir}/${meta.id}/04_subsampling/quick_assembly" }, mode: 'copy', pattern: '{genome_size_estimate.txt,quick_contigs.fasta,spades.log}'

    input:
    tuple val(meta), path(r1), path(r2)

    output:
    tuple val(meta), path("genome_size_estimate.txt"), emit: size
    path "quick_contigs.fasta"
    path "spades.log"
    path "${meta.id}.04_genome_size.summary.tsv", emit: summary

    script:
    """
    set -euo pipefail
    spades.py --isolate -1 $r1 -2 $r2 -o quick -t ${task.cpus} -m ${task.memory.toGiga()} > spades.stdout 2>&1 \\
      || { echo "ERROR: quick SPAdes assembly failed; see spades.log. You can skip this step by giving --genome_size <bp>." >&2; cp quick/spades.log . 2>/dev/null || cp spades.stdout spades.log; exit 1; }
    cp quick/spades.log spades.log
    cp quick/contigs.fasta quick_contigs.fasta
    # SPAdes names contigs NODE_1_length_39937_cov_...; contigs are sorted by length, so the first is the longest
    longest=\$(grep -m1 '>' quick_contigs.fasta | sed -E 's/.*_length_([0-9]+)_.*/\\1/')
    if [ -z "\$longest" ]; then echo "ERROR: quick assembly produced no contigs" >&2; exit 1; fi
    sum10k=\$(grep '>' quick_contigs.fasta | sed -E 's/.*_length_([0-9]+)_.*/\\1/' | awk '\$1>=10000{s+=\$1} END{print s+0}')
    if [ "\$longest" -lt 10000 ]; then
        echo "WARNING: longest quick contig is only \$longest bp; genome size estimate is unreliable. Consider --genome_size." >&2
    fi
    echo "\$longest" > genome_size_estimate.txt
    printf 'genome_size_bp\\t%s\\ngenome_size_source\\tquick_assembly_longest_contig\\nquick_sum_contigs_ge10kb\\t%s\\n' "\$longest" "\$sum10k" > ${meta.id}.04_genome_size.summary.tsv
    """
}

process SUBSAMPLE {
    // seqtk to --coverage x (main set) and --fallback_coverage x (used only if the main assembly fails).
    tag "$meta.id"
    label 'process_low'
    container params.containers.seqtk
    publishDir path: { "${params.outdir}/${meta.id}/04_subsampling" }, mode: 'copy', pattern: '*.{tsv,fastq.gz}'

    input:
    tuple val(meta), path(r1), path(r2), path(cap_info), val(genome_size)

    output:
    tuple val(meta), path("${meta.id}.main_R1.fastq.gz"), path("${meta.id}.main_R2.fastq.gz"),
                     path("${meta.id}.fallback_R1.fastq.gz"), path("${meta.id}.fallback_R2.fastq.gz"),
                     path("${meta.id}.subsample_info.tsv"), emit: reads
    path "${meta.id}.04_subsample.summary.tsv", emit: summary

    script:
    """
    set -euo pipefail
    n=\$(awk '\$1=="pairs_out"{print \$2}' $cap_info)
    len=\$(awk '\$1=="mean_len"{print \$2}' $cap_info)
    G=${genome_size}
    cur=\$(awk -v n=\$n -v l=\$len -v g=\$G 'BEGIN{printf "%d", n*2*l/g}')
    sub() {   # sub <cov> <tag>
        keep=\$(awk -v c=\$1 -v l=\$len -v g=\$G 'BEGIN{printf "%d", c*g/(2*l)}')
        if [ "\$keep" -ge "\$n" ]; then
            ln -sf \$(readlink -f $r1) ${meta.id}.\$2_R1.fastq.gz; ln -sf \$(readlink -f $r2) ${meta.id}.\$2_R2.fastq.gz
            echo "\$n"
        else
            seqtk sample -s ${params.seed} $r1 \$keep | gzip -1 > ${meta.id}.\$2_R1.fastq.gz
            seqtk sample -s ${params.seed} $r2 \$keep | gzip -1 > ${meta.id}.\$2_R2.fastq.gz
            echo "\$keep"
        fi
    }
    main=\$(sub ${params.coverage} main)
    fb=\$(sub ${params.fallback_coverage} fallback)
    same=no; [ "\$main" = "\$fb" ] && same=yes
    printf 'genome_size\\t%s\\npairs_available\\t%s\\nmean_len\\t%s\\ncoverage_available\\t%s\\ntarget_coverage\\t%s\\nmain_pairs\\t%s\\nfallback_coverage\\t%s\\nfallback_pairs\\t%s\\nfallback_same_as_main\\t%s\\n' \\
        "\$G" "\$n" "\$len" "\$cur" ${params.coverage} "\$main" ${params.fallback_coverage} "\$fb" "\$same" > ${meta.id}.subsample_info.tsv
    printf 'coverage_available\\t%s\\nmain_pairs\\t%s\\nmain_coverage_target\\t%s\\n' "\$cur" "\$main" ${params.coverage} > ${meta.id}.04_subsample.summary.tsv
    """
}
