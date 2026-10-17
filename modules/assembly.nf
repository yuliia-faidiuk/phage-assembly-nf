// Steps 5-6 of the homework: assemble (Unicycler in several modes + SPAdes) and draw the graph.

process UNICYCLER {
    // Tries the main (200x) read set first; if Unicycler fails or gives no contigs, retries with the fallback (100x) set.
    tag "${meta.id}:${mode}"
    label 'process_high'
    container params.containers.unicycler
    publishDir path: { "${params.outdir}/${meta.id}/05_assembly" }, mode: 'copy'

    input:
    tuple val(meta), path(m1), path(m2), path(f1), path(f2), path(info), val(mode)

    output:
    tuple val(meta), val("unicycler_${mode}"), path("unicycler_${mode}"), emit: assembly
    path "${meta.id}.05_unicycler_${mode}.summary.tsv", emit: summary

    script:
    """
    set -uo pipefail
    OUT=unicycler_${mode}
    same=\$(awk '\$1=="fallback_same_as_main"{print \$2}' $info)
    run_uni() {   # run_uni <R1> <R2>
        rm -rf "\${OUT:?}"
        unicycler -1 \$1 -2 \$2 -o "\$OUT" --mode $mode --keep 0 -t ${task.cpus} > unicycler.stdout 2>&1
        rc=\$?
        [ \$rc -eq 0 ] && [ -s "\$OUT/assembly.fasta" ] && grep -q '>' "\$OUT/assembly.fasta"
    }
    cov_used=${params.coverage}
    if ! run_uni $m1 $m2; then
        if [ "\$same" = "yes" ]; then
            echo "ERROR: Unicycler ($mode) failed and no smaller read set is available. See \$OUT/unicycler.log" >&2
            mkdir -p "\$OUT"; cp unicycler.stdout "\$OUT/unicycler.stdout"; exit 1
        fi
        echo "Unicycler ($mode) failed at ${params.coverage}x, retrying with the ${params.fallback_coverage}x read set" | tee -a unicycler.stdout
        mv unicycler.stdout unicycler.first_attempt.stdout
        if ! run_uni $f1 $f2; then
            echo "ERROR: Unicycler ($mode) failed at both ${params.coverage}x and ${params.fallback_coverage}x. See \$OUT/unicycler.log" >&2
            mkdir -p "\$OUT"; cp unicycler.stdout unicycler.first_attempt.stdout "\$OUT/"; exit 1
        fi
        cov_used=${params.fallback_coverage}
        cp unicycler.first_attempt.stdout "\$OUT/"
    fi
    echo "\$cov_used" > "\$OUT/coverage_used.txt"
    n=\$(grep -c '>' "\$OUT/assembly.fasta")
    longest=\$(grep '>' "\$OUT/assembly.fasta" | head -1 | sed -E 's/.*length=([0-9]+).*/\\1/')
    circ=\$(grep '>' "\$OUT/assembly.fasta" | head -1 | grep -q 'circular=true' && echo yes || echo no)
    depth=\$(grep '>' "\$OUT/assembly.fasta" | head -1 | sed -E 's/.*depth=([0-9.]+)x.*/\\1/')
    n10k=\$(grep '>' "\$OUT/assembly.fasta" | sed -E 's/.*length=([0-9]+).*/\\1/' | awk '\$1>=10000' | wc -l)
    printf 'unicycler_${mode}:n_contigs\\t%s\\nunicycler_${mode}:longest_bp\\t%s\\nunicycler_${mode}:longest_circular\\t%s\\nunicycler_${mode}:longest_depth\\t%s\\nunicycler_${mode}:contigs_ge10kb\\t%s\\nunicycler_${mode}:coverage_used\\t%s\\n' \\
        "\$n" "\$longest" "\$circ" "\$depth" "\$n10k" "\$cov_used" > ${meta.id}.05_unicycler_${mode}.summary.tsv
    """
}

process SPADES {
    // Plain SPAdes --isolate for comparison. Circular contigs carry a k-mer-length overlap at their ends (see homework step 5).
    tag "$meta.id"
    label 'process_high'
    container params.containers.spades
    publishDir path: { "${params.outdir}/${meta.id}/05_assembly" }, mode: 'copy'

    input:
    tuple val(meta), path(m1), path(m2), path(f1), path(f2), path(info)

    output:
    tuple val(meta), val("spades"), path("spades"), emit: assembly
    path "${meta.id}.05_spades.summary.tsv", emit: summary

    script:
    """
    set -uo pipefail
    same=\$(awk '\$1=="fallback_same_as_main"{print \$2}' $info)
    run_sp() {
        rm -rf sp
        spades.py --isolate -1 \$1 -2 \$2 -o sp -t ${task.cpus} -m ${task.memory.toGiga()} > spades.stdout 2>&1
        rc=\$?
        [ \$rc -eq 0 ] && [ -s sp/contigs.fasta ]
    }
    cov_used=${params.coverage}
    if ! run_sp $m1 $m2; then
        if [ "\$same" = "yes" ]; then echo "ERROR: SPAdes failed, see spades/spades.log" >&2; mkdir -p spades; cp sp/spades.log spades/ 2>/dev/null; exit 1; fi
        mv spades.stdout spades.first_attempt.stdout
        if ! run_sp $f1 $f2; then echo "ERROR: SPAdes failed at both coverages, see spades/spades.log" >&2; mkdir -p spades; cp sp/spades.log spades/ 2>/dev/null; exit 1; fi
        cov_used=${params.fallback_coverage}
    fi
    mkdir -p spades
    cp sp/contigs.fasta spades/assembly.fasta
    cp sp/assembly_graph_with_scaffolds.gfa spades/assembly.gfa
    cp sp/spades.log spades/spades.log
    cp sp/params.txt spades/params.txt
    [ -f spades.first_attempt.stdout ] && cp spades.first_attempt.stdout spades/
    echo "\$cov_used" > spades/coverage_used.txt
    n=\$(grep -c '>' spades/assembly.fasta)
    longest=\$(grep -m1 '>' spades/assembly.fasta | sed -E 's/.*_length_([0-9]+)_.*/\\1/')
    n10k=\$(grep '>' spades/assembly.fasta | sed -E 's/.*_length_([0-9]+)_.*/\\1/' | awk '\$1>=10000' | wc -l)
    printf 'spades:n_contigs\\t%s\\nspades:longest_bp\\t%s\\nspades:contigs_ge10kb\\t%s\\nspades:coverage_used\\t%s\\n' "\$n" "\$longest" "\$n10k" "\$cov_used" > ${meta.id}.05_spades.summary.tsv
    """
}

process BANDAGE_IMAGE {
    // Draws the assembly graph to PNG (what you would see in Bandage) and writes the graph summary table.
    tag "${meta.id}:${name}"
    label 'process_low'
    container params.containers.bandage
    publishDir path: { "${params.outdir}/${meta.id}/06_bandage" }, mode: 'copy'

    input:
    tuple val(meta), val(name), path(asm)

    output:
    path "${meta.id}_${name}.png"
    path "${meta.id}_${name}.graph_info.tsv"

    script:
    """
    export QT_QPA_PLATFORM=offscreen XDG_RUNTIME_DIR=\$PWD/.xdg; mkdir -p .xdg
    Bandage image $asm/assembly.gfa ${meta.id}_${name}.png --names --lengths --depth --colour depth --width 1600
    Bandage info  $asm/assembly.gfa --tsv > ${meta.id}_${name}.graph_info.tsv
    """
}

process CONTIG_TABLE {
    // One table per assembly: contig, length, depth, circular, GC.
    tag "${meta.id}:${name}"
    label 'process_low'
    container params.containers.seqkit
    publishDir path: { "${params.outdir}/${meta.id}/05_assembly" }, mode: 'copy'

    input:
    tuple val(meta), val(name), path(asm)

    output:
    path "${meta.id}_${name}.contigs.tsv"

    script:
    """
    printf 'contig\\tlength\\tgc_pct\\tdepth\\tcircular\\n' > ${meta.id}_${name}.contigs.tsv
    seqkit fx2tab -n -l -g $asm/assembly.fasta | awk -F'\\t' 'BEGIN{OFS="\\t"}{
        h=\$1; d="NA"; c="NA";
        if (match(h,/depth=[0-9.]+x/))  d=substr(h,RSTART+6,RLENGTH-7);
        if (match(h,/_cov_[0-9.]+/))     d=substr(h,RSTART+5,RLENGTH-5);
        if (h ~ /circular=true/) c="yes"; else if (h ~ /circular=/ || "${name}" != "spades") c="no";   # Unicycler only tags circular contigs
        split(h,a," "); print a[1], \$2, \$3, d, c }' >> ${meta.id}_${name}.contigs.tsv
    """
}
