// Steps 7-10 of the homework: map reads back, mixed positions, Pilon, PhageTerm, rotation.

process MAP_READS {
    tag "$meta.id"
    label 'process_medium'
    container params.containers.bwasam
    publishDir path: { "${params.outdir}/${meta.id}/07_mapping" }, mode: 'copy'

    input:
    tuple val(meta), path(r1), path(r2), val(asm_name), path(asm)

    output:
    tuple val(meta), path("${meta.id}.fasta"), path("${meta.id}.fasta.fai"), path("${meta.id}.bam"), path("${meta.id}.bam.bai"), emit: bam
    path "${meta.id}.flagstat.txt", emit: flagstat
    path "${meta.id}.samtools_stats.txt", emit: stats
    path "${meta.id}.coverage.txt"
    path "${meta.id}.depth.tsv.gz"
    path "${meta.id}.igv_session.xml"
    path "${meta.id}.bwa.log"
    tuple val(meta), path("${meta.id}.termini_peaks.tsv"), emit: termini
    path "${meta.id}.07_mapping.summary.tsv", emit: summary

    script:
    """
    set -euo pipefail
    cp $asm/assembly.fasta ${meta.id}.fasta
    bwa index ${meta.id}.fasta 2> ${meta.id}.bwa.log
    samtools faidx ${meta.id}.fasta
    bwa mem -t ${task.cpus} ${meta.id}.fasta $r1 $r2 2>> ${meta.id}.bwa.log \\
      | samtools sort -@ ${task.cpus} -o ${meta.id}.bam -
    samtools index ${meta.id}.bam
    samtools flagstat ${meta.id}.bam > ${meta.id}.flagstat.txt
    samtools stats ${meta.id}.bam > ${meta.id}.samtools_stats.txt
    samtools coverage ${meta.id}.bam > ${meta.id}.coverage.txt
    samtools depth -a ${meta.id}.bam | gzip > ${meta.id}.depth.tsv.gz
    cat > ${meta.id}.igv_session.xml <<XML
<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<Session genome="${meta.id}.fasta" locus="All" version="8">
    <Resources>
        <Resource path="${meta.id}.bam"/>
    </Resources>
</Session>
XML
    mapped=\$(grep -m1 'mapped (' ${meta.id}.flagstat.txt | sed -E 's/.*\\(([0-9.]+)%.*/\\1/')
    proper=\$(grep -m1 'properly paired' ${meta.id}.flagstat.txt | sed -E 's/.*\\(([0-9.]+)%.*/\\1/')
    depth=\$(awk 'NR>1{d+=\$7*\$3; l+=\$3} END{printf "%.0f", (l?d/l:0)}' ${meta.id}.coverage.txt)
    covpct=\$(awk 'NR>1{c+=\$5; l+=\$3} END{printf "%.2f", (l?c*100/l:0)}' ${meta.id}.coverage.txt)
    # strongest pile-up of read 5' ends per strand: physical termini in a ligation library give a peak with a large share of the depth
    printf 'strand\\tcontig\\tpos\\treads\\tdepth\\tfraction_of_depth\\n' > ${meta.id}.termini_peaks.tsv
    termini_peaks.sh ${meta.id}.bam ${meta.id}.fasta.fai 10 ${params.min_phage_contig} >> ${meta.id}.termini_peaks.tsv
    fwd=\$(awk -F'\\t' '\$1=="forward_start"{printf "%s reads at %s:%s (%.0f%% of depth)", \$4, \$2, \$3, \$6*100}' ${meta.id}.termini_peaks.tsv)
    rev=\$(awk -F'\\t' '\$1=="reverse_end"{printf "%s reads at %s:%s (%.0f%% of depth)", \$4, \$2, \$3, \$6*100}' ${meta.id}.termini_peaks.tsv)
    maxfrac=\$(awk -F'\\t' 'NR>1 && \$6>m{m=\$6} END{print m+0}' ${meta.id}.termini_peaks.tsv)
    signal=\$(awk -v m=\$maxfrac 'BEGIN{ if (m>=0.25) print "present (" int(m*100) "% of depth): fixed ends are visible, PhageTerm can work"; else print "weak (max " int(m*100) "% of depth): no fragment-end signal, PhageTerm class is unreliable (tagmentation library or no fixed ends)" }')
    printf 'mapped_to\\t%s\\npct_mapped\\t%s\\npct_properly_paired\\t%s\\nmean_depth\\t%s\\npct_genome_covered\\t%s\\ntermini_peak_forward\\t%s\\ntermini_peak_reverse\\t%s\\ntermini_signal\\t%s\\n' \\
        "$asm_name" "\$mapped" "\$proper" "\$depth" "\$covpct" "\${fwd:-none}" "\${rev:-none}" "\$signal" > ${meta.id}.07_mapping.summary.tsv
    """
}

process VARIANTS {
    // Positions where the reads disagree with the contig, with the fraction of reads carrying the alternative base.
    tag "$meta.id"
    label 'process_low'
    container params.containers.bcftools
    publishDir path: { "${params.outdir}/${meta.id}/08_variants" }, mode: 'copy'

    input:
    tuple val(meta), path(ref), path(fai), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.id}.calls.vcf.gz"), path("${meta.id}.calls.vcf.gz.csi"), emit: vcf
    path "${meta.id}.mixed_positions.tsv"
    path "${meta.id}.08_variants.summary.tsv", emit: summary

    script:
    """
    set -euo pipefail
    bcftools mpileup -f $ref -a AD,DP -d 10000 $bam 2> mpileup.log \\
      | bcftools call --ploidy 1 -mv -Oz -o ${meta.id}.calls.vcf.gz
    bcftools index -f ${meta.id}.calls.vcf.gz
    printf 'contig\\tpos\\tref\\talt\\tqual\\tgenotype\\tallele_depths\\tdepth\\talt_fraction\\n' > ${meta.id}.mixed_positions.tsv
    bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%QUAL\\t[%GT\\t%AD\\t%DP]\\n' ${meta.id}.calls.vcf.gz \\
      | awk -F'\\t' 'BEGIN{OFS="\\t"}{ split(\$7,a,","); af=(\$8>0)?1-a[1]/\$8:0; printf "%s\\t%s\\t%s\\t%s\\t%.0f\\t%s\\t%s\\t%s\\t%.2f\\n",\$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8,af }' >> ${meta.id}.mixed_positions.tsv
    n=\$(( \$(wc -l < ${meta.id}.mixed_positions.tsv) - 1 ))
    nmix=\$(awk -F'\\t' 'NR>1 && \$9>=0.2 && \$9<=0.8' ${meta.id}.mixed_positions.tsv | wc -l)
    printf 'variant_positions\\t%s\\nmixed_positions_20_80pct\\t%s\\n' "\$n" "\$nmix" > ${meta.id}.08_variants.summary.tsv
    """
}

process PILON {
    tag "$meta.id"
    label 'process_medium'
    container params.containers.pilon
    publishDir path: { "${params.outdir}/${meta.id}/08_variants/pilon" }, mode: 'copy'

    input:
    tuple val(meta), path(ref), path(fai), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.id}.pilon.fasta"), emit: fasta
    path "${meta.id}.pilon.changes"
    path "pilon.log"
    path "${meta.id}.08_pilon.summary.tsv", emit: summary

    script:
    """
    set -euo pipefail
    export _JAVA_OPTIONS="-Xmx${task.memory.toGiga() - 1}g"
    pilon --genome $ref --frags $bam --output ${meta.id}.pilon --outdir . --changes --fix all > pilon.log 2>&1
    # Pilon appends "_pilon" to every contig name; remove it so names match the mapping reference
    sed -i 's/_pilon\$//' ${meta.id}.pilon.fasta
    n=\$(wc -l < ${meta.id}.pilon.changes)
    printf 'pilon_changes\\t%s\\n' "\$n" > ${meta.id}.08_pilon.summary.tsv
    """
}

process PHAGETERM {
    // Optional: PhageTermVirome on the trimmed reads vs the polished assembly. Never stops the pipeline.
    tag "$meta.id"
    label 'process_medium'
    label 'optional_step'
    container params.phageterm_image
    publishDir path: { "${params.outdir}/${meta.id}/09_phageterm" }, mode: 'copy', saveAs: { it == 'for_rotation' ? null : it }

    input:
    tuple val(meta), path(r1), path(r2), path(ref)

    output:
    tuple val(meta), path("for_rotation"), optional: true, emit: seq     // class.txt + sequence.fasta, only when PhageTerm re-opened the genome
    path "${meta.id}*.pdf", optional: true
    path "${meta.id}*.csv", optional: true
    path "${meta.id}*sequence*.fasta", optional: true
    path "${meta.id}.phageterm_class.txt", optional: true
    path "phageterm.log"
    path "${meta.id}.09_phageterm.summary.tsv", emit: summary

    script:
    """
    set -uo pipefail
    export MPLCONFIGDIR=\$PWD/.mpl HOME=\$PWD; mkdir -p .mpl
    zcat $r1 > R1.fastq; zcat $r2 > R2.fastq
    # PhageTerm gets only contigs long enough to be a phage genome (it crashes on multi-contig input with low-coverage fragments)
    awk -v m=${params.min_phage_contig} 'BEGIN{RS=">"; ORS=""} NR>1{h=\$0; sub(/\\n.*/,"",h); s=\$0; sub(/^[^\\n]*\\n/,"",s); gsub(/\\n/,"",s); if (length(s)>=m) print ">"h"\\n"s"\\n"}' $ref > ref_long.fasta
    if ! [ -s ref_long.fasta ]; then
        awk 'BEGIN{RS=">"; ORS=""} NR>1{h=\$0; sub(/\\n.*/,"",h); s=\$0; sub(/^[^\\n]*\\n/,"",s); gsub(/\\n/,"",s); if (length(s)>L){L=length(s); H=h; S=s}} END{print ">"H"\\n"S"\\n"}' $ref > ref_long.fasta
    fi
    echo "contigs given to PhageTerm (>= ${params.min_phage_contig} bp, or the longest): \$(grep -c '>' ref_long.fasta)" > phageterm.log
    phageterm -f R1.fastq -p R2.fastq -r ref_long.fasta --report_title ${meta.id} -c ${task.cpus} --nrt >> phageterm.log 2>&1
    rc=\$?
    rm -f R1.fastq R2.fastq
    # --nrt makes PhageTerm write its packaging class (e.g. "DTR (short)", "COS (5')", "Headful (pac)", "UNKNOWN") to nrt.txt
    cls=\$(cat nrt.txt 2>/dev/null | head -1 | xargs || true)
    case "\$cls" in ""|"-") cls="undetermined (PhageTerm found no clear termini)";; esac
    if [ \$rc -ne 0 ]; then
        printf 'phageterm_class\\tFAILED (see 09_phageterm/phageterm.log)\\n' > ${meta.id}.09_phageterm.summary.tsv
        exit \$rc
    fi
    echo "\${cls:-see report}" > ${meta.id}.phageterm_class.txt
    seq=\$(ls ${meta.id}*sequence*.fasta 2>/dev/null | head -1 || true)
    if [ -n "\$seq" ] && [ -s "\$seq" ]; then
        mkdir -p for_rotation; cp "\$seq" for_rotation/sequence.fasta; cp ${meta.id}.phageterm_class.txt for_rotation/class.txt
    fi
    printf 'phageterm_class\\t%s\\nphageterm_reopened_sequence\\t%s\\n' "\${cls:-see report}" "\$([ -n "\$seq" ] && echo yes || echo no)" > ${meta.id}.09_phageterm.summary.tsv
    """
}

process ROTATE {
    // Final genome: PhageTerm's re-opened sequence when fixed termini were found, otherwise dnaapler (start upstream of terminase large subunit).
    tag "$meta.id"
    label 'process_medium'
    container params.containers.dnaapler
    publishDir path: { "${params.outdir}/${meta.id}/10_final" }, mode: 'copy'

    input:
    tuple val(meta), path(polished), path(termini), path(pt, stageAs: 'phageterm')     // pt = PhageTerm's for_rotation/ folder, or [] if none

    output:
    path "${meta.id}.before_rotation.fasta"
    path "${meta.id}.after_rotation.fasta"
    path "rotation_method.txt"
    path "dnaapler", optional: true
    path "${meta.id}.10_final.summary.tsv", emit: summary

    script:
    """
    set -euo pipefail
    export HOME=\$PWD
    cp $polished ${meta.id}.before_rotation.fasta
    # split: contigs >= min_phage_contig are rotated, shorter ones are appended unchanged
    awk -v m=${params.min_phage_contig} 'BEGIN{RS=">"; ORS=""} NR>1{h=\$0; sub(/\\n.*/,"",h); s=\$0; sub(/^[^\\n]*\\n/,"",s); gsub(/\\n/,"",s); f=(length(s)>=m)?"long.fasta":"short.fasta"; print ">"h"\\n"s"\\n" > f}' $polished
    touch long.fasta short.fasta
    if ! [ -s long.fasta ]; then   # nothing long enough: rotate the longest contig
        awk 'BEGIN{RS=">"; ORS=""} NR>1{h=\$0; sub(/\\n.*/,"",h); s=\$0; sub(/^[^\\n]*\\n/,"",s); gsub(/\\n/,"",s); if (length(s)>L){L=length(s); H=h; S=s}} END{print ">"H"\\n"S"\\n"}' $polished > long.fasta
        awk 'BEGIN{RS=">"; ORS=""} NR>1{h=\$0; sub(/\\n.*/,"",h); s=\$0; sub(/^[^\\n]*\\n/,"",s); gsub(/\\n/,"",s); if (length(s)>L){L=length(s)}; a[NR]=">"h"\\n"s"\\n"; l[NR]=length(s)} END{for(i in a) if(l[i]<L) print a[i]}' $polished > short.fasta
    fi
    nshort=\$(grep -c '>' short.fasta || true)
    cls=\$(cat phageterm/class.txt 2>/dev/null | head -1 || true)
    # PhageTerm re-opens the genome at its termini for cos / DTR phages and at the pac site for headful phages;
    # for "UNKNOWN" or Mu-like results its re-opened sequence is not meaningful, so dnaapler decides instead.
    # ... and only when the reads actually show fragment ends (termini_signal present), otherwise PhageTerm's call is a guess
    maxfrac=\$(awk -F'\\t' 'NR>1 && \$6>m{m=\$6} END{print m+0}' $termini)
    signal_ok=\$(awk -v m=\$maxfrac 'BEGIN{print (m>=0.25)?"yes":"no"}')
    if [ -s phageterm/sequence.fasta ] && [ "\$signal_ok" = yes ] && echo "\$cls" | grep -qiE 'COS|DTR|pac|headful'; then
        cat phageterm/sequence.fasta short.fasta > ${meta.id}.after_rotation.fasta
        echo "phageterm: genome re-opened by PhageTerm, class \$cls; \$nshort short contig(s) (< ${params.min_phage_contig} bp) appended unchanged" > rotation_method.txt
    elif dnaapler all -i long.fasta -o dnaapler -p ${meta.id} -t ${task.cpus} -f --autocomplete largest > dnaapler.stdout 2>&1; then
        cat dnaapler/${meta.id}_reoriented.fasta short.fasta > ${meta.id}.after_rotation.fasta
        gene=\$(awk -F'\\t' 'NR==2{print \$2}' dnaapler/${meta.id}_all_reorientation_summary.tsv 2>/dev/null || echo unknown)
        case "\$gene" in
            terL) how="start placed at the terminase large subunit (terL) found by dnaapler" ;;
            autocomplete*) how="no terminase found by dnaapler, start placed at the largest ORF" ;;
            *) how="start placed by dnaapler at gene \$gene" ;;
        esac
        why="PhageTerm class: \${cls:-not run}"; [ -n "\$cls" ] && [ "\$signal_ok" = no ] && why="PhageTerm class '\$cls' ignored because the reads show no fragment-end signal (termini_signal weak)"
        echo "dnaapler: \$how; \$why; \$nshort short contig(s) appended unchanged" > rotation_method.txt
        mv dnaapler.stdout dnaapler/
    else
        cp $polished ${meta.id}.after_rotation.fasta
        mkdir -p dnaapler; mv dnaapler.stdout dnaapler/
        echo "none: dnaapler failed (see 10_final/dnaapler/dnaapler.stdout), sequence left as assembled; PhageTerm class: \${cls:-not run}" > rotation_method.txt
    fi
    n=\$(grep -c '>' ${meta.id}.after_rotation.fasta)
    len=\$(grep -v '>' ${meta.id}.after_rotation.fasta | tr -d '\\n' | wc -c)
    printf 'final_contigs\\t%s\\nfinal_total_bp\\t%s\\nrotation\\t%s\\n' "\$n" "\$len" "\$(cut -d: -f1 rotation_method.txt)" > ${meta.id}.10_final.summary.tsv
    """
}

process IGV_REPORT {
    // Self-contained HTML with the read pileups (igv.js) at every mixed position, the termini peaks,
    // contig ends and a few windows; opens in any browser without IGV, internet or a server.
    tag "$meta.id"
    label 'process_low'
    container params.containers.igvreports
    publishDir path: { "${params.outdir}/${meta.id}/07_mapping" }, mode: 'copy'

    input:
    tuple val(meta), path(ref), path(fai), path(bam), path(bai), path(termini), path(vcf), path(csi)

    output:
    path "${meta.id}.alignment_report.html"
    path "HOW_TO_VIEW_ALIGNMENTS.txt"

    script:
    """
    set -euo pipefail
    export HOME=\$PWD
    make_igv_sites.py --fai $fai --termini $termini --vcf $vcf --min-contig ${params.min_phage_contig} --out sites.bed
    create_report sites.bed --fasta $ref --tracks $bam $vcf --flanking 400 --standalone \\
        --title "${meta.id}: reads mapped to the assembly" --output ${meta.id}.alignment_report.html > create_report.log 2>&1
    cat > HOW_TO_VIEW_ALIGNMENTS.txt <<TXT
How to look at the reads mapped to the assembly (homework steps 7 and 8)

A. Quick look, no installation
   Open ${meta.id}.alignment_report.html in any browser. The table lists the positions worth looking at
   (every position where reads disagree with the contig, the termini peaks, the contig ends and three
   windows); click a row to see the read pileup around it. Each row shows 400 bp on each side.

B. Whole genome, interactive: IGV desktop (https://igv.org/doc/desktop/#DownloadPage)
   On WSL, first copy this whole 07_mapping folder to a Windows folder (e.g. explorer.exe . ), then in IGV:
   1. Genomes > Load Genome from File...   -> ${meta.id}.fasta
   2. File > Load from File...             -> ${meta.id}.bam      (the .bam.bai index must sit next to it)
   3. Type a position in the box at the top (e.g. 1:12000-12500) or zoom with + / -.
   Or in one go: File > Open Session...    -> ${meta.id}.igv_session.xml   (works when the files stayed together)
   Right-click the read track > "Collapsed" to see depth, > "Color alignments by > read strand" for termini,
   > "Sort alignments by > base" at a mixed position.

Files: ${meta.id}.fasta (contigs the reads were mapped to), ${meta.id}.bam/.bai (the alignments),
       ${meta.id}.coverage.txt (depth per contig), ${meta.id}.termini_peaks.tsv (read-start pile-ups),
       ../08_variants/${meta.id}.mixed_positions.tsv (the positions listed in the report).
TXT
    """
}
