#!/usr/bin/env bash
# termini_peaks.sh <sorted.bam> <ref.fasta.fai> [min_dist_from_contig_end] [min_contig_length]
# Reports the strongest pile-up of read 5' ends on each strand (forward read starts, reverse read ends),
# ignoring positions near contig ends. Physical genome termini (cos, DTR, pac) in a library whose adapters
# were ligated to fragment ends give a peak holding a large fraction of the local depth; tagmentation
# libraries or circularly permuted genomes do not. Output: TSV strand contig pos reads depth fraction
set -euo pipefail
BAM=$1; FAI=$2; EDGE=${3:-10}; MINLEN=${4:-10000}
# only contigs long enough to be a phage genome; if none, the longest contig
awk -v m="$MINLEN" '$2>=m' "$FAI" > .long.fai; [ -s .long.fai ] || sort -k2,2nr "$FAI" | head -1 > .long.fai
samtools view -F 0x914 "$BAM" | awk -v OFS='\t' '{c[$3"\t"$4]++} END{for(k in c) print "forward_start", k, c[k]}' > .fwd.tsv
samtools view -f 0x10 -F 0x904 "$BAM" | awk -v OFS='\t' '{
    p=$4; cig=$6
    while (match(cig,/^[0-9]+[MIDNSHP=X]/)) { n=substr(cig,RSTART,RLENGTH-1); op=substr(cig,RLENGTH,1); if (op ~ /[MDN=X]/) p+=n; cig=substr(cig,RLENGTH+1) }
    c[$3"\t"(p-1)]++ } END{for(k in c) print "reverse_end", k, c[k]}' > .rev.tsv
cat .fwd.tsv .rev.tsv | awk -v OFS='\t' -v E="$EDGE" 'NR==FNR{len[$1]=$2; next} ($2 in len) && ($3>E && $3<=len[$2]-E)' .long.fai - | sort -k1,1 -k4,4nr \
  | awk '!seen[$1]++' | while IFS=$'\t' read -r strand contig pos reads; do
      depth=$(samtools depth -a -r "$contig:$pos-$pos" "$BAM" | awk '{print $3+0}')
      printf '%s\t%s\t%s\t%s\t%s\t%.2f\n' "$strand" "$contig" "$pos" "$reads" "${depth:-0}" "$(awk -v r=$reads -v d=${depth:-0} 'BEGIN{print (d>0?r/d:0)}')"
    done
rm -f .fwd.tsv .rev.tsv .long.fai
