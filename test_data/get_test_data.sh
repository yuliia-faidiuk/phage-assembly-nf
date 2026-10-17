#!/usr/bin/env bash
# Downloads two small public phage Illumina runs (paired-end) used by the pipeline's test profile.
#   SRR37741590  Escherichia phage T7       (~146 k pairs, NextSeq 2000, genome 39 937 bp; evolution-study library, no termini signal)
#   DRR406538    Escherichia phage phiEc_1  (~157 k pairs, MiSeq; DTR-short phage, clear termini signal)
# Files are renamed to the Illumina facility pattern so the sample detection is exercised.
set -euo pipefail
cd "$(dirname "$0")"
dl() { [[ -s "$2" ]] && { echo "have $2"; return; }; echo "downloading $2"; curl -sSL --retry 5 -o "$2.part" "https://$1" && mv "$2.part" "$2"; }
dl ftp.sra.ebi.ac.uk/vol1/fastq/SRR377/090/SRR37741590/SRR37741590_1.fastq.gz T7_S1_L001_R1_001.fastq.gz
dl ftp.sra.ebi.ac.uk/vol1/fastq/SRR377/090/SRR37741590/SRR37741590_2.fastq.gz T7_S1_L001_R2_001.fastq.gz
dl ftp.sra.ebi.ac.uk/vol1/fastq/DRR406/DRR406538/DRR406538_1.fastq.gz        phiEc1_S2_L001_R1_001.fastq.gz
dl ftp.sra.ebi.ac.uk/vol1/fastq/DRR406/DRR406538/DRR406538_2.fastq.gz        phiEc1_S2_L001_R2_001.fastq.gz
md5sum -c <<'SUM'
4979d68debf3773e9ff15f5146e6804a  T7_S1_L001_R1_001.fastq.gz
1e09e042c976c13051c0736c314687d3  T7_S1_L001_R2_001.fastq.gz
9e56fc14f1798c23329dc5b45e072639  phiEc1_S2_L001_R1_001.fastq.gz
eafc749fad44798e961001a669c7c490  phiEc1_S2_L001_R2_001.fastq.gz
SUM
echo "test data ready in $(pwd)"
