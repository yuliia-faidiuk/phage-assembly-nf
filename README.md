# phage-assembly-nf — from raw Illumina reads to a finished phage genome, in one command

This pipeline runs every command-line step of the *From Reads to Phage Genome* homework
(Illumina path) automatically, for one or many phages at once, and leaves you the same files
the homework asks you to look at. You still do the *looking* (Bandage, IGV, judging mixed positions,
reading the PhageTerm report); the pipeline does the typing.

```
raw FASTQ ─► FastQC ─► fastp ─► FastQC ─► cap at 500 000 pairs ─► quick assembly (genome size)
         ─► subsample to 200x ─► Unicycler normal / conservative / bold + SPAdes ─► Bandage PNGs
         ─► map reads back (bwa) ─► mixed positions (bcftools) ─► Pilon ─► PhageTerm ─► rotation (dnaapler)
         ─► summary.md + MultiQC
```

## 1. One-time setup (10 minutes)

You need three programs. On Windows, do all of this inside **WSL2 Ubuntu** (not PowerShell).

| program | check it works | install if missing |
|---|---|---|
| Docker | `docker info` | Windows/Mac: install *Docker Desktop*, start it, and in Settings → Resources → WSL integration enable your Ubuntu. Linux: `sudo apt install docker.io && sudo usermod -aG docker $USER` then log out and in. |
| Java 17+ | `java -version` | `sudo apt install openjdk-17-jre-headless` |
| Nextflow | `nextflow -version` | `curl -s https://get.nextflow.io \| bash && mkdir -p ~/.local/bin && mv nextflow ~/.local/bin/` then open a new terminal |

Then get the pipeline and run the pre-flight check; it tells you exactly what is missing:

```bash
git clone https://github.com/yuliia-faidiuk/phage-assembly-nf.git
bash phage-assembly-nf/bin/doctor.sh
```

The first run downloads about 6 GB of tool containers. After that everything works offline.

## 2. Run it

Put the `.fastq.gz` files of your phage(s) in one folder (copy them to a folder under `/home`,
not under `/mnt/c`, which is slow). File names from the sequencing facility are fine:

```
~/reads/run1/
    phageA_S1_L001_R1_001.fastq.gz     phageA_S1_L001_R2_001.fastq.gz
    phageB_S2_L001_R1_001.fastq.gz     phageB_S2_L001_R2_001.fastq.gz
    phageB_S2_L002_R1_001.fastq.gz     phageB_S2_L002_R2_001.fastq.gz     <- lanes are merged automatically
```

```bash
bash phage-assembly-nf/run.sh ~/reads/run1
```

The script checks your computer, shows the samples it found and asks you to confirm,
then runs everything. A phage of 40–150 kb takes about 20–40 minutes per sample on a laptop
(four assemblies are run per sample). Results go to `results/` in the current folder.

Useful options:

```bash
bash phage-assembly-nf/run.sh ~/reads/run1 --genome-size 45000   # skip the genome-size estimate if you know it
bash phage-assembly-nf/run.sh ~/reads/run1 --outdir ~/phageA_results
bash phage-assembly-nf/run.sh ~/reads/run1 --quick               # Unicycler normal only, no SPAdes (3x faster)
bash phage-assembly-nf/run.sh ~/reads/run1 --host ecoli.fasta    # remove host reads before assembly
bash phage-assembly-nf/run.sh ~/reads/run1 --skip-phageterm
bash phage-assembly-nf/run.sh --test                              # two small public phages, ~15 min; proves the setup works
```

If something fails, fix the cause and run **the same command again**: finished steps are
reused, only the failed step and the ones after it are rerun.

## 3. What you get

`results/summary.md` is the first thing to open: one table, one column per phage, with every
number the homework asks for (reads before/after trimming, coverage, contigs, circular or not,
% mapped, mixed positions, Pilon changes, PhageTerm class, final length).

Per sample, the folders follow the homework steps:

| folder | homework step | open this |
|---|---|---|
| `00_raw_reads/` | get the data | `*.raw_stats.tsv` read counts and lengths |
| `01_qc_raw/` | QC of raw reads | `*_fastqc.html` |
| `02_trimmed/` | fastp | `*.fastp.html`; the trimmed reads `*_R1/R2.trim.fastq.gz` |
| `03_qc_trimmed/` | QC again | `*_fastqc.html` (compare with 01) |
| `04_subsampling/` | subsample | `*.subsample_info.tsv` (genome size used, coverage before/after); `quick_assembly/` if the size was estimated |
| `05_assembly/` | assemble | `unicycler_normal/`, `unicycler_conservative/`, `unicycler_bold/`, `spades/` each with `assembly.fasta`, `assembly.gfa`, the log, `coverage_used.txt`; `*.contigs.tsv` tables |
| `06_bandage/` | Bandage | `*.png` one picture per assembly, `*.graph_info.tsv` |
| `07_mapping/` | map reads, IGV | **`*.alignment_report.html`: open in a browser to see the read pileups** at every mixed position, the termini peaks and the contig ends (no IGV needed); `HOW_TO_VIEW_ALIGNMENTS.txt` explains the IGV desktop route with `*.fasta` + `*.bam` (or `*.igv_session.xml`); `*.coverage.txt`, `*.flagstat.txt`; `*.termini_peaks.tsv` (strongest read-start pile-up per strand) |
| `08_variants/` | mixed positions | `*.mixed_positions.tsv` (position, bases, depth, fraction of reads with the other base); `pilon/*.pilon.changes` |
| `09_phageterm/` | PhageTerm | `*_PhageTerm_report.pdf`, `phageterm.log`. Trust the class only when `termini_signal` in `summary.md` says *present*: PhageTerm needs libraries whose adapters were ligated to the fragment ends (TruSeq-type). Tagmentation (Nextera) libraries show no read-start peak at the termini and PhageTerm then guesses. |
| `10_final/` | finishing | `*.before_rotation.fasta`, `*.after_rotation.fasta`, `rotation_method.txt` |

Plus `results/multiqc/multiqc_report.html` (all QC in one page) and
`results/pipeline_info/` (run report, timeline, trace, logs, `work/` with every step's intermediate files).

All four assemblies are kept so you can compare them; the downstream steps (mapping, Pilon,
PhageTerm, rotation) use **Unicycler normal** unless you pass `--primary unicycler_conservative`
or `--primary spades`. No contig is removed automatically: look at `*.contigs.tsv` and the Bandage
picture and decide yourself which contig is the phage (host DNA and second phages show up as extra contigs).

## 4. What it decides for you (and how to change that)

| decision | default | change with |
|---|---|---|
| maximum read pairs given to the assembler | 500 000 | `--cap_pairs` |
| genome size | estimated: longest contig of a quick SPAdes assembly of 100 000 pairs | `--genome-size <bp>` |
| coverage for assembly | 200x; if the assembly fails, 100x is tried automatically | `--coverage`, `--fallback_coverage` |
| fastp | adapter auto-detect, overlap correction, poly-G, window 4 / Q20, min length 50 | `--cut_mean_quality`, `--min_len` |
| assemblers | Unicycler normal + conservative + bold, SPAdes --isolate | `--unicycler_modes normal`, `--run_spades false` |
| genome start | PhageTerm's re-opened sequence when its class is cos, DTR or headful (pac) **and** the reads show a termini signal; otherwise (UNKNOWN / Mu-like, weak signal, or PhageTerm skipped) dnaapler puts the start upstream of the terminase large subunit (largest ORF if none found) | `--skip-phageterm`, or edit the rule in `modules/finishing.nf` (process ROTATE) |
| which contigs are the phage for termini / PhageTerm / rotation | contigs of at least 10 kb (shorter contigs are kept but not rotated); if none, the longest contig | `--min_phage_contig` |
| memory / CPUs | 8 CPUs, 8 GB per assembly, 12 GB total | `--max_cpus`, `--max_memory '24.GB'` |

## 5. How Nextflow and Docker fit together

**Nextflow** (`main.nf`, `modules/*.nf`) is the conductor: each step is a *process* with declared
inputs, outputs, CPUs and memory. Nextflow runs the steps in the right order, several samples in
parallel, caches every finished step (`-resume`) and writes a report of what ran where.

**Docker** is where the tools live. Each process names a container image, e.g.
`quay.io/biocontainers/fastp:1.4.0--h43da1c4_0` (see `params.containers` in `nextflow.config`).
Nextflow pulls the image the first time and mounts the step's folder into it. Nothing is installed
on your computer, every run uses the same tool versions, and the same files work on a cluster with
Singularity/Apptainer (`-with-apptainer`).

PhageTermVirome has no public image, so `run.sh` builds one locally from
`containers/phageterm/Dockerfile` the first time (about 3 minutes).

Running Nextflow directly, without the wrapper (Nextflow can fetch the pipeline from GitHub itself):

```bash
nextflow run yuliia-faidiuk/phage-assembly-nf --input ~/reads/run1 --outdir results -resume
nextflow run phage-assembly-nf/ --samplesheet samples.csv --genome_size 45000 -resume   # samples.csv: sample,fastq_1,fastq_2
nextflow run phage-assembly-nf/ --help
```
The PhageTerm image must exist locally for the direct form (`docker build -t phage-pipeline/phagetermvirome:4.3 phage-assembly-nf/containers/phageterm`), or pass `--skip_phageterm`.

## 6. Files in this folder

```
run.sh                     the one command students use (checks, sample table, PhageTerm image, nextflow)
main.nf                    the workflow: sample detection + the order of steps
nextflow.config            parameters, container versions, CPU/memory, profiles (test, quick)
modules/reads.nf           lane merge, FastQC, fastp, host removal, read cap, quick assembly, subsample
modules/assembly.nf        Unicycler (3 modes, with 100x fallback), SPAdes, Bandage image, contig table
modules/finishing.nf       bwa mapping + IGV session, alignment report (igv.js), bcftools, Pilon, PhageTerm, dnaapler rotation
modules/report.nf          MultiQC, summary table
bin/doctor.sh              pre-flight check of Docker / Java / Nextflow / memory / disk
bin/detect_samples.py      finds R1/R2 pairs and lanes in a folder, writes the samplesheet
bin/summarize.py           merges the per-step numbers into summary.tsv / summary.md
bin/make_igv_sites.py      chooses the positions shown in the alignment report
bin/termini_peaks.sh       strongest read-start pile-up per strand (termini signal)
containers/phageterm/      Dockerfile for PhageTermVirome 4.3
test_data/get_test_data.sh downloads the two public test phages (T7, phiEc_1)
TROUBLESHOOTING.md         what to do when it stops
```
