# Changelog

Versions follow [semantic versioning](https://semver.org): MAJOR.MINOR.PATCH.
MAJOR = outputs or parameters change in a way that breaks existing use; MINOR = new feature, old runs still work;
PATCH = bug fix only. Run a specific version with `nextflow run yuliia-faidiuk/phage-assembly-nf -r v1.1.0`.

## [1.1.1] - 2026-10-18
### Changed
- The version is printed at the start of every run (wrapper and Nextflow log); README explains how to update.

## [1.1.0] - 2026-10-18
### Added
- `07_mapping/<sample>.alignment_report.html`: self-contained read-pileup viewer (igv-reports) at every mixed
  position, the termini peaks, contig ends and three windows per contig; no IGV or internet needed.
- `07_mapping/HOW_TO_VIEW_ALIGNMENTS.txt`: manual steps for IGV desktop (incl. WSL).

## [1.0.0] - 2026-10-18
First release. Raw Illumina FASTQ -> FastQC -> fastp -> read cap -> genome-size estimate -> 200x subsample
(100x fallback) -> Unicycler normal/conservative/bold + SPAdes -> Bandage -> bwa mapping -> bcftools mixed
positions -> Pilon -> PhageTerm (optional) -> rotation (PhageTerm termini or dnaapler terL) -> MultiQC + summary.
Termini-signal check, student wrapper `run.sh`, pre-flight `doctor.sh`, test profile on two public phages.
