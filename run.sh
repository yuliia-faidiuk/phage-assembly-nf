#!/usr/bin/env bash
# One-command runner for the phage assembly pipeline.
#
#   bash phage-assembly-nf/run.sh <folder_with_fastq.gz> [--genome-size 40000] [--outdir results] [--host host.fasta] [--quick] [--yes]
#   bash phage-assembly-nf/run.sh --test   # runs on the two public phage runs in test_data/ (downloads ~40 MB once)
#
# Everything the pipeline does is described in README.md; errors are explained in TROUBLESHOOTING.md (same folder as this script).
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
INPUT=""; OUTDIR=""; EXTRA=(); YES=0; TEST=0; SKIP_PT=0; PROFILE=""
while [ $# -gt 0 ]; do
    case "$1" in
        --test) TEST=1; shift;;
        --genome-size|--genome_size) EXTRA+=(--genome_size "$2"); shift 2;;
        --outdir) OUTDIR="$2"; shift 2;;
        --host) EXTRA+=(--host_fasta "$(readlink -f "$2")"); shift 2;;
        --quick) PROFILE="quick"; shift;;
        --skip-phageterm) SKIP_PT=1; EXTRA+=(--skip_phageterm); shift;;
        --yes|-y) YES=1; shift;;
        -h|--help) sed -n '2,8p' "$0"; echo; nextflow run "$HERE" --help 2>/dev/null | sed -n '/Usage/,$p'; exit 0;;
        --*) EXTRA+=("$1" "$2"); shift 2;;      # any other --param value is passed to nextflow
        *) INPUT="$1"; shift;;
    esac
done

echo "== 1/4  pre-flight check"
bash "$HERE/bin/doctor.sh" "${INPUT:-$HERE}" || exit 1

if [ "$TEST" = 1 ]; then
    INPUT="$HERE/test_data"; OUTDIR="${OUTDIR:-$INPUT/results_test}"
    echo "== test data"; bash "$INPUT/get_test_data.sh" || { echo "test data download failed (internet needed the first time)"; exit 1; }
fi
[ -n "$INPUT" ] || { echo "ERROR: give the folder with your .fastq.gz files, e.g.  bash phage-assembly-nf/run.sh ~/reads/run1"; exit 1; }
[ -d "$INPUT" ] || { echo "ERROR: $INPUT is not a folder"; exit 1; }
INPUT="$(readlink -f "$INPUT")"
OUTDIR="${OUTDIR:-$PWD/results}"; mkdir -p "$OUTDIR/pipeline_info"

echo "== 2/4  samples"
python3 "$HERE/bin/detect_samples.py" "$INPUT" --write "$OUTDIR/pipeline_info/samplesheet.csv" || exit 1
if [ "$YES" = 0 ] && [ -t 0 ]; then
    read -r -p "Run the pipeline on these samples? [Y/n] " ans; case "${ans:-y}" in y|Y) ;; *) echo "stopped"; exit 0;; esac
fi

echo "== 3/4  PhageTerm container"
if [ "$SKIP_PT" = 0 ]; then
    if docker image inspect phage-pipeline/phagetermvirome:4.3 >/dev/null 2>&1; then echo "  image present"
    else
        echo "  building phage-pipeline/phagetermvirome:4.3 (one time, ~3 min)"
        docker build -q -t phage-pipeline/phagetermvirome:4.3 "$HERE/containers/phageterm" > "$OUTDIR/pipeline_info/phageterm_build.log" 2>&1 \
          || { echo "  build failed (see $OUTDIR/pipeline_info/phageterm_build.log); continuing without PhageTerm"; EXTRA+=(--skip_phageterm); }
    fi
fi

echo "== 4/4  nextflow"
cd "$OUTDIR/pipeline_info" || exit 1    # keeps .nextflow.log and work/ next to the results, not in the reads folder
CMD=(nextflow run "$HERE" --samplesheet "$OUTDIR/pipeline_info/samplesheet.csv" --outdir "$OUTDIR" -resume -work-dir "$OUTDIR/pipeline_info/work" "${EXTRA[@]}")
[ -n "$PROFILE" ] && CMD+=(-profile "$PROFILE")
echo "  ${CMD[*]}"; echo
"${CMD[@]}"; rc=$?
echo
if [ $rc -eq 0 ]; then
    echo "Done. Open $OUTDIR/summary.md, then the folders in $OUTDIR/<sample>/."
    echo "To free disk space when you are happy with the results:  rm -r $OUTDIR/pipeline_info/work"
else
    echo "The pipeline stopped with an error (exit $rc)."
    echo "  1. Read the lines above that start with 'ERROR' or 'Caused by'."
    echo "  2. Look up the message in $HERE/TROUBLESHOOTING.md."
    echo "  3. Full log: $OUTDIR/pipeline_info/.nextflow.log"
    echo "  4. After fixing, run exactly the same command again: finished steps are not repeated."
fi
exit $rc
