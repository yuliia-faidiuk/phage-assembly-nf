# When the pipeline stops

## The 4-step routine

1. **Read the end of the screen output.** Nextflow prints the failed step, its exit code,
   the command it ran and the tool's last error lines, e.g.

   ```
   ERROR ~ Error executing process > 'UNICYCLER (phageA:normal)'
   Caused by:  Process `UNICYCLER (phageA:normal)` terminated with an error exit status (1)
   Command error:  ...
   Work dir:  /home/me/results/pipeline_info/work/3f/9c1e...
   ```

2. **Open the work dir.** Every step runs in its own folder, which contains:

   | file | what it is |
   |---|---|
   | `.command.sh` | the exact command that ran (you can copy-paste it to rerun by hand) |
   | `.command.log` | everything the tool printed |
   | `.command.err` | only the error stream |
   | `.exitcode` | the exit code (137 = killed for lack of memory) |
   | the inputs | symlinks to the files the step received |
   | the outputs | whatever the tool managed to write, plus the tool's own log (`unicycler.log`, `spades.log`, `pilon.log`, `phageterm.log`) |

   To rerun a step by hand inside its container:

   ```bash
   cd /home/me/results/pipeline_info/work/3f/9c1e...
   bash .command.run        # same container, same inputs; edit .command.sh first if you want to try different options
   ```

3. **Find the message below.** If it is not there, search the message text on the tool's GitHub issues page.

4. **Fix the cause and rerun the same `run.sh` command.** Finished steps are reused.
   If you want to force a step to rerun, delete its work dir or run with `-resume` removed (everything reruns).

The run report `results/pipeline_info/report.html` shows every task with its memory and time,
and `trace.txt` lists each task's status and work dir in one table.

## Messages and fixes

### Before anything runs

| message | cause | fix |
|---|---|---|
| `[FAIL] docker is installed but not running` | Docker Desktop is closed, or WSL integration is off | start Docker Desktop; Settings → Resources → WSL integration → enable your distro. Linux: `sudo service docker start` |
| `permission denied while trying to connect to the Docker daemon socket` | your user is not in the docker group | `sudo usermod -aG docker $USER`, then log out and in (WSL: `wsl --shutdown` from PowerShell, reopen) |
| `nextflow: command not found` right after installing | `~/.local/bin` is not in PATH yet | open a new terminal, or `export PATH=$HOME/.local/bin:$PATH` |
| `Nextflow requires Java 17 or later` | old Java | `sudo apt install openjdk-17-jre-headless` |
| `ERROR: no *.fastq.gz ... were found` | wrong folder, or files not gzipped, or names without R1/R2 | files must end in `.fastq.gz`/`.fq.gz` and contain `_R1`/`_R2` or `_1`/`_2`; `gzip *.fastq` if needed; rename if the facility used other labels |
| `Sample 'x' has R1 files [...] but R2 files []` | single-end data or a missing file | this pipeline is for paired Illumina data; check the download completed |
| `(ignored: file.fastq.gz — name does not look like ...)` | the name has no read number | rename to `<sample>_R1.fastq.gz` / `<sample>_R2.fastq.gz` |
| `Unable to find image ... / pull access denied / TLS handshake timeout` | no internet, or a proxy | the first run needs internet to download ~6 GB; retry, or pull by hand: `docker pull quay.io/biocontainers/fastp:1.4.0--h43da1c4_0` |
| `no space left on device` | Docker images or `work/` filled the disk | `docker system prune` (removes unused images); delete old `results/*/pipeline_info/work` |

### During the run

| step | message | cause | fix |
|---|---|---|---|
| PREPARE_READS | `R1 has N reads but R2 has M reads` | files from different samples/lanes were paired, or a download was truncated | check the sample table printed at the start; re-download; `gzip -t file.fastq.gz` tests for corruption |
| FASTP | `ERROR: ... is not a valid gzip file` | truncated download or wrong extension | `gzip -t` the file; re-download |
| any | exit code **137** (or `Killed`, `Out of memory`) | not enough RAM. The pipeline retries once with 11 GB; if that also fails, your machine has less | close other programs; WSL: raise memory in `%UserProfile%\.wslconfig` (`[wsl2]` / `memory=12GB`) then `wsl --shutdown`; run with `--max_memory '10.GB'` so Nextflow runs one assembly at a time; use `--quick` |
| QUICK_ASSEMBLY | `quick SPAdes assembly failed` | reads too few, too short, or not a DNA phage; or memory | give `--genome-size <bp>` to skip the estimate; check `02_trimmed/*.fastp.html` for read length and quality |
| QUICK_ASSEMBLY | `WARNING: longest quick contig is only N bp` | fragmented quick assembly, so the genome size (= longest contig) is too small and the subsample too shallow | pass `--genome-size` with the expected size |
| UNICYCLER | `failed at both 200x and 100x` | SPAdes crashed inside Unicycler (see `unicycler.log` in the work dir): usually memory, or reads with very uneven coverage, or almost no phage reads | look at `05_assembly/spades/` (plain SPAdes may have worked); try `--primary spades`; check % reads mapped later; remove host reads with `--host` |
| UNICYCLER | `Error: could not find SPAdes` / `spades.py: command not found` | wrong container | do not change `params.containers.unicycler`; `docker pull` it again |
| UNICYCLER | assembly has many contigs, none circular | low or uneven coverage, host contamination, two phages | compare with `unicycler_conservative`/`bold`; check `06_bandage/*.png`; raise `--coverage 300` if `coverage_available` in `summary.md` allows; use `--host` |
| SPADES | `Error == 0 ... out of memory` | memory | as above (exit 137) |
| BANDAGE_IMAGE | `could not connect to display` / `qt.qpa.plugin` | missing offscreen platform | already set (`QT_QPA_PLATFORM=offscreen`); if it still fails, the graph is in `05_assembly/*/assembly.gfa`, open it in Bandage on your PC |
| MAP_READS | `[E::bwa_idx_load]` / `fail to locate the index` | the index step failed (disk full?) | check disk, rerun |
| IGV_REPORT | `create_report` error | very deep BAM or odd contig names | the BAM and session file are still there; open them in IGV desktop as described in `07_mapping/HOW_TO_VIEW_ALIGNMENTS.txt` |
| VARIANTS | empty table | no disagreement between reads and contig: that is good | nothing |
| PILON | `java.lang.OutOfMemoryError` | memory | run with `--max_memory` higher if your machine has it; or skip Pilon by using the unpolished `07_mapping/*.fasta` |
| PHAGETERM | marked FAILED but the pipeline continued | PhageTerm is optional; typical causes: reference with >1 contig that confuses it, very few reads mapping, or the local image missing | read `09_phageterm/phageterm.log`; run it by hand on the longest contig only; `docker build -t phage-pipeline/phagetermvirome:4.3 phage-assembly-nf/containers/phageterm` if the image is missing |
| ROTATE | `dnaapler failed` | no terminase hit and no ORFs (not a phage? contigs too short) | `10_final/*.before_rotation.fasta` is still a valid genome; rotate by hand with `seqkit restart -i <pos>` as in the homework |
| SUMMARY | `INCOMPLETE: missing ...` in `summary.md` | a step above failed for that sample | fix that step; the summary rebuilds on rerun |

### Interpreting a successful run that looks wrong

| symptom in `summary.md` | likely cause | what to do |
|---|---|---|
| `pct_reads_kept` < 60 | bad run, wrong adapters, or very short inserts | look at `02_trimmed/*.fastp.html` → "adapters" and "insert size" |
| `coverage_available` < 50 | too few reads for this genome size, or the genome-size estimate is wrong | check `04_subsampling/quick_assembly/quick_contigs.fasta` |
| `pct_mapped` < 80 | host DNA or another phage in the sample; the unmapped reads are something else | `--host`; look at the extra contigs in `05_assembly/*/assembly.fasta` |
| `mixed_positions_20_80pct` > 0 | positions where 20–80 % of reads disagree: real mixed population, or a repeat | open `*.igv_session.xml` in IGV at that position (homework step 8) |
| `termini_signal` = weak | the reads carry no fragment-end information (tagmentation library, PCR library) or the phage has no fixed ends (headful) | the PhageTerm class cannot be trusted; the genome is still correct. The pipeline then ignores PhageTerm and rotates to the terminase (see `10_final/rotation_method.txt`). Ask the sequencing facility which library kit was used; a TruSeq-type (ligation) prep is needed to see termini |
| `termini_signal` = present but PhageTerm says UNKNOWN | the peak is at a contig end, or coverage is very uneven | look at `07_mapping/*.termini_peaks.tsv`: the peak positions are the physical ends; the distance between the forward and the reverse peak is the DTR length |
| `longest_circular` = no but one contig of the expected size | ends of a linear genome with fixed termini (cos/DTR) or low coverage at the ends | PhageTerm report; this is expected for many phages |
| `pilon_changes` > 20 | the assembly had errors, often low coverage | check mean depth; consider higher `--coverage` |

## Nextflow housekeeping

```bash
nextflow log                                  # list previous runs (run from results/pipeline_info)
nextflow log <run name> -f name,status,workdir  # which task ran where
nextflow clean -f -k <run name>               # delete that run's work files
rm -r results/pipeline_info/work              # delete all intermediate files once you are done
docker system prune                           # remove unused images (frees several GB)
```
