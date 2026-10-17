#!/usr/bin/env python3
"""Merge the per-step key/value summary files into one table per sample (summary.tsv + summary.md)."""
import argparse, csv, glob, json, os, re
from collections import OrderedDict

ORDER = ["raw_pairs", "raw_mean_len", "trimmed_pairs", "pct_reads_kept", "mean_len_after_trim", "q30_after_trim_pct", "adapter_detected",
         "host_pairs_removed", "host_pct_removed", "pairs_after_cap", "read_cap_applied", "genome_size_bp", "genome_size_source", "quick_sum_contigs_ge10kb",
         "coverage_available", "main_pairs", "main_coverage_target"]
for _a in ["unicycler_normal", "unicycler_conservative", "unicycler_bold", "spades"]:
    ORDER += ["%s:%s" % (_a, k) for k in ["n_contigs", "longest_bp", "longest_circular", "longest_depth", "contigs_ge10kb", "coverage_used"]]
ORDER += ["mapped_to", "pct_mapped", "pct_properly_paired", "mean_depth", "pct_genome_covered", "termini_peak_forward", "termini_peak_reverse", "termini_signal",
          "variant_positions", "mixed_positions_20_80pct", "pilon_changes",
          "phageterm_class", "phageterm_reopened_sequence", "final_contigs", "final_total_bp", "rotation"]
HIDE = {"raw_files"}

# steps that must have produced a summary file; (name used in the status column, file suffix)
STEPS = [("raw reads", "00_raw"), ("read cap", "04_cap"), ("subsample", "04_subsample"), ("unicycler", "05_unicycler_"),
         ("mapping", "07_mapping"), ("variants", "08_variants"), ("pilon", "08_pilon"), ("final", "10_final")]   # phageterm is optional

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--summaries", required=True); ap.add_argument("--fastp", required=True)
    ap.add_argument("--samplesheet", required=True); ap.add_argument("--tsv", required=True); ap.add_argument("--md", required=True)
    ap.add_argument("--skip-phageterm", action="store_true")
    a = ap.parse_args()

    samples = []
    with open(a.samplesheet) as fh:
        for row in csv.DictReader(fh):
            if row["sample"] not in samples: samples.append(row["sample"])

    data = {s: OrderedDict() for s in samples}
    seen = {s: set() for s in samples}
    for f in sorted(glob.glob(os.path.join(a.summaries, "*.summary.tsv"))):
        base = os.path.basename(f)
        m = re.match(r"(.+?)\.(\d\d_[A-Za-z_]+?)\.summary\.tsv$", base)
        if not m: continue
        sid, step = m.group(1), m.group(2)
        if sid not in data: continue
        seen[sid].add(step)
        for line in open(f):
            if "\t" in line:
                k, v = line.rstrip("\n").split("\t", 1)
                data[sid][k] = v

    for f in glob.glob(os.path.join(a.fastp, "*.fastp.json")):
        sid = os.path.basename(f)[:-len(".fastp.json")]
        if sid not in data: continue
        try:
            j = json.load(open(f)); b, af = j["summary"]["before_filtering"], j["summary"]["after_filtering"]
            data[sid]["trimmed_pairs"] = str(af["total_reads"] // 2)
            data[sid]["pct_reads_kept"] = "%.1f" % (100.0 * af["total_reads"] / max(1, b["total_reads"]))
            data[sid]["mean_len_after_trim"] = str(af.get("read1_mean_length", "NA"))
            data[sid]["q30_after_trim_pct"] = "%.1f" % (100.0 * af["q30_rate"])
            ad = j.get("adapter_cutting", {}).get("read1_adapter_sequence", "none")
            data[sid]["adapter_detected"] = "none detected" if ad in ("unspecified", "none", "") else ad
        except Exception as e:
            data[sid]["trimmed_pairs"] = "NA (%s)" % e

    # status: which steps are missing
    for s in samples:
        if "phageterm_class" not in data[s]:
            data[s]["phageterm_class"] = "skipped (--skip-phageterm)" if a.skip_phageterm else "FAILED (see 09_phageterm/phageterm.log)"
        missing = [name for name, suf in STEPS if not any(st.startswith(suf) for st in seen[s])]
        data[s]["status"] = "OK" if not missing else "INCOMPLETE: missing " + ", ".join(missing)

    keys = [k for k in ORDER if any(k in data[s] for s in samples)]
    keys += sorted(k for s in samples for k in data[s] if k not in keys and k != "status" and k not in HIDE)
    keys = list(OrderedDict.fromkeys(keys)) + ["status"]

    with open(a.tsv, "w") as out:
        out.write("sample\t" + "\t".join(keys) + "\n")
        for s in samples:
            out.write(s + "\t" + "\t".join(data[s].get(k, "") for k in keys) + "\n")

    with open(a.md, "w") as out:
        out.write("# Pipeline summary\n\n")
        out.write("| metric | " + " | ".join(samples) + " |\n|---|" + "---|" * len(samples) + "\n")
        for k in keys:
            out.write("| " + k.replace(":", " / ") + " | " + " | ".join(data[s].get(k, "") for s in samples) + " |\n")
        out.write("\nFolders per sample: 00_raw_reads, 01_qc_raw, 02_trimmed, 03_qc_trimmed, 04_subsampling, 05_assembly, 06_bandage, "
                  "07_mapping, 08_variants (+pilon), 09_phageterm, 10_final (before/after rotation).\n")

if __name__ == "__main__":
    main()
