#!/usr/bin/env python3
"""Find paired Illumina FASTQ files in a folder and print/write a samplesheet (same naming rule as main.nf)."""
import argparse, os, re, sys
from collections import defaultdict

PAT = re.compile(r"^(.+?)(?:_S\d+)?(?:_L\d{3})?[._]R?([12])(?:_001)?$")

def detect(folder):
    samples = defaultdict(lambda: {"1": [], "2": []}); ignored = []
    for f in sorted(os.listdir(folder)):
        if not re.search(r"\.(fastq|fq)\.gz$", f): continue
        base = re.sub(r"\.(fastq|fq)\.gz$", "", f)
        m = PAT.match(base)
        if not m: ignored.append(f); continue
        samples[m.group(1)][m.group(2)].append(os.path.abspath(os.path.join(folder, f)))
    return samples, ignored

def main():
    ap = argparse.ArgumentParser(); ap.add_argument("folder"); ap.add_argument("--write", help="write samplesheet CSV here")
    a = ap.parse_args()
    if not os.path.isdir(a.folder): sys.exit("ERROR: %s is not a folder" % a.folder)
    samples, ignored = detect(a.folder)
    problems = []
    print("\nSamples found in %s:" % os.path.abspath(a.folder))
    print("  %-20s %6s  %s" % ("sample", "lanes", "files"))
    for s, d in samples.items():
        n1, n2 = len(d["1"]), len(d["2"])
        names = ", ".join(os.path.basename(x) for x in d["1"] + d["2"])
        flag = ""
        if n1 == 0 or n2 == 0: flag = "  <-- missing R1 or R2"; problems.append(s)
        elif n1 != n2: flag = "  <-- %d R1 but %d R2 files" % (n1, n2); problems.append(s)
        print("  %-20s %6d  %s%s" % (s, n1, names, flag))
    for f in ignored: print("  (ignored: %s — name does not look like <sample>_R1/_R2.fastq.gz)" % f)
    if not samples: sys.exit("\nERROR: no *.fastq.gz / *.fq.gz files with R1/R2 (or _1/_2) in their names were found.")
    if problems: sys.exit("\nERROR: fix the samples marked above (this pipeline needs paired Illumina reads).")
    if a.write:
        with open(a.write, "w") as out:
            out.write("sample,fastq_1,fastq_2\n")
            for s, d in samples.items():
                for r1, r2 in zip(d["1"], d["2"]): out.write("%s,%s,%s\n" % (s, r1, r2))
    print()

if __name__ == "__main__":
    main()
