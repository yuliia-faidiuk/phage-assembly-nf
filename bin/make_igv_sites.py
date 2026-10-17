#!/usr/bin/env python3
"""Build the list of positions shown in the alignment report (BED with a name column).
Sites: every variant / mixed position (VCF), the two termini peaks, start and end of each long contig,
and three evenly spaced windows per long contig so the report is never empty."""
import argparse, gzip

ap = argparse.ArgumentParser()
ap.add_argument("--fai", required=True); ap.add_argument("--termini", required=True)
ap.add_argument("--vcf", required=True); ap.add_argument("--min-contig", type=int, default=10000)
ap.add_argument("--max-variants", type=int, default=100); ap.add_argument("--out", required=True)
a = ap.parse_args()

sites = []  # (contig, pos0, name)
lengths = {}
for line in open(a.fai):
    c, l = line.split("\t")[:2]; lengths[c] = int(l)
long_contigs = [c for c, l in lengths.items() if l >= a.min_contig] or [max(lengths, key=lengths.get)]

n = 0
for line in gzip.open(a.vcf, "rt"):
    if line.startswith("#"): continue
    f = line.split("\t")
    n += 1
    if n <= a.max_variants:
        sites.append((f[0], int(f[1]) - 1, "variant_%s_%s>%s" % (f[1], f[3], f[4])))
if n > a.max_variants:
    print("note: %d variants, only the first %d are in the report" % (n, a.max_variants))

for line in open(a.termini):
    f = line.rstrip("\n").split("\t")
    if len(f) >= 6 and f[0] in ("forward_start", "reverse_end") and f[1] in lengths:
        sites.append((f[1], int(f[2]) - 1, "termini_peak_%s_%s_of_depth" % (f[0], f[5].replace(".", "p"))))

for c in long_contigs:
    L = lengths[c]
    sites.append((c, 0, "contig_%s_start" % c))
    sites.append((c, L - 1, "contig_%s_end" % c))
    for i in (1, 2, 3):
        sites.append((c, L * i // 4, "contig_%s_window_%d" % (c, i)))

with open(a.out, "w") as out:
    for c, p, name in sites:
        out.write("%s\t%d\t%d\t%s\n" % (c, max(0, p), max(0, p) + 1, name))
print("%d sites written" % len(sites))
