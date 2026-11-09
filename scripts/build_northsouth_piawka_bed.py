#!/usr/bin/env python3
"""Assemble results/staging/real/northsouth_piawka.bed (input of the gene-wise Fst rank tables and the jackknife site counts).

  piawka calc --persite ... -s fst | build_northsouth_piawka_bed.py GENES.bed FST.bed OUT.bed

FST.bed is the `default` piawka run of run_a3_northsouth.sh (rows `fst`, copied unchanged). Every row of the --persite output is one site
that contributed to a gene's Fst (stdin); the sites falling in each GENES.bed interval are counted and written as `segr2` rows
(North, South, value = numerator = number of segregating sites usable by the pair, multiallelic included, denominator 1; 0 for genes
without such sites). Pure standard library.
The `lines2` rows (sites usable by the pair, invariant included) of the original file came from a build that is no longer available;
nothing in this repository reads them.
"""
import collections, csv, sys

genes_bed, fst_bed, out_bed = sys.argv[1:4]
by_locus = collections.defaultdict(list)           # (chrom, locus) -> [(start, end)]
genes = []
for r in csv.reader(open(genes_bed), delimiter="\t"):
    if r and not r[0].startswith("#"):
        genes.append((r[0], int(r[1]), int(r[2]), r[3])); by_locus[(r[0], r[3])].append((int(r[1]), int(r[2])))
count = collections.Counter(); seen = set()
for r in csv.reader(sys.stdin, delimiter="\t"):
    if r[0].startswith("#") or r[6] != "fst": continue
    chrom, s, e, locus = r[0], int(r[1]), int(r[2]), r[3]
    if (chrom, s, locus) in seen: continue      # overlapping same-name intervals print their shared sites once per interval
    seen.add((chrom, s, locus))
    for a, b in by_locus.get((chrom, locus), ()):
        if a <= s and e <= b: count[(chrom, locus, a, b)] += 1
with open(out_bed, "w") as out:
    fst = open(fst_bed); out.write(fst.readline())
    n = 0
    for line in fst:
        if line.split("\t")[6] == "fst": out.write(line); n += 1
    for chrom, a, b, locus in genes:
        c = count[(chrom, locus, a, b)]
        out.write(f"{chrom}\t{a}\t{b}\t{locus}\tNorth\tSouth\tsegr2\t{c}\t{c}\t1\n")
print(f"{n} fst rows, {len(genes)} segr2 rows -> {out_bed}", file=sys.stderr)
