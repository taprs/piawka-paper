#!/usr/bin/env python3
"""A3 follow-up: per-site contributions to piawka's Hudson Fst in the genes where piawka and pixy percentiles differ
(non-overlapping jackknife CIs, GO:2000028, piawka percentile >= 75). Pure standard library.

  analyze_a3_skew_sites.py select OUT.bed   write the gene BED (same rows as results/staging/real/genes.bed) for scripts/run_a3_skew_sites.sh
  analyze_a3_skew_sites.py                  read results/a3_skew_sites/{persite.bed,sites.csq.vcf} and write
    results/tables/northsouth_skew_sites_piawka_only.tsv  one row per site piawka uses but pixy cannot (REF/ALT not a single A/C/G/T:
                                                          multi-ALT, indel/MNP, spanning deletion '*'), sorted by chromosome and
                                                          position, all analysed genes; alt_observed = the ALT allele(s) seen in
                                                          North/South (used to underline them in the xlsx)
    results/tables/northsouth_skew_gene_summary.tsv       per gene: Fst, Fst restricted to pixy-usable sites, pixy Fst, csq split

fst_effect = gene Fst - gene Fst without that site (leave-one-out; positive = the site raises Fst). The pixy-usable subset (single A/C/G/T REF and ALT) reproduces pixy's Fst and no_snps
in these genes but is not pixy's own filter code.
"""
import csv, collections, pathlib, re, sys
ROOT = pathlib.Path(__file__).resolve().parents[1]
TB = ROOT / "results" / "tables"; O = ROOT / "results" / "a3_skew_sites"
JK = TB / "northsouth_jackknife_ci_go2000028_genes.tsv"
BASE = re.compile("^[ACGT]$")

def select():
    out = []
    for r in csv.DictReader(open(JK), delimiter="\t"):
        pk_lo, pk_hi, px_lo, px_hi = (float(r[k]) for k in ("piawka_ci_lo", "piawka_ci_hi", "pixy_ci_lo", "pixy_ci_hi"))
        if float(r["piawka_percentile"]) >= 75 and (pk_lo > px_hi or px_lo > pk_hi): out.append(r)
    return out

if len(sys.argv) > 1 and sys.argv[1] == "select":
    genes = select(); rows = []
    for r in genes:
        chrom, span = r["interval"].split(":"); s, e = span.split("-")
        rows.append((chrom, s, e, r["gene_lyrata"]))
    pathlib.Path(sys.argv[2]).write_text("".join("\t".join(x) + "\n" for x in rows))
    print("selected", [x[3] for x in rows]); sys.exit(0)

VALID = re.compile(r"^[A-Za-z]$")  # piawka only uses single-character alleles (not '*')

def observed_alts(inf):
    """ALT allele(s) with AC>0 in North or South that piawka can use (single character), in VCF order."""
    ac = [[int(x) for x in inf.get(k, "").split(",") if x] for k in ("AC_North", "AC_South")]
    return [a for i, a in enumerate(inf["_alts"]) if VALID.match(a) and any(len(c) > i and c[i] > 0 for c in ac)]

def polymorphic(ref, inf):
    """True if the record has >= 2 piawka-usable alleles observed (REF counts when single-character and seen)."""
    ref_seen = VALID.match(ref) and int(inf.get("AN", 0)) - sum(int(x) for x in inf.get("AC", "0").split(",")) > 0
    return len(observed_alts(inf)) + bool(ref_seen) >= 2

def observed_csq(inf):
    """csq category(ies) for the ALT allele(s) with AC>0 in North or South; see module docstring."""
    seen_alts = set(observed_alts(inf))
    cats = []
    for e in inf.get("BCSQ", "").split(","):
        if not e or e[0] == "@": continue
        f = e.split("|")
        if len(f) >= 7:  # allele-specific entry: dna_change like 15494078G>GGACCA
            if f[6].split("+")[0].split(">")[-1] not in seen_alts: continue
        label = f[0]
        if f[0] == "missense" and len(f) > 5:  # BCSQ amino_acid_change "16D>16E" -> "missense (D16E)"
            aa = [a + p + b for p, a, q, b in re.findall(r"(\d+)([A-Z*])>(\d+)([A-Z*])", f[5])]
            if aa: label += " (" + "+".join(aa) + ")"
        if label not in cats: cats.append(label)
    return ",".join(cats) if cats else "none"
chrom_key = lambda c: [int(x) if x.isdigit() else x for x in re.split(r"(\d+)", c)]
genes = {r["gene_lyrata"]: r for r in select()}
info = collections.defaultdict(list)  # (chrom,pos) -> [(ref, alt, INFO dict)], several records can share a position
for l in open(O / "sites.csq.vcf"):
    if l[0] == "#": continue
    f = l.rstrip("\n").split("\t")
    d = dict(x.split("=", 1) for x in f[7].split(";") if "=" in x); d["_alts"] = f[4].split(",")
    info[(f[0], int(f[1]))].append((f[3], f[4], d))
sites = collections.defaultdict(list)
for r in csv.reader(open(O / "persite.bed"), delimiter="\t"):
    if r[0][0] == "#" or r[6] != "fst" or r[7] in ("NA", "nan", ""): continue
    pos = int(r[2])
    # piawka --persite writes one row per position; pick the VCF record that is polymorphic in usable alleles
    cand = [x for x in info[(r[0], pos)] if polymorphic(x[0], x[2])]
    assert len(cand) == 1, ("ambiguous or missing VCF record", r[0], pos, [x[:2] for x in cand])
    ref, alt, inf = cand[0]
    sites[r[3]].append(dict(chrom=r[0], pos=pos, ref=ref, alt=alt, inf=inf, num=float(r[8]), den=float(r[9]),
                            pixy=bool(BASE.match(ref)) and bool(BASE.match(alt))))
site_rows, summ_rows = [], []
for g, rs in sites.items():
    N = sum(r["num"] for r in rs); D = sum(r["den"] for r in rs); F = N / D
    pu = [r for r in rs if r["pixy"]]; Fp = sum(r["num"] for r in pu) / sum(r["den"] for r in pu)
    # sanity: per-site sums must reproduce the jackknife table's piawka Fst
    assert abs(F - float(genes[g]["piawka_fst"])) < 5e-4, (g, F, genes[g]["piawka_fst"])
    mass = collections.Counter()
    for r in rs:
        r["eff"] = F - (N - r["num"]) / (D - r["den"])
        r["csq"] = observed_csq(r["inf"])
        if not r["pixy"]: mass[re.sub(r" \([^)]*\)", "", r["csq"])] += r["num"]  # category only, amino acid change dropped
    for r in (r for r in rs if not r["pixy"]):
        site_rows.append([g, genes[g]["common_name"], r["chrom"], r["pos"], r["ref"], r["alt"],
                          "%.5f" % r["num"], "%.5f" % r["den"], "%+.5f" % r["eff"], r["csq"], ",".join(observed_alts(r["inf"]))])
    summ_rows.append([g, genes[g]["common_name"], len(rs), len(pu), "%.4f" % F, "%.4f" % Fp, genes[g]["pixy_fst"],
                      "; ".join("%s=%.2f" % kv for kv in mass.most_common())])
site_rows.sort(key=lambda r: (chrom_key(r[2]), r[3]))
with open(TB / "northsouth_skew_sites_piawka_only.tsv", "w") as f:
    f.write("\t".join(["gene_lyrata", "common_name", "chrom", "pos", "ref", "alt",
                       "hudson_num", "hudson_den", "fst_effect", "csq", "alt_observed"]) + "\n")
    f.writelines("\t".join(map(str, r)) + "\n" for r in site_rows)
with open(TB / "northsouth_skew_gene_summary.tsv", "w") as f:
    f.write("\t".join(["gene_lyrata", "common_name", "piawka_sites", "pixy_usable_sites", "piawka_fst", "piawka_fst_pixy_usable_sites",
                       "pixy_fst", "numerator_of_piawka_only_sites_by_csq"]) + "\n")
    f.writelines("\t".join(map(str, r)) + "\n" for r in summ_rows)
print(len(site_rows), "piawka-only sites in", len(summ_rows), "genes")
