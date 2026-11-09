#!/usr/bin/env python3
"""North/South sample lists for the gene-wise Fst analysis (Figure 3), rebuilt from the public metadata table.

Individuals: Admixture population == "East Siberia" and ploidy == 2 (124 individuals; the 6 tetraploid and 2 "unsure" ones are dropped).
Groups: "South" = latitude < 57 (7 individuals, all at 56.76-56.82 N), "North" = the rest (117, 63.1-75.3 N).
The VCF calls BAM12.1-1 "BAM_12.1-1"; that name is written to the outputs.

  make_northsouth_groups.py [METADATA.tsv|URL] [OUT_DIR]

METADATA defaults to https://lyrata.mpipz.mpg.de/data/lyrata_metadata.tsv (latin-1 encoded), OUT_DIR to results/staging/real.
Writes OUT_DIR/northsouth_groups_124.tsv (sample, group; input of piawka/pixy) and scripts/northsouth_vcf/samples_ns.txt
(sample list for the GT-only VCF extraction) when OUT_DIR is the default. Pure standard library.
"""
import csv, io, pathlib, sys, urllib.request

ROOT = pathlib.Path(__file__).resolve().parents[1]
URL = "https://lyrata.mpipz.mpg.de/data/lyrata_metadata.tsv"
LAT_SOUTH = 57.0
VCF_NAME = {"BAM12.1-1": "BAM_12.1-1"}

src = sys.argv[1] if len(sys.argv) > 1 else URL
out = pathlib.Path(sys.argv[2]) if len(sys.argv) > 2 else ROOT / "results" / "staging" / "real"
raw = urllib.request.urlopen(src).read() if src.startswith("http") else pathlib.Path(src).read_bytes()
rows = list(csv.DictReader(io.StringIO(raw.decode("latin-1")), delimiter="\t"))
es = [r for r in rows if r["Admixture population"] == "East Siberia" and r["ploidy"] == "2"]
groups = [(VCF_NAME.get(r["Sample_name"], r["Sample_name"]), "South" if float(r["latitude"]) < LAT_SOUTH else "North") for r in es]
n_south = sum(g == "South" for _, g in groups)
assert (len(groups), n_south) == (124, 7), (len(groups), n_south)
out.mkdir(parents=True, exist_ok=True)
(out / "northsouth_groups_124.tsv").write_text("".join(f"{s}\t{g}\n" for s, g in groups))
if len(sys.argv) <= 2:
    (ROOT / "scripts" / "northsouth_vcf" / "samples_ns.txt").write_text("".join(f"{s}\n" for s, _ in groups))
print(f"{len(groups)} individuals: {len(groups) - n_south} North, {n_south} South (latitude < {LAT_SOUTH:g})")
