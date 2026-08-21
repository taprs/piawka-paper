#!/usr/bin/env python3
"""Bundle Supplementary Tables 1 & 2 (the North/South gene-wise Fst rank
tables produced by scripts/analyze_northsouth_fst_ranks.py) into a single
.xlsx workbook, one sheet per table.

Usage:
  conda run -n piawka-paper python scripts/collect_fst_tables.py [output.xlsx]

Defaults to results/tables/piawka_supplementary_data.xlsx.
"""
import csv
import pathlib
import sys

from openpyxl import Workbook

ROOT = pathlib.Path(__file__).resolve().parents[1]
TABLES = ROOT / "results" / "tables"

SHEETS = [
    ("Supp Table 1 - GWAS SNP genes", TABLES / "northsouth_fst_ranks_gwas_snp_genes.tsv"),
    ("Supp Table 2 - GO2000028 genes", TABLES / "northsouth_fst_ranks_go2000028_genes.tsv"),
]

DEFAULT_OUTPUT = TABLES / "piawka_supplementary_data.xlsx"


def main():
    output = pathlib.Path(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT_OUTPUT

    wb = Workbook()
    wb.remove(wb.active)

    for sheet_name, tsv_path in SHEETS:
        if not tsv_path.exists():
            raise FileNotFoundError(f"missing input table: {tsv_path}")
        ws = wb.create_sheet(title=sheet_name)
        with tsv_path.open(newline="") as fh:
            for row in csv.reader(fh, delimiter="\t"):
                ws.append(row)
        for col_cells in ws.columns:
            width = max(len(str(c.value)) for c in col_cells if c.value is not None)
            ws.column_dimensions[col_cells[0].column_letter].width = min(width + 2, 40)

    output.parent.mkdir(parents=True, exist_ok=True)
    wb.save(output)
    print(f"wrote {output}")


if __name__ == "__main__":
    main()
