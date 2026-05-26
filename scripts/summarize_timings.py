#!/usr/bin/env python3
import csv
import pathlib
import sys


def parse_time_file(path: pathlib.Path):
    rows = []
    with path.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            rows.append(row)
    return rows


def main():
    if len(sys.argv) != 3:
        print("Usage: summarize_timings.py <results_dir> <out_tsv>", file=sys.stderr)
        return 1

    results_dir = pathlib.Path(sys.argv[1])
    out_tsv = pathlib.Path(sys.argv[2])
    out_tsv.parent.mkdir(parents=True, exist_ok=True)

    out_rows = []
    for time_file in sorted(results_dir.rglob("*.time.tsv")):
        for row in parse_time_file(time_file):
            out_rows.append(
                {
                    "run_id": time_file.stem.replace(".time", ""),
                    "tool": row.get("tool", ""),
                    "elapsed_sec": row.get("elapsed_sec", ""),
                    "cpu_user_sec": row.get("cpu_user_sec", ""),
                    "cpu_sys_sec": row.get("cpu_sys_sec", ""),
                    "cpu_total_sec": row.get("cpu_total_sec", ""),
                    "max_rss_kb": row.get("max_rss_kb", ""),
                    "time_file": str(time_file),
                }
            )

    with out_tsv.open("w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=[
                "run_id",
                "tool",
                "elapsed_sec",
                "cpu_user_sec",
                "cpu_sys_sec",
                "cpu_total_sec",
                "max_rss_kb",
                "time_file",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(out_rows)

    print(f"Wrote {out_tsv} ({len(out_rows)} rows)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
