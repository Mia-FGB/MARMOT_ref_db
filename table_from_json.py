#!/usr/bin/env python3

# usage python table_from_json.py <input.json> -o <output.tsv>


import sys
import json
import argparse
from pathlib import Path
from typing import Iterable, Dict, Any, List, Optional

# Fields we want in the output
OUT_HEADERS = ["Organism name", "TaxID", "Accession number"]

# Mapping from JSON keys to output headers
KEY_MAP = {
    "organism_name": "Organism name",
    "taxid": "TaxID",
    "assembly_accession": "Accession number",
}

def load_json(path: Optional[str]) -> Any:
    """Load JSON from a file path or stdin."""
    if path:
        with open(path, "r", encoding="utf-8") as f:
            return json.load(f)
    # read from stdin
    return json.load(sys.stdin)

def normalise_records(data: Any) -> Iterable[Dict[str, Any]]:
    """
    Accept either a single object, a list of objects, or newline delimited JSON.
    Yield dicts.
    """
    if isinstance(data, dict):
        yield data
    elif isinstance(data, list):
        for item in data:
            if isinstance(item, dict):
                yield item
            else:
                raise ValueError("List contains a non object item")
    else:
        raise ValueError("Top level JSON must be an object or an array of objects")

def row_from_record(rec: Dict[str, Any]) -> List[str]:
    """
    Build a TSV row in the output order.
    Missing fields become empty strings.
    """
    out = {h: "" for h in OUT_HEADERS}
    for src_key, dst_header in KEY_MAP.items():
        if src_key in rec and rec[src_key] is not None:
            out[dst_header] = str(rec[src_key])
    return [out[h] for h in OUT_HEADERS]

def write_tsv(rows: Iterable[List[str]], outpath: Optional[str]) -> None:
    import csv
    if outpath:
        out_file = open(outpath, "w", newline="", encoding="utf-8")
    else:
        out_file = sys.stdout
    writer = csv.writer(out_file, delimiter="\t", lineterminator="\n")
    writer.writerow(OUT_HEADERS)
    for r in rows:
        writer.writerow(r)
    if out_file is not sys.stdout:
        out_file.close()

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Convert genome JSON to a TSV with Organism name, TaxID, Accession number."
    )
    parser.add_argument("input", nargs="?", help="Path to JSON file. Omit to read from stdin.")
    parser.add_argument("-o", "--output", help="Path to write TSV. Omit to write to stdout.")
    args = parser.parse_args()

    data = load_json(args.input)
    records = list(normalise_records(data))
    records_sorted = sorted(records, key=lambda r: r.get("organism_name", ""))
    rows = [row_from_record(r) for r in records_sorted]
    write_tsv(rows, args.output)

if __name__ == "__main__":
    main()
