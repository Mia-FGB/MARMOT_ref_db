#!/usr/bin/env python3

# usage python table_from_json.py <input.json> -o <output.tsv> --genbank-summary assembly_summary_genbank.txt


import sys
import json
import argparse
from pathlib import Path
from typing import Iterable, Dict, Any, List, Optional

# Fields we want in the output
OUT_HEADERS = ["Organism name", "TaxID", "Accession number", "Assembly level"]
def build_genbank_lookup(summary_path: str) -> Dict[str, str]:
    """
    Build a lookup from both GenBank and RefSeq accessions to assembly_level (first word only) from the assembly_summary_genbank.txt file.
    """
    lookup = {}
    with open(summary_path, "r", encoding="utf-8") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 18:
                continue
            assembly_level = fields[11].split()[0] if fields[11] else ""
            genbank_acc = fields[0]
            refseq_acc = fields[17] if len(fields) > 17 else ""
            if genbank_acc:
                lookup[genbank_acc] = assembly_level
            if refseq_acc and refseq_acc != 'na':
                lookup[refseq_acc] = assembly_level
    return lookup

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

def row_from_record(rec: Dict[str, Any], genbank_lookup: Dict[str, str]) -> List[str]:
    """
    Build a TSV row in the output order, including only Assembly level from GenBank summary.
    Missing fields become empty strings.
    """
    out = {h: "" for h in OUT_HEADERS}
    for src_key, dst_header in KEY_MAP.items():
        if src_key in rec and rec[src_key] is not None:
            out[dst_header] = str(rec[src_key])
    acc = rec.get("assembly_accession", "")
    out["Assembly level"] = genbank_lookup.get(acc, "")
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
        description="Convert genome JSON to a TSV with Organism name, TaxID, Accession number, Assembly level, and RefSeq category from GenBank summary."
    )
    parser.add_argument("input", nargs="?", help="Path to JSON file. Omit to read from stdin.")
    parser.add_argument("-o", "--output", help="Path to write TSV. Omit to write to stdout.")
    parser.add_argument("--genbank-summary", default="assembly_summary_genbank.txt", help="Path to assembly_summary_genbank.txt file.")
    args = parser.parse_args()

    genbank_lookup = build_genbank_lookup(args.genbank_summary)
    data = load_json(args.input)
    records = list(normalise_records(data))
    records_sorted = sorted(records, key=lambda r: r.get("organism_name", ""))
    rows = [row_from_record(r, genbank_lookup) for r in records_sorted]
    write_tsv(rows, args.output)

if __name__ == "__main__":
    main()
