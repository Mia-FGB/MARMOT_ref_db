import sys
import json
import csv

# get details from the download.json to build a table 
# python scripts/table_from_json.py Pathogen_Database_0825/0825_download_input.json Pathogen_Database_0825/0825_database_info.tsv

def load_json(path):
    with open(path, "r", encoding="utf-8") as f:
        data = json.load(f)
    if isinstance(data, dict):
        return [data]
    return data

def main():
    if len(sys.argv) != 3:
        print("Usage: python table_from_json.py <input.json> <output.tsv>")
        sys.exit(1)

    input_json = sys.argv[1]
    output_tsv = sys.argv[2]

    records = load_json(input_json)
    if not records:
        print("No records found in JSON.")
        sys.exit(1)

    # Only keep selected columns
    headers = ["species_name", "selected_taxid", "assembly_accession", "type"]

    with open(output_tsv, "w", newline="", encoding="utf-8") as f_out:
        writer = csv.DictWriter(f_out, fieldnames=headers, delimiter="\t")
        writer.writeheader()
        for rec in records:
            filtered = {h: rec.get(h, "") for h in headers}
            writer.writerow(filtered)

if __name__ == "__main__":
    main()