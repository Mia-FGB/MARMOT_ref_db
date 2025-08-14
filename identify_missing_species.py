# Script to compare the species present in DDMMYY_unique_species_python.csv
# and the generated database_info.tsv file from the finished database (created with table_from_json.py)

# Usage: python scripts/identify_missing_species.py <csv_file> <species_col> <json_file> <output_file>
# python scripts/identify_missing_species.py Pathogen_Database_042025_v2/phibase_042025.csv Pathogen_species Pathogen_Database_042025_v2/042025_download_input.json Pathogen_Database_042025_v2/missing_species.csv

import csv
import sys

def read_species_csv(csv_path, species_col):
	species = set()
	with open(csv_path, newline='', encoding='utf-8') as f:
		reader = csv.DictReader(f)
		for row in reader:
			name = row.get(species_col, '').strip()
			if name:
				species.add(name)
	return species

import json

def read_species_json(json_path):
	species = set()
	with open(json_path, encoding='utf-8') as f:
		data = json.load(f)
		if isinstance(data, dict):
			name = data.get("organism_name", "").strip()
			if name:
				species.add(name)
		elif isinstance(data, list):
			for rec in data:
				name = rec.get("organism_name", "").strip()
				if name:
					species.add(name)
	return species

def main():

	# Hardcoded column name and index

	if len(sys.argv) != 5:
		print("Usage: python identify_missing_species.py <csv_file> <species_col> <json_file> <output_file>")
		sys.exit(1)

	csv_file = sys.argv[1]
	csv_species_col = sys.argv[2]
	json_file = sys.argv[3]
	output_file = sys.argv[4]

	csv_species = read_species_csv(csv_file, csv_species_col)
	json_species = read_species_json(json_file)


	only_in_csv = csv_species - json_species

	# Write missing species to CSV, preserving all columns from the original CSV
	with open(csv_file, newline='', encoding='utf-8') as f_in, open(output_file, 'w', newline='', encoding='utf-8') as f_out:
		reader = csv.DictReader(f_in)
		writer = csv.DictWriter(f_out, fieldnames=reader.fieldnames)
		writer.writeheader()
		count = 0
		for row in reader:
			name = row.get(csv_species_col, '').strip()
			if name in only_in_csv:
				writer.writerow(row)
				count += 1
	print(f"Wrote {count} missing species to {output_file}")

if __name__ == "__main__":
	main()


