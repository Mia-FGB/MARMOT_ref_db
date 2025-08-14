import pandas as pd
import json
import csv
import argparse
from ete3 import NCBITaxa
import requests

# Example:
# python Make_Pathogen_Database_one_per_species.py \
#   --phibase "$PHIBASE_CSV" \
#   --risk_register "$RISK_REGISTER_CSV" \
#   --output "$OUTDIR/${DATE_TAG}_"

REQUEST_TIMEOUT = 30  # seconds

# Optional synonym map (left empty so nothing is forced)
KNOWN_SYNONYMS = {
    # "Some tricky name": 12345,
}

# Initialise NCBI Taxa
ncbi = NCBITaxa()
# If you need to refresh the local taxonomy, run once manually:
# ncbi.update_taxonomy_database()

# --------------------------
# Helpers
# --------------------------

def normalise_name(name: str) -> str:
    s = str(name)
    s = s.replace("’", "'").replace("‘", "'").replace("`", "'")
    s = s.replace("–", "-").replace("—", "-")
    s = " ".join(s.split())
    return s.strip().strip('"').strip("'")


def download_file(url, filename):
    resp = requests.get(url, timeout=REQUEST_TIMEOUT)
    resp.raise_for_status()
    with open(filename, 'wb') as fh:
        fh.write(resp.content)


def get_taxid(species_name):
    """Resolve a species to a TaxID with fallbacks and clear logging, without forcing any species."""
    try:
        name_raw = str(species_name)
        name = normalise_name(name_raw)

        # exact
        tx = ncbi.get_name_translator([name])
        if name in tx and tx[name]:
            tid = int(tx[name][0])
            print(f"Getting TaxaID for {name} -> {tid}")
            return tid

        # variants
        variants = {name, name.title(), name.lower(), name.upper()}
        for v in variants:
            tx = ncbi.get_name_translator([v])
            if v in tx and tx[v]:
                tid = int(tx[v][0])
                print(f"Getting TaxaID for {name} -> {tid} (via {v})")
                return tid

        # known synonyms or pins if you ever choose to add them
        if name in KNOWN_SYNONYMS:
            tid = int(KNOWN_SYNONYMS[name])
            print(f"Getting TaxaID for {name} -> {tid} (via KNOWN_SYNONYM)")
            return tid

        print(f"Getting TaxaID for {name} -> FAILED")
        return None
    except Exception as e:
        print(f"Error fetching TaxID for {species_name}: {e}")
        return None


def expand_to_descendant_taxa(taxid):
    """Include all descendants such as subspecies, forma specialis and strains."""
    try:
        descendants = ncbi.get_descendant_taxa(int(taxid), intermediate_nodes=True) or []
        expanded = {int(taxid)}
        expanded.update(int(t) for t in descendants)
        return expanded
    except Exception as e:
        print(f"Could not expand descendants for {taxid}: {e}")
        return {int(taxid)}


def to_https(url: str) -> str:
    url = str(url)
    if url.startswith("ftp://"):
        return "https://" + url[len("ftp://"):]
    return url


def download_and_parse_assembly_stats(ftp_path):
    """Download the assembly stats file and return total length."""
    base = str(ftp_path).rstrip("/")
    base = to_https(base)
    asm = base.split("/")[-1]
    stats_url = f"{base}/{asm}_assembly_stats.txt"
    try:
        response = requests.get(stats_url, stream=True, timeout=REQUEST_TIMEOUT)
        response.raise_for_status()
        for line in response.text.splitlines():
            if "all\tall\tall\tall\ttotal-length" in line:
                size = int(line.split('\t')[-1])
                return size
    except Exception as e:
        print(f"Failed to process {stats_url}: {e}")
    return None


def get_longest_accession(df):
    """Pick the longest assembly by stats, tie break by seq_rel_date."""
    results = []
    print("finding the longest genome for provided taxaID")
    for _, row in df.iterrows():
        ftp_path = row['ftp_path']
        refseq_category = row.get('refseq_category', 'na')
        assembly_level = row.get('assembly_level', 'na')
        size = download_and_parse_assembly_stats(ftp_path)
        if size is not None:
            results.append((row, size))
            print(f"Processed {ftp_path}: size={size}, refseq_category={refseq_category}, assembly_level={assembly_level}")
        else:
            print(f"Could not read assembly stats for {ftp_path}")

    if not results:
        print("No assembly stats available, falling back to first candidate")
        return df.iloc[0]

    def sort_key(item):
        row, size = item
        date = row.get('seq_rel_date', '')
        return (size, date)

    sorted_rows = sorted(results, key=sort_key, reverse=True)
    return sorted_rows[0][0]


def select_best_assembly(candidates_df):
    """Select best assembly using ref or representative preference, assembly level, then length and date."""
    priority_map = {
        'Complete Genome': 1,
        'Chromosome': 2,
        'Scaffold': 3,
        'Contig': 4
    }
    if candidates_df.empty:
        return None

    df = candidates_df.copy()
    df['priority'] = df['assembly_level'].map(priority_map).fillna(5)

    # Prefer reference or representative
    ref_like = df[df['refseq_category'].isin(['reference genome', 'representative genome'])]
    if not ref_like.empty:
        df = ref_like

    # Sort by priority and restrict to top tier
    df = df.sort_values(by='priority', kind='stable')
    top_priority = df.iloc[0]['priority']
    df_top = df[df['priority'] == top_priority]

    # If multiple remain, use assembly stats to pick the longest and newest
    if len(df_top) > 1:
        return get_longest_accession(df_top)
    return df_top.iloc[0]


def generate_download_links_row(ftp_path):
    base = to_https(str(ftp_path).rstrip("/"))
    asm = base.split("/")[-1]
    genomic_url = f"{base}/{asm}_genomic.fna.gz"
    md5_url = f"{base}/md5checksums.txt"
    return genomic_url, md5_url


def save_to_json(df, output_path):
    records = df.to_dict(orient='records')
    with open(output_path, 'w') as json_file:
        json.dump(records, json_file, indent=4)


def save_summary_and_missing(accessions_df, output_path_prefix, missing_species_list):
    summary = {
        "total_species_with_selection": len(accessions_df),
        "assembly_level_counts": accessions_df['type'].value_counts(dropna=False).to_dict(),
        "refseq_category_counts": accessions_df['source_db'].value_counts(dropna=False).to_dict(),
        "missing_species_count": len(missing_species_list),
    }
    summary_path = output_path_prefix + "_summary.json"
    with open(summary_path, 'w') as summary_file:
        json.dump(summary, summary_file, indent=4)
    print(f"Summary file saved to {summary_path}")

    missing_path = output_path_prefix + "missing_species.json"
    with open(missing_path, 'w') as fh:
        json.dump(missing_species_list, fh, indent=4)
    print(f"Missing species list saved to {missing_path}")

# --------------------------
# Main
# --------------------------

def main():
    parser = argparse.ArgumentParser(description="Process PHIbase and Risk Register input files, selecting one genome per species.")
    parser.add_argument("-p", "--phibase", required=True, help="Path to the PHIbase input CSV file")
    parser.add_argument("-r", "--risk_register", required=True, help="Path to the Risk Register input CSV file")
    parser.add_argument("-o", "--output", default="download", help="Output file prefix for the download JSON")
    args = parser.parse_args()

    # Read sources
    print("Reading in", args.risk_register)
    risk_register = pd.read_csv(args.risk_register)
    risk_register['Type of pest'] = risk_register['Type of pest'].fillna("")
    remove = ["Insect", "Mite", "Nematode", "Plant"]
    risk_register = risk_register[~risk_register['Type of pest'].isin(remove)]
    risk_register['Pest Name'] = risk_register['Pest Name'].astype(str).str.replace("'", "", regex=False).str.strip()

    print("Reading in", args.phibase)
    phibase = pd.read_csv(args.phibase)
    phibase['Pathogen_species'] = phibase['Pathogen_species'].astype(str).str.strip()

    # Build species list
    risk_register = risk_register[['Pest Name']].rename(columns={'Pest Name': 'species_name'})
    phibase = phibase[['Pathogen_species']].rename(columns={'Pathogen_species': 'species_name'})
    species_df = pd.concat([phibase, risk_register]).drop_duplicates().reset_index(drop=True)
    species_df['species_name'] = species_df['species_name'].astype(str).str.strip()

    print("Number of unique before taxaID species:", len(species_df))

    # Resolve taxids for species
    species_df['species_taxid'] = species_df['species_name'].apply(get_taxid)
    print("Finished getting TaxaIDs")

    # Report failures
    failed = species_df[species_df['species_taxid'].isna()]
    if not failed.empty:
        print("Names with no TaxID:", list(failed['species_name']))

    # Keep only resolved
    species_df = species_df.dropna(subset=['species_taxid']).copy()
    species_df['species_taxid'] = species_df['species_taxid'].astype(int)

    # Save the resolved list
    species_df.to_csv(args.output + "unique_species_python.csv", index=False)

    # Download assembly tables - only every so often 
    # print("Downloading RefSeq dataframe")
    # download_file("https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt", "assembly_summary_refseq.txt")
    # print("Downloading GenBank dataframe")
    # download_file("https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt", "assembly_summary_genbank.txt")

    # Read them
    print("Reading in RefSeq dataframe")
    refseq = pd.read_csv("assembly_summary_refseq.txt", sep='\t', skiprows=1, header=0, dtype='object', low_memory=False)
    refseq = refseq.loc[:, refseq.columns.notna()]
    refseq = refseq.rename(columns={'#assembly_accession': 'assembly_accession'})

    print("Reading in GenBank dataframe")
    genbank = pd.read_csv("assembly_summary_genbank.txt", sep='\t', skiprows=1, header=0, dtype='object', low_memory=False)
    genbank = genbank.loc[:, genbank.columns.notna()]
    genbank = genbank.rename(columns={'#assembly_accession': 'assembly_accession'})

    # Merge, removing GenBank entries mirrored in RefSeq
    refseq_set = set(refseq['gbrs_paired_asm'])
    genbank_filtered = genbank[~genbank['assembly_accession'].isin(refseq_set)]
    ref_gen = pd.concat([refseq, genbank_filtered], ignore_index=True)

    # Keep any ftp or https entry, we rewrite for downloads later
    ref_gen = ref_gen[ref_gen['ftp_path'].astype(str).str.contains("://", na=False)].copy()

    # Dtypes
    ref_gen['taxid'] = ref_gen['taxid'].astype(int)

    # For each species, expand to descendants and pick one best assembly
    accessions_rows = []
    missing_species = []

    print(f"Total input species with TaxIDs: {len(species_df)}")
    for _, srow in species_df.iterrows():
        species_name = srow['species_name']
        species_taxid = int(srow['species_taxid'])

        # Expand to descendants so subspecies and formae speciales are included
        expanded_taxids = expand_to_descendant_taxa(species_taxid)

        # Candidate assemblies for this species including descendants
        cand = ref_gen[ref_gen['taxid'].isin(expanded_taxids)].copy()
        print(f"[SELECT] {species_name} (taxid {species_taxid}) candidates: {len(cand)}")

        # Pick the single best assembly
        best = select_best_assembly(cand)
        if best is None:
            # record as missing
            missing_species.append({
                "species_name": species_name,
                "species_taxid": species_taxid
            })
            print(f"[MISS] No assembly selected for {species_name}")
            continue

        # Build download links and source db flag
        dl, md5 = generate_download_links_row(best['ftp_path'])
        source_db = "RefSeq" if "refseq" in str(best['ftp_path']).lower() else "GenBank"

        accessions_rows.append({
            "species_name": species_name,                    # input species
            "species_taxid": species_taxid,                  # input species taxid
            "selected_taxid": int(best['taxid']),            # taxid of chosen assembly (may be descendant)
            "organism_name": best['organism_name'],          # organism label from assembly table
            "assembly_accession": best['assembly_accession'],
            "ftp_path": best['ftp_path'],
            "type": best['assembly_level'],
            "source_db": source_db,
            "dlLink": dl,
            "dlLinkMD5": md5
        })

    # One row per species by construction
    accessions_df = pd.DataFrame(accessions_rows)

    # Write the download list expected by downstream tooling
    dl_cols = ["species_name", "species_taxid", "selected_taxid", "organism_name",
               "assembly_accession", "dlLink", "dlLinkMD5", "type", "source_db"]
    output_path = args.output + "download_input.json"
    save_to_json(accessions_df[dl_cols], output_path)
    print(f"Wrote one-genome-per-species download list to {output_path}")

    # Save summary and missing species list
    save_summary_and_missing(accessions_df, args.output, missing_species)

    # Friendly tail line
    if missing_species:
        names = ", ".join(ms['species_name'] for ms in missing_species[:10])
        extra = " …" if len(missing_species) > 10 else ""
        print(f"Species with no identified URLs ({len(missing_species)}): {names}{extra}")
    else:
        print("All species had a selected genome")

if __name__ == "__main__":
    main()
