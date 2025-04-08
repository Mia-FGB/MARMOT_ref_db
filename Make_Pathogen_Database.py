import pandas as pd 
import json
import csv
import argparse
from ete3 import NCBITaxa
import requests

# Script to create a database of pathogens from PHIbase and DEFRA Risk Register
# In conda environment pathogen-database
# Usage: python scripts/Make_Pathogen_Database.py \
#  --phibase Pathogen_Database_Test/phibase_test.csv \
#  --risk_register Pathogen_Database_Test/Risk_Register_Test.csv \
#  --output Download_MMYY_

# Output of this script can then be used in download.py

# Initialize NCBI Taxa object ===
# this is instead of accessionTaxa.sql in R script
ncbi = NCBITaxa()
# Uncomment the following line only if you need to refresh the database (Downloads the taxonomy database)
#ncbi.update_taxonomy_database()

#Functions ===

def download_file(url, filename):
    response = requests.get(url)
    response.raise_for_status()  # Check that the request was successful
    with open(filename, 'wb') as file:
        file.write(response.content)

# Function to get TaxID for a species 
def get_taxid(species_name):
    try:
        # Use NCBITaxa to fetch TaxID
        taxid = ncbi.get_name_translator([species_name])
        print("Getting TaxaID for", species_name)
        return taxid[species_name][0] if species_name in taxid else None
    except Exception as e:
        print(f"Error fetching TaxID for {species_name}: {e}")
        return None

#Function to get the best genome for a given taxid 
def get_best_genome_for_taxid(taxid, ref_gen_df):
    priority_map = {     # Define priority mapping
        'Complete Genome': 1,
        'Chromosome': 2,
        'Scaffold': 3,
        'Contig': 4
    }

    # Filter for the given taxid
    genomes = ref_gen_df[ref_gen_df['taxid'] == taxid].copy()  # Create a copy here
    if genomes.empty:
        return "No genome found for taxid", taxid

    # Assign priorities
    genomes['priority'] = genomes['assembly_level'].map(priority_map).fillna(5)

    # Step 1: Check for a reference genome
    ref_genome = genomes[genomes['refseq_category'] == 'reference genome']
    if not ref_genome.empty:
        # Directly return the single reference genome
        ref_genome = ref_genome.iloc[0]
        return (
            int(ref_genome['taxid']),
            ref_genome['organism_name'],
            ref_genome['assembly_accession'],
            ref_genome['ftp_path'],
            ref_genome['assembly_level']
        )
    
    # If no reference genomes then sort by priority and filter the top priority group
    genomes = genomes.sort_values(by='priority')
    top_priority = genomes.iloc[0]['priority']  # Get the top priority value (e.g. complete genome / chromosme...)
    top_priority_genomes = genomes[genomes['priority'] == top_priority]  # Filter only rows with top priority

    # Select the best genome (longest if ties exist)
    if len(top_priority_genomes) > 1:
        best_genome = get_longest_accession(top_priority_genomes)
    else:
        best_genome = top_priority_genomes.iloc[0]

    return (
        int(best_genome['taxid']),
        best_genome['organism_name'],
        best_genome['assembly_accession'],
        best_genome['ftp_path'],
        best_genome['assembly_level']
    )


def download_and_parse_assembly_stats(ftp_path):
    """
    Downloads the assembly stats file and parses the total size.
    """
    # Construct the URL for the stats file
    assembly_name = ftp_path.strip('/').split('/')[-1]
    stats_url = f"{ftp_path}/{assembly_name}_assembly_stats.txt"
    
    try:
        # Download the file
        response = requests.get(stats_url, stream=True)
        response.raise_for_status()
        
        # Read and parse the lines
        lines = response.text.splitlines()
        for line in lines:
            if "all\tall\tall\tall\ttotal-length" in line:
                # Extract the size (last field in the tab-separated line)
                size = int(line.split('\t')[-1])
                return size
    except Exception as e:
        print(f"Failed to process {stats_url}: {e}")
        return None

def get_longest_accession(df):
    """
    Finds the row for the longest genome for each taxon without modifying the dataframe structure.
    """
    results = []
    print("finding the longest genome for provided taxaID")
    
    for _, row in df.iterrows():
        ftp_path = row['ftp_path']
        #for debugging adding these
        refseq_category = row['refseq_category']
        assembly_level = row['assembly_level']
        
        size = download_and_parse_assembly_stats(ftp_path)
        if size is not None:
            results.append((row, size))  # Keep the original row and the computed size
            #for debugging
            print(f"Processed {ftp_path}: size={size}, refseq_category={refseq_category}, assembly_level={assembly_level}")
    
    # Sort rows by size and release date, and pick the largest and newest for each taxon
    sorted_rows = sorted(results, key=lambda x: (x[1], x[0]['seq_rel_date']), reverse=True)  # Sort by size and date descending
    longest_row = sorted_rows[0][0]  # Get the row with the largest size and newest date

    return longest_row

#Add strings to generate URLs for downloading the genome and MD5 checksum files 
def generate_download_links(accessions_df):
    # Create the 'dlLink' and 'dlLinkMD5' columns
    accessions_df['dlLink'] = accessions_df['ftp_path'].apply(lambda x: f"{x}/{x.split('/')[-1]}_genomic.fna.gz")
    accessions_df['dlLinkMD5'] = accessions_df['ftp_path'].apply(lambda x: f"{x.replace('ftp:', '')}/md5checksums.txt")
    
    # Select the relevant columns
    selected_columns = accessions_df[['taxid', 'organism_name', 'assembly_accession', 'dlLink', 'dlLinkMD5',]]

    return selected_columns

#Create output download file
def save_to_json(df, output_path):
    # Convert the DataFrame to a list of dictionaries
    records = df.to_dict(orient='records')
    
    # Save to a JSON file
    with open(output_path, 'w') as json_file:
        json.dump(records, json_file, indent=4)

def main():
    #Add command line arguments ===
    parser = argparse.ArgumentParser(description="Process PHIbase and Risk Register input files.")
    parser.add_argument("-p", "--phibase", required=True, help="Path to the PHIbase input CSV file")
    parser.add_argument("-r", "--risk_register", required=True, help="Path to the Risk Register input CSV file")
    parser.add_argument("-o", "--output", default="download", help="Output file prefix for the download JSON")

    # Parse the command line arguments
    args = parser.parse_args() 


    # DEFRA Risk Register database ===
    # https://planthealthportal.defra.gov.uk/pests-and-diseases/uk-plant-health-risk-register/downloadEntireRiskRegister.cfm 
    #Download the most up to date version each time
    #read in risk register
    print("Reading in", args.risk_register)
    risk_register = pd.read_csv(args.risk_register)
    risk_register['Type of pest'] = risk_register['Type of pest'].fillna("")  # Fill NaN values with empty strings
    # Define the list of items to remove
    remove = ["Insect", "Mite", "Nematode", "Plant"] 
    # Filter rows where 'Type of pest' is not in the remove list
    risk_register = risk_register[~risk_register['Type of pest'].isin(remove)] 
    # Remove single quotes from 'Pest Name' column if present
    risk_register['Pest Name'] = risk_register['Pest Name'].str.replace("'", "")

    # PHIbase https://github.com/PHI-base/data/tree/master/releases ===
    #Download the most up to date version each time 
    #read in phibase 
    print("Reading in", args.phibase)
    phibase = pd.read_csv(args.phibase)

    # Check for the correct host column (handle potential typo)
    if "Host_description" in phibase.columns:
        host_column = "Host_description"
    elif "Host_descripton" in phibase.columns:  # Note typo here
        host_column = "Host_descripton"
    else:
        raise ValueError("Neither 'Host_description' nor 'Host_descripton' column found in PHIbase data.")

    # Define the list of categories to keep
    keep = ["monocots", "eudicots", "flowering plants",
            "basidiomycetes", "seed plants", "eukaryotes"]

    # Filter the dataframe based on the host column values
    phibase = phibase[phibase[host_column].isin(keep)]

    #Merge the two datasets into one, only keeping species name ===
    # Rename the columns to a common name ('species_name')
    risk_register = risk_register[['Pest Name']].rename(columns={'Pest Name': 'species_name'})
    phibase = phibase[['Pathogen_species']].rename(columns={'Pathogen_species': 'species_name'})

    # Combine the species names from both DataFrames only retain unique values
    unique_species_df = pd.concat([phibase, risk_register]).drop_duplicates().reset_index(drop=True)
    #Note - could add extra species here not in the databases if wanted
    print("Number of unique before taxID species:", len(unique_species_df))
    
    # Get the TaxID for each species in the dataframe ===
    unique_species_df['taxid'] = unique_species_df['species_name'].apply(get_taxid)
    print("Finished getting TaxaIDs")
    # Drop rows with missing TaxID
    unique_species_df = unique_species_df.dropna(subset=['taxid'])
    #Convert taxid to integer
    unique_species_df['taxid'] = unique_species_df['taxid'].astype(int)

    unique_species_df.to_csv("unique_species_python.csv", index=False)

    # Download Refseq & Genbank tables ===
    # Only need to do this occassionally to keep up to date as it takes a long time
    print("Downloading RefSeq dataframe")
    download_file("https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt", "assembly_summary_refseq.txt")
    print("Downloading GenBank dataframe")
    download_file("https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt", "assembly_summary_genbank.txt")

    # Read in RefSeq dataframe
    print("Reading in RefSeq dataframe")
    refseq = pd.read_csv("assembly_summary_refseq.txt", sep='\t', skiprows=1, header=0, dtype='object')
    # Remove columns with NaN names
    refseq = refseq.loc[:, refseq.columns.notna()]
    refseq = refseq.rename(columns={'#assembly_accession': 'assembly_accession'})

    # Read in GenBank dataframe
    print("Reading in GenBank dataframe")
    genbank = pd.read_csv("assembly_summary_genbank.txt", sep='\t', skiprows=1, header=0, dtype='object')
    # Remove columns with NaN names
    genbank = genbank.loc[:, genbank.columns.notna()]
    genbank = genbank.rename(columns={'#assembly_accession': 'assembly_accession'})

    #Merge the assembly databases ====
    #First Remove entries from genbank that are present in refseq to avoid redundancies 
    # Convert refseq's 'gbrs_paired_asm' column to a set for faster comparison
    refseq_set = set(refseq['gbrs_paired_asm'])
    # Filter genbank DataFrame by excluding entries where 'assembly_accession' is in refseq_set
    genbank_filtered = genbank[~genbank['assembly_accession'].isin(refseq_set)]
    # Merge the two DataFrames (combine them into one)
    ref_gen = pd.concat([refseq, genbank_filtered], ignore_index=True)
    ref_gen['taxid'] = ref_gen['taxid'].astype(int)


    # Filter ref_gen for species present among the PHIbase and DEFRA pathogens - using taxid
    ref_gen_match = ref_gen[ref_gen['taxid'].isin(unique_species_df['taxid'])]

    #Only take rows with an https entry
    ref_gen_match = ref_gen_match[ref_gen_match['ftp_path'].str.contains("https://")]

    # Create a list to store the data
    accessions_list = []
    # Iterate over the unique taxids in ref_gen_match and call the function
    for taxid in sorted(ref_gen_match['taxid'].unique()):
        print(f"Processing taxid {taxid}")
        try:
            # Call your function to get the best genome
            taxid, organism_name, assembly_accession, ftp_path, genome_type = get_best_genome_for_taxid(taxid, ref_gen_match)
            
            # Append the result as a dictionary to the list
            accessions_list.append({
                'taxid': taxid,
                'organism_name': organism_name,
                'assembly_accession': assembly_accession,
                'ftp_path': ftp_path,
                'type': genome_type
            })
        except Exception as e:
            print(f"Error processing taxid {taxid}: {e}")

    # Convert the list of dictionaries to a DataFrame
    accessions_df = pd.DataFrame(accessions_list)

    accessions_df_links = generate_download_links(accessions_df)

    # Define the output path for the JSON file
    output_path = args.output + ".json"

    # Save the results to JSON
    save_to_json(accessions_df_links, output_path)

if __name__ == "__main__":
    main()
