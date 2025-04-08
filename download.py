#!/usr/bin/env python3

#run in Conda environment = pathogen_database
# Uses the output from Make_Pathogen_Database.py
# This script downloads FASTA files from NCBI, checks their checksums, modifies headers, and concatenates them into a single database.
# python scripts/download.py --input Download_MMYY_ --date MMYYYY


import os
import hashlib
import logging
import requests
from urllib.parse import urlparse
import gzip
import argparse
import json
from tqdm import tqdm
import shutil

def download_file(fasta_url, filename):
    logging.info("Downloading %s...", filename)
    response = requests.get(fasta_url, stream=True)
    with open(filename, 'wb') as file:
        for chunk in response.iter_content(chunk_size=8192): #8 bytes
            file.write(chunk)

# Function to calculate checksum for the newly downloaded file
def calculate_checksum(file_path):
    sha256_hash = hashlib.md5() #initialize the hash object
    with open(file_path, 'rb') as f:
        for byte_block in iter(lambda: f.read(4096), b""): #4kb
            sha256_hash.update(byte_block)
    return sha256_hash.hexdigest()

# Function to extract checksum for a specific file from MD5 URL content
def extract_checksum(md5_url, filename):
    response = requests.get(md5_url)
    for line in response.text.splitlines():
        if filename in line:
            return line.split()[0]
    return None

# Function to immediately check checksum after download and log errors
def check_and_log_checksum(error_log, filename, md5_url, organism_id, taxid,fasta_url):
    expected_checksum = extract_checksum(md5_url, filename)
    if expected_checksum:
        calculated_checksum = calculate_checksum(filename)
        if expected_checksum != calculated_checksum:
            logging.error(f"Checksum mismatch for {filename}. Redownloading...")
            download_file(fasta_url, filename)
            new_calculated_checksum = calculate_checksum(filename)
            if expected_checksum != new_calculated_checksum:
                logging.error(f"Failed to download {filename} correctly.")
                with open(error_log, 'a') as err_log:
                    err_log.write(f"Organism ID: {organism_id}, TaxID: {taxid}, MD5 URL: {md5_url}, FASTA URL: {fasta_url}\n")
        else:
            logging.info(f"Checksum for {filename} matches.")

# Function to modify fasta headers
def modify_fasta_headers(filename, taxid, organism_name):
    # Read the content of the file
    with gzip.open(filename, 'rt') as f:
        lines = f.readlines()
    logging.info(f"Modifying headers for {filename}...")

    # Modify the headers
    modified_lines = []
    for line in lines:
        if line.startswith(">"):
            modified_lines.append(f">taxid|{taxid}|{organism_name}|{line[1:]}")
        else:
            modified_lines.append(line)

    # Write the modified content back to the file
    with gzip.open(filename, 'wt') as f:
        f.writelines(modified_lines)

# Function to unzip a file
def unzip_file(input_filename, output_filename):
    with gzip.open(input_filename, 'rb') as infile:
        with open(output_filename, 'wb') as outfile:
            shutil.copyfileobj(infile, outfile)

# Function to concatenate files into one large database
def concatenate_files(file_list, output_filename):
    with open(output_filename, 'wb') as outfile:
        for filename in file_list:
            temp_unzipped_filename = filename.replace('.gz', '')
            unzip_file(filename, temp_unzipped_filename)
            with open(temp_unzipped_filename, 'rb') as infile:
                shutil.copyfileobj(infile, outfile)
            os.remove(temp_unzipped_filename)  # Remove the temporary unzipped file


def main():
    # Ensure the download directory exists
    os.makedirs("download", exist_ok=True)
    os.chdir("download")

    #Add command line arguments ===
    parser = argparse.ArgumentParser(description="Provide Date as a prefix and input file")
    parser.add_argument("-i", "--input", required=True, help="Input file with list of URLs")
    parser.add_argument("-d", "--date", required=True, help="Date in MMYYYY format")

    # Parse the command line arguments
    args = parser.parse_args() 

    # Set up logging
    output_log = '../logs/' + args.date + '_output_log.txt'
    logging.basicConfig(filename=output_log, level=logging.INFO, format='%(asctime)s - %(message)s')
    error_log = '../logs/' + args.date + '_error_log.txt'

    all_files = []

#    Open the JSON file and process it
    with open('../' + args.input + '.json', 'r') as f:
        data = json.load(f)
        for entry in tqdm(data, desc="Processing entries"):
            taxid = entry['taxid']
            organism_name = entry['organism_name']
            assembly_accession = entry['assembly_accession']
            fasta_url = entry['dlLink']
            md5_url = entry['dlLinkMD5']
            filename = fasta_url.strip('/').split('/')[-2] + "_genomic.fna.gz"

            all_files.append(filename)

#           Check if the file already exists in download directory
            if not os.path.isfile(filename):
                download_file(fasta_url, filename)
                check_and_log_checksum(error_log, filename, md5_url, organism_name, taxid, fasta_url)
                modify_fasta_headers(filename, taxid, organism_name)
            else:
                logging.info(f"{filename} already exists, skipping download.")
        
    # Concatenate the downloaded files into one large database - save where the script is run
    output_filename = f"pathogen_database_{args.date}.fa"
    concatenate_files(all_files, output_filename)
    logging.info(f"Concatenated files into {output_filename}")
         

if __name__ == "__main__":
    main()