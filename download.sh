#!/bin/bash
#Note this script may need further imprivements if there are files inside download that aren't in download.txt
#no longer use this script - instead use download.py which is much fater

mkdir -p download
cd download

# Check if the correct number of arguments is provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 MMYYYY"
    exit 1
fi

# Extract the date argument
DATE=$1

# Construct the file paths using the provided date
FILE_PATH="../pathogen_database_${DATE}.fa"
ERROR_LOG="../error_log_${DATE}.txt"
OUTPUT_LOG="../output_log_${DATE}.txt"

# Print the constructed file paths (for debugging purposes)
echo "Using file path: $FILE_PATH" >> $OUTPUT_LOG
echo "Using error log: $ERROR_LOG"  >> $OUTPUT_LOG
echo "Using output log: $OUTPUT_LOG"  >> $OUTPUT_LOG
echo "in this directory $(pwd)"  >> $OUTPUT_LOG

# Function to extract accession number from URL
get_accession_number() {
    url="$1"
    accession=$(echo "$url" | awk -F'/' '{print $(NF-1)}')
    echo "$accession"
}

# Function to calculate checksum for a file
calculate_checksum() {
    file="$1"
    md5sum "$file" | awk '{print $1}'
}

# Function to download file
download_file() {
    url="$1"
    wget -c "$url"
}

# Function to extract checksum for a specific file from MD5 URL content
extract_checksum() {
    md5_url="$1"
    filename="$2"
    checksum_line=$(wget -qO- "$md5_url" | grep "$filename")
    checksum=$(echo "$checksum_line" | awk '{print $1}')
    echo "$checksum" 
}
downloaded_files=()

# Download files that are not already present
echo "Downloading missing files..."  >> $OUTPUT_LOG
while IFS=$'\t' read -r organism_id taxid version md5_url fasta_url; do
    accession=$(get_accession_number "$fasta_url")
    filename="${accession}_genomic.fna.gz"

    if [ ! -f "$filename" ]; then
        echo "Downloading $filename..." >> $OUTPUT_LOG
        download_file "$fasta_url"
        downloaded_files+=("$filename")  # Add the downloaded filename to the list
    else
        echo "$filename already exists, skipping download."  >> $OUTPUT_LOG
    fi
done < ../download.txt

# Run checksums for downloaded files and attempt to redownload corrupted files

echo "Calculating checksums of new files..."  >> $OUTPUT_LOG
while IFS=$'\t' read -r organism_id taxid version md5_url fasta_url; do
    accession=$(get_accession_number "$fasta_url")
    filename="${accession}_genomic.fna.gz"

    if [[ " ${downloaded_files[@]} " =~ " ${filename} " ]]; then
        echo "Calculating checksum for $filename..."  >> $OUTPUT_LOG
        expected_checksum=$(extract_checksum "$md5_url" "$filename")
        calculated_checksum=$(calculate_checksum "$filename")
        
        if [ "$expected_checksum" != "$calculated_checksum" ]; then
            echo "Error: Checksum mismatch for $filename. Downloaded file may be corrupted. Attempting to redownload..."  >> $OUTPUT_LOG
            download_file "$fasta_url"
            new_calculated_checksum=$(calculate_checksum "$filename")
            if [ "$expected_checksum" != "$new_calculated_checksum" ]; then
                echo "Failed to download $filename correctly. Adding to error log: $ERROR_LOG"  >> $OUTPUT_LOG
                echo "Organism ID: $organism_id, TaxID: $taxid, Version: $version, MD5 URL: $md5_url, FASTA URL: $fasta_url" >> "$ERROR_LOG"
            fi
        fi
    fi
done < ../download.txt

# Check taxids before modifying the fasta headers
echo "Checking taxids..."  >> $OUTPUT_LOG
if [ -f taxids.txt ]; then rm taxids.txt; fi
while IFS=$'\t' read -r organism_id taxid version md5_url fasta_url; do
    accession=$(get_accession_number "$fasta_url")
    filename="${accession}_genomic.fna.gz"
    echo "$filename $taxid" >> taxids.txt
done < ../download.txt

# Modify fasta headers by assigning the taxid number in front of the description and concatenate to a single fasta file
echo "Modifying fasta headers..."  >> $OUTPUT_LOG
if [ -f ../$FILE_PATH ]; then rm ../$FILE_PATH; fi
while IFS=$'\t' read -r organism_id taxid version md5_url fasta_url; do
    accession=$(get_accession_number "$fasta_url")
    filename="${accession}_genomic.fna.gz"
    gunzip -c "$filename" | sed "s/^>/>taxid|${taxid}|/g" >> ../$FILE_PATH
done < ../download.txt

exit
