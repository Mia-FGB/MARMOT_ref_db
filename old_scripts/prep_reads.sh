#!/bin/bash

#Script that takes raw reads and outputs concatenated fast reads for the barcode
#Also calculates contig stats for the barcode - output into different folder
#Filters the pass reads based on their length
#Run it within results directory -> ./prep_reads.sh barcode raw_read_location

# Check for the correct number of arguments
if [ "$#" -ne 4 ]; then
    echo "Usage: $0 <barcode> <raw_read_location> <filter_length> <scratch_dir>" 
    exit 1
fi

# Extract the arguments
barcode_number="$1"
location="$2"
filter_length="$3"
scratch_dir="$4"

# Check if the specified location exists
if [ ! -d "$location" ]; then
    echo "Error: Location '$location' does not exist."
    exit 1
fi

# Specify barcode directory
barcode_dir="./barcode${barcode_number}"

# Create the concatenated read output directory if it doesn't exist (-p flag)
#This should be scratch area 
mkdir -p "$scratch_dir"

# Concatenate the pass files based on the provided barcode number and location 
# only unzips when necessary and outputs the new concatenated file into scratch area
if [ -d "$location/fastq_pass/barcode${barcode_number}" ]; then
    for file in "$location/fastq_pass/barcode${barcode_number}"/*; do
        if [[ "$file" == *.fastq || "$file" == *.fq ]]; then
            cat "$file" >> "$scratch_dir/${barcode_number}_barcode.fastq"
        elif [[ "$file" == *.gz ]]; then
            zcat "$file" >> "$scratch_dir/${barcode_number}_barcode.fastq"
        else
            echo "Warning: Skipping unrecognized file format '$file'."
        fi
    done
elif [ -d "$location/fastq/barcode${barcode_number}" ]; then
    for file in "$location/fastq/barcode${barcode_number}"/*; do
        if [[ "$file" == *.fastq || "$file" == *.fq ]]; then
            cat "$file" >> "$scratch_dir/${barcode_number}_barcode.fastq"
        elif [[ "$file" == *.gz ]]; then
            zcat "$file" >> "$scratch_dir/${barcode_number}_barcode.fastq"
        else
            echo "Warning: Skipping unrecognized file format '$file'."
        fi
    done
else
    echo "Error: Neither fastq or fastq_pass directory exists for location ${location} and barcode ${barcode_number}."
    exit 1
fi

# Run the Perl script (get_contig_stats) on the concatenated FASTQ file - not length filtered
# Outputs the contig stats into the barcode directory not scratch area
if ! get_contig_stats.pl -q -i "$scratch_dir/${barcode_number}_barcode.fastq" > "$barcode_dir/${barcode_number}_contig_stats.txt"; then
    echo "Error: Failed to run get_contig_stats.pl on the concatenated FASTQ file."
    exit 1
fi

# Check if the fastq_fail directory exists
if [ -d "$location/fastq_fail/barcode${barcode_number}" ]; then
    # Unzip fail barcodes
    zcat "$location/fastq_fail/barcode${barcode_number}"/* > "$scratch_dir/fail_barcode${barcode_number}.fastq"
    if [ ! -s "$scratch_dir/fail_barcode${barcode_number}.fastq" ]; then
        echo "Warning: No fail barcodes found for barcode ${barcode_number}."
    else
        # Count the number of fail barcodes
        if ! wc -l "$scratch_dir/fail_barcode${barcode_number}.fastq" | tail -1 | awk '{print $1/4}' > "$barcode_dir/${barcode_number}_num_fail.txt"; then
            echo "Error: Failed to count the number of fail barcodes."
        fi

        #Removed the section to delete this file as it will live in scratch
    fi
else
    echo "Warning: fastq_fail directory for barcode ${barcode_number} does not exist. Skipping fail barcode processing."
fi

# Filter pass reads on length and create a new fastq file
if ! awk -f /ei/projects/7/724b6a9a-6416-47eb-be58-d3737023e29b/scratch/getBigReads.awk -v min="${filter_length}" "$scratch_dir/${barcode_number}_barcode.fastq" > "$scratch_dir/${barcode_number}_barcode_${filter_length}bp.fastq"; then
    echo "Error: Failed to filter pass reads based on length."
    exit 1
fi

# Calculate the percentage of reads retained after filtering
if ! echo "scale=2; (100 * $(wc -l < "$scratch_dir/${barcode_number}_barcode_${filter_length}bp.fastq") / $(wc -l < "$scratch_dir/${barcode_number}_barcode.fastq"))" | bc > "$barcode_dir/${barcode_number}_barcode_percent_retained.txt"; then
    echo "Error: Failed to calculate the percentage of reads retained after filtering."
    exit 1
fi
