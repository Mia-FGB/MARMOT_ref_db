#!/bin/bash

# Check for the correct number of arguments
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <barcode_number> <reference_database> <filter_length>"
    exit 1
fi

# Extract the arguments
barcode_number="$1"
reference_database="$2"
filter_length="$3"

barcode_dir="./barcode${barcode_number}"

# Run Minimap2
sbatch --mem 60G -p ei-medium -o "$barcode_dir/${barcode_number}_minimap2.out" --wrap "minimap2 -x map-ont ${reference_database} \"$barcode_dir/${barcode_number}_barcode_${filter_length}bp.fastq\" > \"$barcode_dir/${barcode_number}_mapped.paf\""
echo "Minimap batch ID: $minimap_job_id" >> $barcode_dir/${barcode_number}_minimap_job_id.txt