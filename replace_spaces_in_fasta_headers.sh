#!/bin/bash

# Script to replace spaces in FASTA headers with underscores.
# To be used on the reference database FASTA files.

# Check arguments
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 input.fasta output.fasta"
    exit 1
fi

input="$1"
output="$2"

# Process the FASTA
awk '{
    if ($0 ~ /^>/) {
        gsub(" ", "_", $0)
        print
    } else {
        print
    }
}' "$input" > "$output"

echo "Done. Modified FASTA saved to: $output"
