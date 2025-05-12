#!/bin/bash
#script that adds different variables to files shared by the pipeline

# Check for the correct number of arguments
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <barcode>" 
    exit 1
fi

# Extract the arguments
barcode_number="$1"
barcode_dir="barcode${barcode_number}"

# Append percentage reads retained after filtering on length to file
echo "Barcode_${barcode_number}:" "$(cat "${barcode_dir}/${barcode_number}_barcode_percent_retained.txt")" >> "percent_reads_retained_length_filter.txt"

# Check if the number of failed reads file is empty and append appropriate message
if [ ! -s "${barcode_dir}/${barcode_number}_num_fail.txt" ]; then
    echo "Barcode_${barcode_number}: Fail reads not analysed" >> "no_fail_reads.txt"
else
    echo "Barcode_${barcode_number}:" "$(cat "${barcode_dir}/${barcode_number}_num_fail.txt")" >> "no_fail_reads.txt"
fi

# Append number of ignored reads after paf parse to file
echo "Barcode_${barcode_number}:" "$(cat "${barcode_dir}/${barcode_number}_number_ignored_reads.tsv")" >> "no_reads_ignored_parse_filter.txt"

# Append the taxaID_count file to all taxaID_count file from paf parse
echo "$(cat "${barcode_dir}/${barcode_number}_taxaID_counts.tsv")" >> "all_taxaID_count.tsv"

# Add a barcode column to the lcaparse_summary file
awk 'BEGIN {OFS="\t"} {print "'${barcode_number}'", $0}' "${barcode_dir}/${barcode_number}_lcaparse_summary.txt" >> "lcaparse_summary.txt"

# Add a barcode column to the lcaparse_perread file
awk 'BEGIN {OFS="\t"} {print "'${barcode_number}'", $0}' "${barcode_dir}/${barcode_number}_lcaparse_perread.txt" >> "lcaparse_perread.txt"