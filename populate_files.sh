#!/bin/bash
#Script to populate the txt file outputs to combine stats from different barcodes
#Ordinarilly single_barcode_proccess should do this

barcode_list=("17" "18" "20" "21" "22" "25" "26" "27" "28")

for filename in "percent_reads_retained_first_filter.txt" "number_fail_reads.txt" "number_reads_ignored_parse_filter.txt" "all_taxaID_count.tsv"; do
    touch -c "$filename"
done

percent_results_file="percent_reads_retained_length_filter.txt"
fail_number_file="no_fail_reads.txt"
ignored_reads_number_file="no_reads_ignored_parse_filter.txt"
all_taxaID_count="all_taxaID_count.tsv"

for barcode_number in "${barcode_list[@]}"; 
    do 
    
    # Calculate the barcode directory path
    barcode_dir="barcode${barcode_number}"
    
    echo "Barcode_${barcode_number}: $(cat "$barcode_dir/${barcode_number}_barcode_percent_retained.txt")" >> "$percent_results_file"
    echo "Barcode_${barcode_number}: $(cat "$barcode_dir/${barcode_number}_num_fail.txt")" >> "$fail_number_file"
    echo "barcode${barcode_number}: $(cat "$barcode_dir/${barcode_number}_number_ignored_reads.tsv")" >> "$ignored_reads_number_file"
    echo "$(cat "$barcode_dir/${barcode_number}_taxaID_counts.tsv")" >> "$all_taxaID_count"

done
