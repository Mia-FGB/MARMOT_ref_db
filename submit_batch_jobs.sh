#!/bin/bash
#Older script which has been replaced with one that submits arrays

# Define an array of barcode numbers
barcode_list=($(seq -w 01 87))

# Define raw read location - level above fastqpass
location=("/ei/projects/9/9742f7cc-c169-405d-bf27-cd520e26f0be/data/raw/CF_2023_combined_sup_runs/")

#Define length to filter reads to 
filter_length=("300")

#Define mapping database
reference_database=(
    /ei/projects/9/9742f7cc-c169-405d-bf27-cd520e26f0be/data/results/phibase/pathogen_database_080524.fa
)

# Create txt files if they don't exist
touch ./percent_reads_retained_length_filter.txt
touch ./no_fail_reads.txt
touch ./no_reads_ignored_parse_filter.txt
touch ./all_taxaID_count.tsv

# Submit SLURM jobs in a loop for each barcode
for barcode_number in "${barcode_list[@]}"; do 
    
    # Calculate the barcode directory path
    barcode_dir="barcode${barcode_number}"
    #Create log directory if it doesn't exist
    mkdir -p "$barcode_dir/logs"
    
    #Submit Jobs
    sub_job_id=$(sbatch \
                --mem "2G" \
                -c 1 \
                -o "$barcode_dir/logs/${barcode_number}_submission_log.txt" \
                --error "$barcode_dir/logs/${barcode_number}_submission.err" \
                --job-name="${barcode_number}_submission" \
                --wrap "/ei/projects/9/9742f7cc-c169-405d-bf27-cd520e26f0be/data/results/nanopore_PHIbase_analysis_scripts/single_barcode_process.sh $barcode_number $location $filter_length $reference_database" | awk '{print $4}')
    
 # Wait for the error file to be created and check if it's not empty
    if [ -s "$barcode_dir/logs/${barcode_number}_submission.err" ]; then
        echo "Error found in ${barcode_number}_submission.err. Cancelling job $sub_job_id"
        scancel "$sub_job_id"
        echo "Error message:" >> "$barcode_dir/logs/${barcode_number}_submission_log.txt"
        cat "$barcode_dir/logs/${barcode_number}_submission.err" >> "$barcode_dir/logs/${barcode_number}_submission_log.txt"
    fi
done

