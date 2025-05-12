#!/bin/bash
#SBATCH -J emerge_pipeline                    # Job name
#SBATCH -o logs/emerge_pipeline_sub_%a.out     # Output for each job
#SBATCH -e logs/emerge_pipeline_sub_%a.err     # Error for each job
#SBATCH -p ei-medium                            # may need to make this ei-long with full pipeline, wait and see
#SBATCH -c 1                                    # Number of CPU cores
#SBATCH --mem=2G                                # Memory size - optimise later
#SBATCH --mail-type=END,FAIL                    # Notifications for job done & fail 
#SBATCH --mail-user=mia.berelson@earlham.ac.uk  # Email address 
#SBATCH --array=01-88%20                        # Array range and concurrency limit (specify this based on the number of barcodes)

#Set up variables
location="/ei/projects/9/9742f7cc-c169-405d-bf27-cd520e26f0be/data/raw/CF_2023_combined_sup_runs/" #raw read location - level above fastqpass
filter_length="300" #length to filter reads to be greater than
reference_database="/ei/projects/9/9742f7cc-c169-405d-bf27-cd520e26f0be/data/results/phibase/pathogen_database_080524.fa"
scratch_dir="/ei/projects/9/9742f7cc-c169-405d-bf27-cd520e26f0be/scratch/church_farm_2023" #where the concatenated reads will be output

# Create txt files if they don't exist
touch ./percent_reads_retained_length_filter.txt
touch ./no_fail_reads.txt
touch ./no_reads_ignored_parse_filter.txt
touch ./all_taxaID_count.tsv
touch ./lcaparse_summary.txt
touch ./lcaparse_perread.txt

# Add a header line to lcaparse_summary.txt if it doesn't already exist
grep -q "Barcode\tRead_Count\tPercentage_of_Reads\tTaxon_ID\tTaxon_Path\tTaxon_Rank" ./lcaparse_summary.txt || printf "Barcode\tRead_Count\tPercentage_of_Reads\tTaxon_ID\tTaxon_Path\tTaxon_Rank\n" > ./lcaparse_summary.txt
# Add a header line to lcaparse_perread.txt if it doesn't already exist
grep -q "Barcode\tRead_ID\tTaxon_ID\tTaxon_Name\tTaxon_Rank\tMean_Identity\tMaxMeanIdentity" ./lcaparse_perread.txt || printf "Barcode\tRead_ID\tTaxon_ID\tTaxon_Name\tTaxon_Rank\tMean_Identity\tMaxMeanIdentity\n" > ./lcaparse_perread.txt

# Define the list of barcode numbers
barcode_list=($(seq -w 01 88))

# Get the barcode corresponding to this task ID
barcode_number="${barcode_list[$((SLURM_ARRAY_TASK_ID - 1))]}"

# Now barcode_number holds the specific barcode to use for this task
echo "Processing barcode: $barcode_number"

# Set up directory for the current barcode
barcode_dir="barcode${barcode_number}"

# Execute the main processing script - commented out for now while I work on the lookup table
 /ei/projects/9/9742f7cc-c169-405d-bf27-cd520e26f0be/data/results/nanopore_PHIbase_analysis_scripts/single_barcode_process.sh \
     $barcode_number $location $filter_length $reference_database $scratch_dir

# Check if the job script encountered an error and handle cancellation
if [ $? -ne 0 ]; then
    echo "Error detected for barcode ${barcode_number}. Cancelling job ${SLURM_JOB_ID}."
    scancel "$SLURM_JOB_ID"
fi
