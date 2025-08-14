#!/bin/bash
#script to run lcaparse on the minimap output
# Check for the correct number of arguments
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <barcode>" 
    exit 1
fi

# Extract the arguments
barcode_number="$1"
barcode_dir="barcode${barcode_number}"

# Load marti - will need to update with new marti releases in the future
source package /tgac/software/testing/bin/marti-0.9.19

# Run the lcaparse command
lcaparse -input $barcode_dir/${barcode_number}_mapped.paf \
-output $barcode_dir/${barcode_number}_lcaparse \
-taxonomy /ei/projects/9/9742f7cc-c169-405d-bf27-cd520e26f0be/data/results/taxonomy/taxdmp_2024-04-01 \
-mapfile /ei/projects/9/9742f7cc-c169-405d-bf27-cd520e26f0be/data/results/lcaparse/mapfiles/accession_map_2024-04_eukaryota.txt \
-format PAF \
-minidentity 85 \
-minlength 150 
