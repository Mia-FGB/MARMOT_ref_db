#!/bin/bash

# Script that goes through the steps to build a reference database for pathogen detection.
# Usage example:
# ./scripts/build_reference_database.sh \
#   -p Pathogen_Database_Test/phibase_test.csv \
#   -r Pathogen_Database_Test/Risk_Register_Test.csv \
#   -d 042024 \
#   -o Pathogen_Database_Test

# Parse Arguments
while getopts ":p:r:d:o:" opt; do
  case ${opt} in
    p )
      PHIBASE_CSV="$OPTARG"
      ;;
    r )
      RISK_REGISTER_CSV="$OPTARG"
      ;;
    d )
      DATE_TAG="$OPTARG"
      ;;
    o )
      OUTDIR="$OPTARG"
      ;;
    \? )
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
    : )
      echo "Option -$OPTARG requires an argument." >&2
      exit 1
      ;;
  esac
done

# Check arguments provided
if [[ -z "$PHIBASE_CSV" || -z "$RISK_REGISTER_CSV" || -z "$DATE_TAG" ]]; then
  echo "Usage: $0 -p path/to/phibase.csv -r path/to/risk_register.csv -d MMYYYY -o path/to/output_directory"
  exit 1
fi

# Display arguments
echo "PHI-base CSV:       $PHIBASE_CSV"
echo "Risk Register CSV:  $RISK_REGISTER_CSV"
echo "Date Tag:           $DATE_TAG"


# Check if the input files exist
if [[ ! -f "$PHIBASE_CSV" ]]; then
  echo "Error: PHI-base CSV file not found: $PHIBASE_CSV"
  exit 1
fi
if [[ ! -f "$RISK_REGISTER_CSV" ]]; then
  echo "Error: Risk register CSV file not found: $RISK_REGISTER_CSV"
  exit 1
fi

# Create output directory if it doesn't exist
echo "Creating output directory: $OUTDIR"
mkdir -p "$OUTDIR"
mkdir -p "$OUTDIR/logs"

# Set up logging
LOGFILE="$OUTDIR/logs/build_reference_database.log"
exec > >(tee -a "$LOGFILE") 2>&1 # Redirect stdout and stderr to log file
echo "Logging to $LOGFILE"

# Activate conda environment
echo "Activating conda environment..."
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate pathogen_database
if [[ $? -ne 0 ]]; then
  echo "Error: Failed to activate conda environment."
  exit 1
fi
echo "Conda environment activated."


# Generate risk table 
python scripts/generate_risk_table.py -i $RISK_REGISTER_CSV -o $OUTDIR/risk_table.csv
echo "Risk table generated: risk_table.csv"

# Run Make_Pathogen_Database.py - Takes a long time 
echo "Generating pathogen json list ..."

python scripts/Make_Pathogen_Database.py \
  --phibase "$PHIBASE_CSV" \
  --risk_register "$RISK_REGISTER_CSV" \
  --output "$OUTDIR/download_input_${DATE_TAG}"

echo "Pathogen database generated: download_input_${DATE_TAG}"

# Download genomes & Build database
echo "Downloading genomes and building database pathogen_database_${DATE_TAG}.fa..."
# Overwrite the existing output database if there is one
> "$OUTDIR/pathogen_database_${DATE_TAG}.fa"

python scripts/download.py \
  -i "$OUTDIR/download_input_${DATE_TAG}" \
  -d "$DATE_TAG" \
  -o "$OUTDIR"
echo "Genomes downloaded and database built, output: pathogen_database_${DATE_TAG}.fa"

# Generate lengths table 
echo "Generating genome lengths table..."
python scripts/genome_lengths_from_fasta.py $OUTDIR/pathogen_database_${DATE_TAG}.fa $OUTDIR/${DATE_TAG}

echo "Generated genome lengths table: ${DATE_TAG}_genome_lengths.tsv"