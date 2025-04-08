# MARMOT_ref_db

A repository for scripts to generate the reference database and associated files, to be run on local not HPC
Linked to /Users/berelsom/Library/CloudStorage/OneDrive-NorwichBioScienceInstitutes/Pathogen_Database/scripts

## Step 1: Create New Directory

**Location:**  
```bash
cd /Users/berelsom/Library/CloudStorage/OneDrive-NorwichBioScienceInstitutes/Pathogen_Database
mkdir Pathogen_Database_DDMMYY
```

> Replace `DDMMYY` with the current date (e.g., `080524`)

---

## Step 2: Download Latest CSV Files

### 🧬 PHI-base:
- Go to: https://github.com/PHI-base/data/tree/master/releases
- Download `phibase_current.csv`
- **Rename** to: `phibase_vers.csv` (e.g., `phibase_4.17.csv`)
- Use the release **month & date** in your filenames

### 🌿 DEFRA Risk Register:
- Source: **UK Plant Health Risk Register**
- Save as: `Risk_Register_MMYY.csv`

---

## Step 3: Run the Scripts

**Navigate to working directory:**
```bash
cd /Users/berelsom/Library/CloudStorage/OneDrive-NorwichBioScienceInstitutes/Pathogen_Database
```

**Activate conda environment:**
```bash
conda activate pathogen-database
```

### Make the JSON:
```bash
python scripts/Make_Pathogen_Database.py \
  --phibase Pathogen_Database_Test/phibase_test.csv \
  --risk_register Pathogen_Database_Test/Risk_Register_Test.csv \
  --output test_download
```

### Download and Build Database:
```bash
python scripts/download.py --input test_download --date MMYYYY
```

---

## HPC Upload Instructions

Once complete, upload results to HPC at:
```text
/ei/projects/9/9742f7cc-c169-405d-bf27-cd520e26f0be/data/results/phibase
```

### Naming:
- Database file: `pathogen_database_MMYYYY.fa`
- Genome lengths: `MMYYYY_genome_lengths.tsv`

Move **older versions** into the `old_pathogen_database/` folder.

