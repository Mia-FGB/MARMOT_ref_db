import pandas as pd
import argparse

#Script to generate a simplified table from the DEFRA risk Register
# python generate_risk_table.py -i /path/to/input.csv -o /path/to/output.csv

# Set up argument parsing
parser = argparse.ArgumentParser(description="Generate a risk table from a CSV file.")
parser.add_argument("-i", "--input", required=True, help="Path to the input DEFRA CSV file")
parser.add_argument("-o", "--output", default="risk_table.csv", help="Path to the output CSV file (default: risk_table.csv)")
args = parser.parse_args()

# Read the input CSV file
risk = pd.read_csv(args.input)

# Fill NaN values with empty strings
risk['Type of pest'] = risk['Type of pest'].fillna("")  
# Define the list of items to remove
remove = ["Insect", "Mite", "Nematode", "Plant"] 
# Filter rows where 'Type of pest' is not in the remove list
risk = risk[~risk['Type of pest'].isin(remove)] 
# Remove single quotes from 'Pest Name' column if present
risk['Pest Name'] = risk['Pest Name'].str.replace("'", "")

# Selecting columns
columns_keep = [
    'Type of pest', 'Pest Name', 'EU and EPPO listing', 'UK', 'Pathways', 'Likelihood', 'Impact ',
    'UK Relative Risk Rating (unmitigated)', 'Regulation', 'Likelihood.1', 'Impact .1',
     'UK Relative Risk Rating (mitigated)', 'Scenario for Risk Register'
    ]
risk = risk[columns_keep]

# Create a new column 'Species' with only the Genus and Species
risk['Species'] = risk['Pest Name'].apply(lambda x: ' '.join(x.split()[:2]))

# Create a new column 'Regulated' based on the condition
risk['Regulated'] = risk['EU and EPPO listing'].apply(
    lambda x: 'Yes' if isinstance(x, str) and 'regulated quarantine pest' in x.lower() else 'No'
)

# Create a new column 'Natural Spread' based on the condition
risk['Natural_Spread'] = risk['Pathways'].apply(
    lambda x: 'Yes' if isinstance(x, str) and 'natural spread' in x.lower() else 'No'
)

# Drop the old columns
risk = risk.drop(columns=['EU and EPPO listing', 'Pathways'])

 # Rename some columns 
risk = risk.rename(columns={
    'Type of pest': 'Type_of_pest',
    'Pest Name': 'Pest_Name',
    'Likelihood': 'Likelihood_unmitigated',
    'Likelihood.1': 'Likelihood_mitigated',
    'Impact ': 'Impact_unmitigated',
    'Impact .1': 'Impact_mitigated',
    'UK Relative Risk Rating (unmitigated)': 'Risk_Rating_unmitigated',
    'UK Relative Risk Rating (mitigated)': 'Risk_Rating_mitigated',
    'Scenario for Risk Register': 'Scenario_for_Risk_Register'
})

# Output the DataFrame to a CSV file
risk.to_csv(args.output, index=False)
print(f"Risk table generated and saved to {args.output}")