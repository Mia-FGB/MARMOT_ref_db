import json
import pandas as pd

#A script to look at the species present in the download json output 

# Read the new JSON file
with open("test161224_download.json", "r") as file:
    json_data = json.load(file)

# Convert JSON data to a pandas DataFrame
new_json = pd.DataFrame(json_data)

#read older json file
with open("test131224_download.json", "r") as file:
    old_json_data = json.load(file)

# Convert JSON data to a pandas DataFrame
old_json = pd.DataFrame(old_json_data)

# Extract species from the JSON DataFrame (set means only unique)
taxid_new_json = set(new_json['organism_name'])
taxid_old_json = set(old_json['organism_name'])

#Script to compare original download.txt and download.json
# Read the download.txt file
# df_txt = pd.read_csv("download.txt", sep='\t', header=None, names=['species', 'taxid', 'assembly_accession', 'md5_url', 'fasta_url'])
# # Extract species from the download.txt DataFrame and trim
# df_txt['species'] = df_txt['species'].apply(lambda x: ' '.join(x.split()[:3]))
# taxid_txt = set(df_txt['species'])

# Find species that are in new JSON but not older
taxa_only_in_new_json = taxid_new_json - taxid_old_json
#and reverse
taxa_only_in_old_json = taxid_old_json - taxid_new_json

# Find species that are in download.txt but not in JSON
# taxa_only_in_txt = taxid_txt - taxid_json

# Print the results
print("Species only in new JSON:")
for species in taxa_only_in_new_json:
    print(species)

print(len(taxa_only_in_new_json))

print("\nSpecies only in old JSON:")
for species in taxa_only_in_old_json:
    print(species)
print(len(taxa_only_in_old_json))