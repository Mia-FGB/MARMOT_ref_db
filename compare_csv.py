import pandas as pd

def compare_csv(file1, file2):
    # Read the CSV files
    dfR = pd.read_csv(file1)
    dfPy = pd.read_csv(file2)

    # Rename the id column to taxid in file1
    dfR.rename(columns={'id': 'taxid'}, inplace=True)
    
    # Find common and unique taxids
    common_taxids = pd.merge(dfR,dfPy, on='taxid')
    unique_to_R = dfR[~dfR['taxid'].isin(dfPy['taxid'])]
    unique_to_Python = dfPy[~dfPy['taxid'].isin(dfR['taxid'])]
    
    # Print results
    print("Number of shared taxids:")
    print(len(common_taxids))
    print("\nUnique to R:")
    print(unique_to_R)
    print("\n Number of taxid unique to R:")
    print(len(unique_to_R))
    print("\nUnique to Python:")
    print(unique_to_Python)
    print("\n Number of taxid unique to Python:")
    print(len(unique_to_Python))

if __name__ == "__main__":
    file1 = 'unique_species_R.csv'
    file2 = 'unique_species_python.csv'
    compare_csv(file1, file2)