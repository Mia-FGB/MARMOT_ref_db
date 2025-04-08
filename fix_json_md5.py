import json

def fix_json_file(input_path, output_path):
    # Load the JSON data from the file
    with open(input_path, 'r') as f:
        json_data = json.load(f)
    
    # Loop through each item in the JSON list and fix the dlLinkMD5
    for item in json_data:
        if "dlLinkMD5" in item:
            item["dlLinkMD5"] = item["dlLinkMD5"].replace("https:https://", "https://")
    
    # Save the modified JSON back to a file
    with open(output_path, 'w') as f:
        json.dump(json_data, f, indent=4)

    print(f"Fixed JSON saved to {output_path}")

# Example usage: Provide your input and output file paths
input_file_path = "test161224_download.json"
output_file_path = "fix161224_download.json"

fix_json_file(input_file_path, output_file_path)


