import os
import gzip
import shutil

# This script gzips uncompressed FASTA files in the specified directory, need them to be zipped for download.py to work.

# Set the directory where files are stored
download_dir = "download"

# File extensions to target
uncompressed_exts = [".fna", ".fa", ".fasta"]

# Loop through the download directory
for filename in os.listdir(download_dir):
    filepath = os.path.join(download_dir, filename)

    # Check if it's an uncompressed FASTA file
    if any(filename.endswith(ext) for ext in uncompressed_exts):
        gzipped_path = filepath + ".gz"

        # Skip if gzipped version already exists
        if os.path.exists(gzipped_path):
            print(f"Skipping {filename} (gzipped version already exists)")
            # delete the unzipped version 
            os.remove(filepath)
            print(f"Removed unzipped version: {filename}")
            continue

        # Gzip the file
        print(f"Gzipping {filename}...")
        with open(filepath, 'rb') as f_in, gzip.open(gzipped_path, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)

        # Optionally delete original
        os.remove(filepath)
        print(f"Compressed and removed: {filename}")
