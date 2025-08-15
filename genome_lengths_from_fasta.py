#this script doesn't work on the hpc without a singualrity container for BIO
import sys
from Bio import SeqIO  # Biopython library for parsing FASTA files

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python genome_lengths_from_fasta.py <input_fasta_file> <output_prefix>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_prefix = sys.argv[2]

def generate_genome_length_table(fasta_file, output_prefix):
    genome_lengths = {}
    output_file = f"{output_prefix}_genome_lengths.tsv"
    
    # Parse the FASTA file
    for record in SeqIO.parse(fasta_file, "fasta"):
        header_parts = record.description.split('|')
        if len(header_parts) >= 2 and header_parts[0] == "taxid":
            taxa_id = header_parts[1]
            seq_length = len(record.seq)
            # Aggregate length for this taxaID
            if taxa_id not in genome_lengths:
                genome_lengths[taxa_id] = 0
            genome_lengths[taxa_id] += seq_length
        else:
            print(f"SKIPPED: No taxaID/seq length identified for header: {record.description}")

    # Write the output to a file
    with open(output_file, "w") as out_f:
        out_f.write("taxaID\tgenome_length\n")
        for taxa_id, total_length in genome_lengths.items():
            out_f.write(f"{taxa_id}\t{total_length}\n")



generate_genome_length_table(input_file, output_prefix)
