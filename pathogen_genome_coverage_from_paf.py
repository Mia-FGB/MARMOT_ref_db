#!/usr/bin/python

#A script to calculate genome coverage for specific genera from minimap output
#Filters on identity and coverage
#Should be run from directory which contains barcode directories
#Uses a genome_lengths_file generated with genome_lengths_from_fasta.py script using the reference database used for mapping

import sys, os, csv

#Opening a tab delimited file with taxaID and reference genome length
def load_genome_lengths(genome_lengths_file):
    genome_lengths = {} #Dictionary to store genome lengths
    with open(genome_lengths_file, 'r') as genome_lengths_f:
        reader = csv.reader(genome_lengths_f, delimiter='\t')
        next(reader)  # Skip header row
        for row in reader:
            taxaID = row[0]
            length = int(row[1])  # Genome length in base pairs
            genome_lengths[taxaID] = length #populate dictionary 
    return genome_lengths

def process_paf_file(pafFilename, genome_lengths): #Function to process paf file one row at a time
    taxa_mapped_bases = {}
    filtered_reads = []
    multi_taxa_reads = {}
    processed_reads = set()
    taxa_read_ids = {}

    with open(pafFilename, 'r') as paf_file:
        reader = csv.reader(paf_file, delimiter='\t')
        for row in reader:
            read_id = row[0]
            q_length = int(row[1])
            a_start = int(row[2])
            a_end = int(row[3])
            taxaID = row[5].split("|")[1]
            matching_bases = int(row[9])
            a_length = int(row[10])
            identity = (matching_bases / a_length) * 100
            coverage = ((a_length) / q_length) * 100

            # See how many alignments have coverage > 110
            # if coverage > 110:
            #     print (f"Warning: coverage for read {read_id} is greater than 110% ({coverage:.2f}%)")

            if coverage >= 80 and identity >= 80: #change as needed
                if (read_id, taxaID) not in processed_reads: #only count each read aligning to same taxa once
                    if taxaID not in taxa_mapped_bases:
                        taxa_mapped_bases[taxaID] = 0 #initialize taxaID in dictionary
                        taxa_read_ids[taxaID] = set() 
                    taxa_mapped_bases[taxaID] += q_length
                    taxa_read_ids[taxaID].add(read_id)
                    processed_reads.add((read_id, taxaID))

                    if read_id not in multi_taxa_reads:
                        multi_taxa_reads[read_id] = []
                    multi_taxa_reads[read_id].append((taxaID, identity, coverage))
            else:
                filtered_reads.append((read_id, taxaID, q_length, a_length, matching_bases, identity, coverage))

    return taxa_mapped_bases, filtered_reads, multi_taxa_reads, taxa_read_ids

def write_genome_coverage(genome_coverage_file, taxa_mapped_bases, genome_lengths, taxa_read_ids):
    with open(genome_coverage_file, "w") as mapped_file:
        mapped_file.write("taxaID\tmapped_bases\tgenome_length\tcoverage_percentage\tnum_reads\tread_ids\n")
        for taxaID, bases in taxa_mapped_bases.items():
            genome_length = genome_lengths.get(taxaID, 0)
            if genome_length > 0:
                coverage_percentage = (bases / genome_length) * 100
            else:
                print(f"Warning: taxaID {taxaID} not found in genome lengths table.")
                coverage_percentage = 'N/A'
            num_reads = len(taxa_read_ids[taxaID])
            read_ids = ",".join(taxa_read_ids[taxaID])
            mapped_file.write(f"{taxaID}\t{bases}\t{genome_length}\t{coverage_percentage:.4f}\t{num_reads}\t{read_ids}\n")

def write_filtered_reads(filtered_reads_file, filtered_reads):
    with open(filtered_reads_file, "w") as filtered_file:
        filtered_file.write("read_id\ttaxaID\tread_length\talignment_length\tmatching_bases\tidentity\tcoverage\n")
        for read in filtered_reads:
            filtered_file.write("\t".join(map(str, read)) + "\n")

def write_multi_taxa_reads(multi_taxa_reads_file, multi_taxa_reads):
    with open(multi_taxa_reads_file, "w") as multi_taxa_file:
        multi_taxa_file.write("read_id\ttaxaID\tidentity\tcoverage\n")
        for read_id, taxa_info in multi_taxa_reads.items():
            if len(taxa_info) > 1:
                for taxaID, identity, coverage in taxa_info:
                    multi_taxa_file.write(f"{read_id}\t{taxaID}\t{identity:.2f}\t{coverage:.2f}\n")

def main():
    if len(sys.argv) != 3:
        print("Usage: python pathogen_genome_coverage_from_paf.py <barcode_number> <genome_lengths_file>")
        sys.exit(1)

    barcode_number = sys.argv[1]
    genome_lengths_file = sys.argv[2]

    barcode_dir = f"./barcode{barcode_number}"
    pafFilename = os.path.join(barcode_dir, f"{barcode_number}_mapped.paf")

    if not pafFilename:
        print("pafFilename is not valid. Exiting.")
        sys.exit(2)

    genome_coverage_file = os.path.join(barcode_dir, f"{barcode_number}_genome_coverage.txt")
    filtered_reads_file = os.path.join(barcode_dir, f"{barcode_number}_coverage_excluded_reads.txt")
    multi_taxa_reads_file = os.path.join(barcode_dir, f"{barcode_number}_coverage_multi_taxa_reads.txt")

    genome_lengths = load_genome_lengths(genome_lengths_file)
    taxa_mapped_bases, filtered_reads, multi_taxa_reads, taxa_read_ids = process_paf_file(pafFilename, genome_lengths)

    write_genome_coverage(genome_coverage_file, taxa_mapped_bases, genome_lengths, taxa_read_ids)
    write_filtered_reads(filtered_reads_file, filtered_reads)
    write_multi_taxa_reads(multi_taxa_reads_file, multi_taxa_reads)

    print(f"Processed genome coverage for barcode{barcode_number}")

if __name__ == "__main__":
    main()