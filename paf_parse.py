#!/usr/bin/python

#Now use lca parse instead of this script as it provides a better output
#A script to parse paf files from minimap output
#Filters for mapping quality 
#attempts to choose best mapping for reads which map to multiple taxa
#Should be run from directory which contains barcode directories

import sys, getopt, errno, os

#A class which contains all the info for each alignment 
class alignment_info:
    def __init__(self, q_name, taxaID, match_bases, map_length, MQ):
        self.q_name = q_name
        self.taxaID = taxaID
        self.match_bases = match_bases
        self.map_length = map_length
        self.MQ = MQ
    #Function within class that will print query name & MQ 
    def print_alignment(self):
        print(self.q_name, self.MQ)

#This block of code means it can take any paf file 
# When using the -b flag in the terminal will take barcodes
try:
    opts, args = getopt.getopt(sys.argv[1:],"b:")
except getopt.GetoptError:
    print("Option not recognised.")
    print("python my_script.py -b <barcode_number>")
    sys.exit(2)
for opt, arg in opts:
    if opt == "-b":
        input_arg = arg
    else:
        print("python my_script.py -i <barcode_number>")
        sys.exit(2)
barcode_number = input_arg

def main(barcode_number):
    pafFilename = None  # Define pafFilename with a default value
    barcode_dir = None 

    # Set the directory where output files should be saved
    barcode_dir = "./barcode{}".format(barcode_number)


    # Combine the barcode directory and the PAF file name to create the full path
    pafFilename = os.path.join(barcode_dir, "{}_mapped.paf".format(barcode_number))

    # Check if pafFilename is an empty string before attempting to open the file
    if not pafFilename:
        print("pafFilename is not valid. Exiting.")
        sys.exit(2)
    
    # Call parse_paf_file function to populate queries and taxa_count dictionaries
    queries, taxa_count = parse_paf_file(pafFilename)

    # Call the retain function with necessary arguments using queries & taxa_count from parse_paf_file function
    retain(pafFilename, barcode_number, barcode_dir, queries, taxa_count)

def parse_paf_file(pafFilename):

    queries = dict() #dictionary with queries as keys 
    taxa_count = dict() #dictionary which wil have taxaID for keys and counts as values 

    try:
        with open(pafFilename, 'r') as pafFile: 
            for line in pafFile:
                line = line.strip()
                fields = line.split('\t')
                q_name = fields[0].strip()
                #only want to get the numbers not the whole section
                taxaID =  fields[5].split("|")[1] #for PHIbase 
                #taxaID =  fields[5] #for other fasta references
                match_bases = int(fields[9].strip())
                map_length = int(fields[10].strip())
                MQ = int(fields[11].strip())
                #starting the dictionary to populate with taxaIDs
                taxa_count[taxaID] = 0
                #only want to include MQ which are 0 or >=5 so filter on this line before creating the dictionary 
                if MQ >= 5 or MQ == 0:
                    #Now calling the class alignment_info and inputting this info
                    ai = alignment_info(q_name, taxaID, match_bases, map_length, MQ)
                    #whilst looping through each line in the paf file
                    if q_name not in queries.keys():
                        #adding a new key for new queries
                        queries[q_name] = list()
                    #then appending the class to the dictionary based on queryName
                    queries[q_name].append(ai)
    
    except (OSError, IOError) as e: 
        if getattr(e, 'errno', 0) == errno.ENOENT:
            print ("Could not find file " + pafFilename)
            sys.exit(2)
        else:
            print("An error occurred while parsing the PAF file.")
            sys.exit(2)
    
    # Return the updated queries and taxa_count dictionaries
    return queries, taxa_count
    


#Function to select which alignments are retained
def retain(pafFilename, barcode_number, barcode_dir, queries, taxa_count):
    single_query = []
    multiple_query = []
    #counting how many reads are ignored
    ignored_reads = 0
    #looping through the queries dictionary 
    for q_name in queries.keys():
        #creating a list that contains the info from the dictionary 
        ai_list = queries[q_name]
        #If there is only one input for a query name then the read is unique 
        if len(ai_list) == 1:
            #so it can be added to the unique list 
            single_query.append(q_name)
            #this ai list only contains one thing so we can get it with [0]
            ai = ai_list[0]
            #getting the taxaID out of the class we set up using the ai_list query name
            taxa_id = ai.taxaID
            #updating the count for the specific taxaID in the dictionary
            #taxa_ID is the key
            taxa_count[taxa_id] += 1
        else:
            #If there is more than one entry then there must be multiple reads
            multiple_query.append(q_name)
            #making a set as this can only contain unique values 
            mapping_qualities = set()
            multiple_taxa_ID = set()
            #looping through list of dict info
            for ai in ai_list:
                #adding the mapping quality to the set
                mapping_qualities.add(ai.MQ)
                multiple_taxa_ID.add(ai.taxaID)

            #1x taxa, add as 1 alignment
            if len(multiple_taxa_ID) == 1:
                #extracting from the set to use it in dict
                taxa_id = list(multiple_taxa_ID)[0]
                taxa_count[taxa_id] +=1
            
            #1xMQ multi x taxa
            #currenlty going to ignore these reads but want to print them to see how big of a problem this is 
            #In the future could look at LCA 
            if len(mapping_qualities) == 1 and len(multiple_taxa_ID) > 1:
                # Create the path to the ignored_reads.txt file within the barcode_dir
                ignored_file_path = os.path.join(barcode_dir, "{}_ignored_reads_query_ID.tsv".format(barcode_number))
                
                # Open the file in 'a' (append) mode
                with open(ignored_file_path, 'a') as ignored_file:
                    ignored_file.write("This query: {} has the same MQ {} but maps to different taxaIDs: {} so has been ignored.\n".format(q_name, mapping_qualities, multiple_taxa_ID))
            ignored_reads += 1

            #multi x MQ and multi x taxaID, just want to take taxaID from highest MQ
            if len(mapping_qualities) > 1 and len(multiple_taxa_ID) > 1:
                #sort and choose the highest MQ to add to taxa_count dictionary
                maxMQ = 0
                taxaIDForMaxMQ = ""
                for ai in ai_list:
                    #looping though the ai list which contains all info
                    #Still within the q_name loop
                    if ai.MQ > maxMQ:
                        maxMQ = ai.MQ
                        taxaIDForMaxMQ = ai.taxaID
                        taxa_count[taxaIDForMaxMQ] +=1
   
    # Print the dictionary as a list of taxaIDs & counts to a separate file
    # Create the path to the ignored_reads.txt file within the barcode_dir
    taxa_file_path = os.path.join(barcode_dir, "{}_taxaID_counts.tsv".format(barcode_number))
    # Open the file in 'a' (append) mode
    with open(taxa_file_path, 'a') as taxa_file:
        for taxa_id, count in taxa_count.items():
            taxa_file.write('{}\t{}\t{}\n'.format(taxa_id, count, barcode_number))

    # Print the number of ignored reads to a file within the barcode_dir
    ignored_summary_file_path = os.path.join(barcode_dir, "{}_number_ignored_reads.tsv".format(barcode_number))
    with open(ignored_summary_file_path, 'w') as ignored_summary_file:
        ignored_summary_file.write('{}\t{}'.format(ignored_reads, barcode_number))


if __name__ == "__main__":
    main(barcode_number)
