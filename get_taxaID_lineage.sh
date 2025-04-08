#!/bin/bash
#A script to run on the all_taxaID_count.tsv file, but need to run off HPC to download taxonkit
#conda install -c bioconda taxonkit 
#https://bioinf.shenwei.me/taxonkit/usage/

#Name of sample, will change this with each analysis
sample=PHIbase_24hr

#Prep the taxaIDs_counts file
cut -f 1 -d $'\t' ${sample}_all_taxaID_count.tsv > ${sample}_taxaID.txt

#Get lineages
taxonkit lineage ${sample}_taxaID.txt > ${sample}_taxaID_lineage.txt

#Fill in blanks 
taxonkit reformat ${sample}_taxaID_lineage.txt -r Unassigned | cut -f 1,3 > ${sample}_taxaID_lineage_clean.txt

#Change delimiter
sed 's/;/\t/g' ${sample}_taxaID_lineage_clean.txt | awk -F'\t' 'BEGIN {OFS=","} { print $1, $2, $3, $4, $5, $6, $7, $8 }' > ${sample}_taxaID_lineage_sep.csv

#Add headers
echo "taxid,kindom,phylum,class,order,family,genus,species" > header.txt
cat header.txt ${sample}_taxaID_lineage_sep.csv > ${sample}_taxaID_lineage_sep_head.csv

rm ${sample}_taxaID.txt ${sample}_taxaID_lineage.txt ${sample}_taxaID_lineage_clean.txt  ${sample}_taxaID_lineage_sep.csv  header.txt
#Can now use this file in RStudio or pandas to plot
