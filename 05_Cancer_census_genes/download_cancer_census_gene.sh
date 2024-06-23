#! /bin/bash

BASE_PATH="/mnt/c/Users/CEEL-PC-005/Desktop/Joon/Final_scripts/Cancer_genes_analysis_of_APOBEC_motifs_test/" # CHANGE HERE

# Download cosmic caner gene census 
AUTHSTRING=$(echo "joon-hyun.song@stonybrook.edu:Syl97er%\$812" | base64) # CHANGE HERE
downloadlink=$( curl -H "Authorization: Basic ${AUTHSTRING}" "https://cancer.sanger.ac.uk/api/mono/products/v1/downloads/scripted?path=grch38/cosmic/v99/Cosmic_CancerGeneCensus_Tsv_v99_GRCh38.tar&bucket=downloads" )
echo "Downloading cancer_gene_census.csv"
curl -o ${BASE_PATH}cancer_gene_census.csv ${downloadlink}