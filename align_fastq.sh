#!/bin/bash

# Define the array of URLs
URL_LIST=(
    "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR106/031/SRR10691631/SRR10691631_1.fastq.gz"
    "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR106/031/SRR10691631/SRR10691631_2.fastq.gz"
)

# Define the fastq files, output and STAR genome directories
FASTQ_DIR=/fastq_folder
OUTPUT_DIR=/output_folder
STAR_GENOME=/yourstargenomefolder

# Create the directory if it doesn't exist
mkdir -p "$FASTQ_DIR"

for URL in "${URL_LIST[@]}"; do
    # Extract filename from URL
    FILE_NAME=$(basename "$URL")
    
    # Download the fastq file
    curl -L "$URL" -o "$FASTQ_DIR/$FILE_NAME"
    
    # Check if the download was successful

	if [ $? -eq 0 ]; then
    echo "Download of $FILE_NAME completed successfully."
else
    echo "Error in downloading $FILE_NAME."
fi
done

for R1 in ${FASTQ_DIR}/*_1.fastq.gz
 do
 	# Derive R2 file by replacing "_1.fastq.gz" with "_2.fastq.gz"
 	R2=${R1/_1.fastq.gz/_2.fastq.gz}
	# Extract the sample name (e.g., SRR10691631 from SRR10691631_1.fastq.gz)
 	sample_name=$(basename ${R1} _1.fastq.gz)
	# Define the output prefix
     output_prefix="${OUTPUT_DIR}/${sample_name}."
	# Run STAR alignment
 	STAR --runThreadN $ncpus \
      	--genomeDir ${STAR_GENOME} \
      	--outFileNamePrefix ${output_prefix} \
      	--readFilesIn ${R1} ${R2} \
      	--readFilesCommand zcat \
      	--quantMode GeneCounts \
      	--twopassMode Basic
 done