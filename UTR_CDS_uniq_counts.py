'''
This script processes multiple BAM files to count uniquely mapped reads for each transcript specified in a BED file.
The results are compiled into a count matrix and saved as a CSV file.

Steps:
1. Define a function to count uniquely mapped reads in a specified region of a BAM file.
2. Load the BED file containing transcript regions.
3. Initialize a dictionary to store read counts for each transcript across multiple BAM files.
4. Read the BED file and initialize the dictionary with transcript names.
5. Process each BAM file in the specified folder:
    a. Open the BAM file.
    b. Read the BED file and count reads for each transcript in the current BAM file.
    c. Close the BAM file.
6. Ensure all lists in the dictionary have the same length by appending zeros where necessary.
7. Create a DataFrame from the read counts dictionary.
8. Set the column names to be the BAM file names (without extensions).
9. Save the count matrix to a CSV file.

Usage:
- Ensure you have `pysam` and `pandas` installed: `pip install pysam pandas`
- Update the `bedfile` and `bam_folder` variables with the appropriate paths.
- Run the script to generate the count matrix CSV file.
'''


import os
import pysam
import pandas as pd

# Function to count uniquely mapped reads in a region
def count_unique_reads(bamfile, chrom, start, end, mapq_threshold=30):
    count = 0
    for read in bamfile.fetch(chrom, start, end):
        if not read.is_unmapped and read.mapping_quality >= mapq_threshold:
            count += 1
    return count

# Load BED file
#bedfile = "regions.bed"
bedfile = "longest.transcripts_with_5UTR_CDS.bed"

# Initialize a dictionary to store read counts for each transcript across multiple BAM files
read_counts = {}

# Read BED file and initialize the dictionary with transcript names
with open(bedfile, 'r') as bed:
    for line in bed:
        fields = line.strip().split()
        transcript_name = fields[0]
        read_counts[transcript_name] = []

# Folder containing BAM files
bam_folder = "./"
#bam_files = [f for f in os.listdir(bam_folder) if f.endswith(".No.rDNA.transcript.bam")]
bam_files = sorted([f for f in os.listdir(bam_folder) if f.endswith(".No.rDNA.transcript.bam")])


# Process each BAM file in the folder
for bam_filename in bam_files:
    bam_filepath = os.path.join(bam_folder, bam_filename)
    print(bam_filepath)
    bamfile = pysam.AlignmentFile(bam_filepath, "rb")

    # Read BED file and count reads for each transcript in the current BAM file
    with open(bedfile, 'r') as bed:
        for line in bed:
            fields = line.strip().split()
            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            transcript_name = fields[0]
            read_count = count_unique_reads(bamfile, chrom, start, end)
            read_counts[transcript_name].append(read_count)

    # Close the BAM file
    bamfile.close()

# Ensure all lists in the dictionary have the same length
for transcript in read_counts:
    while len(read_counts[transcript]) < len(bam_files):
        read_counts[transcript].append(0)

# Create a DataFrame from the read counts dictionary
#count_matrix = pd.DataFrame.from_dict(read_counts, orient='index', columns=[os.path.splitext(f)[0] for f in bam_files])

count_matrix = pd.DataFrame.from_dict(read_counts, orient='index', columns=[f.split('.')[0] for f in bam_files])
simplified_headers = [f.split('.')[0] for f in bam_files]
print(simplified_headers)

# Save the count matrix to a CSV file
count_matrix.to_csv("5UTR_CDS_count_matrix.csv")

print("Count matrix created and saved to 5UTR_CDS_count_matrix.csv")

