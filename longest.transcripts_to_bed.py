# Read input data from a file and convert it to BED format with an additional 0 at the start
input_file = "longest.transcripts.info.txt"
output_file = "longest.transcripts.bed"

with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
    # Skip the header line
    next(infile)
    for line in infile:
        fields = line.strip().split()
        transcript = fields[1]
        start = 0
        end = fields[9]
        name = fields[4]
        gene_id = fields[3]
        strand = fields[2]
        outfile.write(f"{transcript}\t{start}\t{end}\t{name}\t{gene_id}\t{strand}\t\n")

print(f"BED file created: {output_file}")


