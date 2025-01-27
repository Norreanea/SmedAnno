# Define input and output file paths
input_file = "/mnt/e/Stringtie_anno/SM_anno/ncRNA/ncRNA_seq_dedup_underscore_short_header_linear.fa"
output_file = "/mnt/d/Elac2/final_results/tRF_target/RNAcentral_ncRNA_lengths.txt"

# Initialize variables
ncRNA_lengths = []

# Process the FASTA file
with open(input_file, "r") as infile:
    current_name = None
    current_seq = ""
    for line in infile:
        if line.startswith(">"):
            # Process the previous sequence
            if current_name:
                ncRNA_lengths.append((current_name, len(current_seq)))
            # Start a new sequence
            current_name = line[1:].strip()  # Remove ">" and whitespace
            current_seq = ""
        else:
            # Accumulate sequence lines
            current_seq += line.strip()
    # Process the last sequence
    if current_name:
        ncRNA_lengths.append((current_name, len(current_seq)))

# Write results to the output file
with open(output_file, "w") as outfile:
    for name, length in ncRNA_lengths:
        outfile.write(f"{name}\t{length}\n")

print(f"Sequence lengths have been written to {output_file}")
