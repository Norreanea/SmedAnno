import pandas as pd
import re
import os
import sys

def classify_ncRNA(query_id):
    """Classify ncRNA type based on query ID patterns."""
    patterns = [
        (r'tRNA|transfer|UNDET_|leader', 'tRNA'),
        (r'rRNA|ribosomal|LS|SS|ITS|28|Spacer', 'rRNA'),
        (r'miRNA|Sme-|micro|miR|sme-lin|sme-let|bantam', 'miRNA'),
        (r'long|lncRNA|LNC|lincRNA', 'lncRNA'),
        (r'piRNA|piwi', 'piRNA'),
        (r'snoRNA|nucleolar', 'snoRNA'),
        (r'snRNA|nuclear|spliceosomal', 'snRNA'),
        (r'recognition', 'SRP_RNA')
    ]
    for pattern, rna_type in patterns:
        if re.search(pattern, query_id, re.IGNORECASE):
            return rna_type
    return 'other'

# Coverage thresholds (adjust as needed)
COVERAGE_THRESHOLDS = {
    'tRNA': 0.95,
    'rRNA': 0.90,
    'miRNA': 0.90,
    'lncRNA': 0.70,
    'piRNA': 0.80,
    'snoRNA': 0.80,
    'snRNA': 0.80,
    'SRP_RNA': 0.80,
    'other': 0.80
}

# Filtering logic with debugging
def filter_blast(input_path, output_path):
    # Debug: Print current working directory
    print(f"Current working directory: {os.getcwd()}")
    
    # Debug: Check if input file exists
    if not os.path.isfile(input_path):
        print(f"Error: Input file '{input_path}' does not exist.")
        sys.exit(1)
    else:
        print(f"Input file '{input_path}' found. Proceeding with filtering.")
    
    # Attempt to open the input and output files
    try:
        with open(input_path, 'r') as fin, open(output_path, 'w') as fout:
            total_lines = 0
            written_lines = 0
            for line_number, line in enumerate(fin, start=1):
                total_lines += 1
                cols = line.strip().split('\t')
                
                # Debug: Check for correct number of columns
                if len(cols) < 13:
                    print(f"Warning: Line {line_number} skipped due to insufficient columns ({len(cols)} columns found).")
                    continue  # Skip malformed lines
                
                # Parse key columns with error handling
                try:
                    query_id = cols[0]
                    pident = float(cols[2])
                    alen = int(cols[3])
                    evalue = float(cols[10])
                    qlen = int(cols[12])
                except ValueError as ve:
                    print(f"Warning: Line {line_number} skipped due to value error: {ve}")
                    continue  # Skip lines with invalid numerical values
                
                # Basic filters: E-value ≤1e-20, %ID ≥95
                if evalue > 1e-20 or pident < 95:
                    # Debug: Line filtered out by basic criteria
                    # Uncomment the next line to enable verbose filtering messages
                    # print(f"Info: Line {line_number} filtered out (E-value: {evalue}, %ID: {pident}).")
                    continue
                
                # Prevent division by zero
                if qlen == 0:
                    print(f"Warning: Line {line_number} skipped due to qlen=0.")
                    continue
                
                # Calculate coverage of the query
                coverage = alen / qlen
                
                # Classify ncRNA and get threshold
                rna_type = classify_ncRNA(query_id)
                threshold = COVERAGE_THRESHOLDS.get(rna_type, 0.80)
                
                # Apply type-specific coverage filter
                if coverage >= threshold:
                    fout.write(line)
                    written_lines += 1
                else:
                    # Debug: Line filtered out by coverage criteria
                    # Uncomment the next line to enable verbose filtering messages
                    # print(f"Info: Line {line_number} filtered out (Coverage: {coverage:.2f} < Threshold: {threshold}).")
                    pass
                
                # Periodically print progress
                if line_number % 100000 == 0:
                    print(f"Processed {line_number} lines. Written {written_lines} lines so far.")
                    
            # Final debug information
            print(f"Filtering completed. Total lines processed: {total_lines}. Total lines written: {written_lines}.")
    
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        sys.exit(1)

# Run the filter with specified paths
if __name__ == "__main__":
    input_file = '/mnt/d/Elac2/final_results/tRF_target/blast_results/prefiltered_blast_with_lengths.txt'
    output_file = '/mnt/d/Elac2/final_results/tRF_target/blast_results/filtered_blast_type_specific.txt'
    
    filter_blast(input_file, output_file)
