#!/usr/bin/env python3

"""
Script: remove_random_genes.py
Description: Removes a specified number of random genes (and their transcripts/exons) from a GTF file.
Usage: python remove_random_genes.py -i genome.gtf -n 2 -o genome_filtered.gtf
"""

import argparse
import random
import sys

def parse_arguments():
    parser = argparse.ArgumentParser(description="Randomly remove a specified number of genes from a GTF file.")
    parser.add_argument('-i', '--input', required=True, help='Path to the input GTF file (e.g., genome.gtf)')
    parser.add_argument('-n', '--number', type=int, required=True, help='Number of genes to remove')
    parser.add_argument('-o', '--output', help='Path to the output filtered GTF file (default: input_filtered.gtf)')
    return parser.parse_args()

def extract_gene_ids(gtf_path):
    gene_ids = set()
    try:
        with open(gtf_path, 'r') as gtf:
            for line in gtf:
                if line.startswith('#'):
                    continue
                fields = line.strip().split('\t')
                if len(fields) < 9:
                    continue
                feature = fields[2]
                if feature == 'gene':
                    attributes = fields[8]
                    for attr in attributes.split(';'):
                        attr = attr.strip()
                        if attr.startswith('gene_id'):
                            gene_id = attr.split('"')[1]
                            gene_ids.add(gene_id)
                            break
    except FileNotFoundError:
        print(f"Error: File '{gtf_path}' not found.")
        sys.exit(1)
    except Exception as e:
        print(f"Error while reading GTF file: {e}")
        sys.exit(1)
    return list(gene_ids)

def main():
    args = parse_arguments()
    input_gtf = args.input
    num_remove = args.number
    output_gtf = args.output if args.output else f"{'.'.join(input_gtf.split('.')[:-1])}_filtered.gtf"

    # Extract gene IDs
    gene_ids = extract_gene_ids(input_gtf)
    total_genes = len(gene_ids)
    print(f"Total number of genes in the annotation: {total_genes}")

    # Validate number to remove
    if num_remove < 1 or num_remove > total_genes:
        print(f"Error: Number of genes to remove must be between 1 and {total_genes}.")
        sys.exit(1)

    # Select random genes to remove
    selected_genes = set(random.sample(gene_ids, num_remove))
    print(f"Selected genes to remove ({num_remove}):")
    for gene in selected_genes:
        print(gene)

    # Process GTF and write to output
    try:
        with open(input_gtf, 'r') as infile, open(output_gtf, 'w') as outfile:
            for line in infile:
                if line.startswith('#'):
                    outfile.write(line)
                    continue
                fields = line.strip().split('\t')
                if len(fields) < 9:
                    outfile.write(line)
                    continue
                attributes = fields[8]
                gene_id = None
                for attr in attributes.split(';'):
                    attr = attr.strip()
                    if attr.startswith('gene_id'):
                        gene_id = attr.split('"')[1]
                        break
                if gene_id and gene_id in selected_genes:
                    continue  # Skip lines with selected gene IDs
                outfile.write(line)
    except Exception as e:
        print(f"Error during GTF processing: {e}")
        sys.exit(1)

    print(f"Filtered GTF file created at '{output_gtf}'.")

if __name__ == "__main__":
    main()
