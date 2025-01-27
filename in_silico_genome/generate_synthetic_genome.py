#!/usr/bin/env python3

import random
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

# Configuration Parameters
CHROMOSOME_NAME = "chr1"
CHROMOSOME_LENGTH = 1000000  # 1,000,000 bp
NUM_GENES = 100
MIN_GENE_LENGTH = 1000  # Minimum length of a gene
MAX_GENE_LENGTH = 5000  # Maximum length of a gene
MIN_EXONS = 2
MAX_EXONS = 10
MIN_EXON_LENGTH = 100
MAX_EXON_LENGTH = 1000
MIN_INTRON_LENGTH = 100  # Minimum intron length
MAX_INTRON_LENGTH = 1000  # Maximum intron length
INTERGENIC_MIN = 500  # Minimum distance between genes

# Nucleotide bases
NUCLEOTIDES = ['A', 'T', 'C', 'G']

def generate_random_sequence(length):
    return ''.join(random.choices(NUCLEOTIDES, k=length))

def generate_genome():
    # Initialize the genome with random nucleotides
    genome_seq = generate_random_sequence(CHROMOSOME_LENGTH)
    genome_record = SeqRecord(Seq(genome_seq), id=CHROMOSOME_NAME, description="")
    return genome_record

def assign_genes():
    genes = []
    current_pos = 1  # 1-based coordinates

    for gene_id in range(1, NUM_GENES + 1):
        # Assign number of exons first based on gene length constraints
        # To ensure that gene_length >= (MIN_EXON_LENGTH * num_exons) + (MIN_INTRON_LENGTH * (num_exons -1))
        # We need to select num_exons such that gene_length >= 200*num_exons -100

        # Initially assign a random gene length
        gene_length = random.randint(MIN_GENE_LENGTH, MAX_GENE_LENGTH)
        
        # Determine the maximum possible number of exons for the current gene_length
        max_possible_exons = min(
            MAX_EXONS,
            (gene_length + MIN_INTRON_LENGTH) // (MIN_EXON_LENGTH + MIN_INTRON_LENGTH)
        )
        # Ensure that max_possible_exons is at least MIN_EXONS
        if max_possible_exons < MIN_EXONS:
            max_possible_exons = MIN_EXONS
            # Recalculate gene_length to satisfy the minimum requirement
            min_required_length = (MIN_EXON_LENGTH * MIN_EXONS) + (MIN_INTRON_LENGTH * (MIN_EXONS - 1))
            if gene_length < min_required_length:
                gene_length = min_required_length
        
        # Now, assign num_exons within the feasible range
        num_exons = random.randint(MIN_EXONS, max_possible_exons)
        
        # Recalculate max_possible_exons based on the potentially adjusted gene_length
        max_possible_exons = min(
            MAX_EXONS,
            (gene_length + MIN_INTRON_LENGTH) // (MIN_EXON_LENGTH + MIN_INTRON_LENGTH)
        )
        num_exons = min(num_exons, max_possible_exons)  # Ensure num_exons does not exceed feasible number
        
        # Re-verify gene_length is sufficient
        min_required_length = (MIN_EXON_LENGTH * num_exons) + (MIN_INTRON_LENGTH * (num_exons - 1))
        if gene_length < min_required_length:
            gene_length = min_required_length

        # Now, proceed to assign exon lengths
        # Generate initial exon lengths
        exon_lengths = [random.randint(MIN_EXON_LENGTH, MAX_EXON_LENGTH) for _ in range(num_exons)]
        total_exon_length = sum(exon_lengths)

        # Calculate minimum required intron lengths
        min_total_intron = MIN_INTRON_LENGTH * (num_exons - 1)

        # Ensure that total_exon_length + min_total_intron <= gene_length
        if total_exon_length + min_total_intron > gene_length:
            # Scale down exon lengths proportionally
            scale = (gene_length - min_total_intron) / total_exon_length
            exon_lengths = [max(int(l * scale), MIN_EXON_LENGTH) for l in exon_lengths]
            total_exon_length = sum(exon_lengths)
        
        # Calculate remaining length for introns
        remaining_length = gene_length - total_exon_length

        # Assign intron lengths ensuring minimum intron length
        if num_exons > 1:
            intron_lengths = []
            for i in range(num_exons - 1):
                # Calculate the maximum intron length possible for this intron
                # Ensure that the remaining introns can still have at least MIN_INTRON_LENGTH
                remaining_introns = (num_exons - 1) - i - 1
                max_intron = remaining_length - (MIN_INTRON_LENGTH * remaining_introns)
                if max_intron < MIN_INTRON_LENGTH:
                    # Adjust the last intron to satisfy the minimum intron length
                    max_intron = MIN_INTRON_LENGTH
                # Assign intron length randomly between MIN_INTRON_LENGTH and max_intron
                intron_length = random.randint(MIN_INTRON_LENGTH, max_intron)
                intron_lengths.append(intron_length)
                remaining_length -= intron_length
        else:
            intron_lengths = []

        # Assign exon positions within the gene
        exons = []
        exon_start = current_pos
        for exon_num, exon_length in enumerate(exon_lengths, start=1):
            exon_end = exon_start + exon_length - 1
            exons.append((exon_start, exon_end, exon_num))
            if exon_num < num_exons:
                intron_length = intron_lengths[exon_num - 1]
                exon_start = exon_end + 1 + intron_length
        
        # Update gene_end based on the last exon
        gene_end = exons[-1][1]
        
        # Final check to ensure gene_end does not exceed chromosome length
        if gene_end > CHROMOSOME_LENGTH:
            print(f"Gene{gene_id} exceeds chromosome length. Adjusting gene_end.")
            gene_end = CHROMOSOME_LENGTH
            exons[-1] = (exons[-1][0], gene_end, exons[-1][2])
        
        strand = random.choice(['+', '-'])
        
        gene = {
            'gene_id': f"Gene{gene_id}",
            'gene_name': f"Gene{gene_id}",
            'transcript_id': f"Gene{gene_id}.1",
            'start': current_pos,
            'end': gene_end,
            'strand': strand,
            'exons': exons
        }
        genes.append(gene)
        
        # Update current_pos for the next gene
        current_pos = gene_end + INTERGENIC_MIN + 1  # +1 to avoid overlap
        
        # Check if there's enough space left on the chromosome for another gene
        if current_pos + MIN_GENE_LENGTH > CHROMOSOME_LENGTH:
            print(f"Reached chromosome end. Stopping gene assignment at Gene{gene_id}.")
            break

    return genes

def create_gtf(genes, output_gtf):
    with open(output_gtf, 'w') as gtf:
        for gene in genes:
            # Write gene feature
            gtf.write(f"{CHROMOSOME_NAME}\tsource\tgene\t{gene['start']}\t{gene['end']}\t.\t{gene['strand']}\t.\tgene_id \"{gene['gene_id']}\"; gene_name \"{gene['gene_name']}\";\n")
            # Write transcript feature
            gtf.write(f"{CHROMOSOME_NAME}\tsource\ttranscript\t{gene['start']}\t{gene['end']}\t.\t{gene['strand']}\t.\tgene_id \"{gene['gene_id']}\"; transcript_id \"{gene['transcript_id']}\";\n")
            # Write exon features
            for exon in gene['exons']:
                exon_start, exon_end, exon_num = exon
                gtf.write(f"{CHROMOSOME_NAME}\tsource\texon\t{exon_start}\t{exon_end}\t.\t{gene['strand']}\t.\tgene_id \"{gene['gene_id']}\"; transcript_id \"{gene['transcript_id']}\"; exon_number \"{exon_num}\";\n")

def generate_rrna():
    num_rrna = 5  # Number of rRNA sequences
    rRNA_records = []
    for i in range(1, num_rrna + 1):
        rRNA_length = random.randint(1000, 2000)  # Length of rRNA
        rRNA_seq = generate_random_sequence(rRNA_length)
        rRNA_record = SeqRecord(Seq(rRNA_seq), id=f"rRNA{i}", description=f"Synthetic rRNA {i}")
        rRNA_records.append(rRNA_record)
    # Write rRNA.fa
    SeqIO.write(rRNA_records, "rRNA.fa", "fasta")

def main():
    # Optionally, set a random seed for reproducibility
    # random.seed(42)
    
    # Generate genome
    genome_record = generate_genome()
    # Write genome.fa
    SeqIO.write(genome_record, "genome.fa", "fasta")
    print("Generated genome.fa")
    
    # Assign genes
    genes = assign_genes()
    print(f"Assigned {len(genes)} genes")
    
    # Create genome.gtf
    create_gtf(genes, "genome.gtf")
    print("Generated genome.gtf")
    
    # Generate rRNA.fa with multiple rRNA sequences
    generate_rrna()
    print("Generated rRNA.fa")

if __name__ == "__main__":
    main()
