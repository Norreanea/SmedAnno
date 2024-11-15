# Genome-annotation-pipeline

Conda environments:
- **stringtie211:** Contains StringTie version 2.1.1 for steps using short reads.
- **stringtie221:** Contains StringTie version 2.2.1 for steps using mixed reads and merging annotations.

Please make sure that Conda is installed 

Make the script executable:
```bash
chmod +x rna_seq_pipeline.sh
```
Required inputs:
- **rRNA Reference FASTA:** Path to rrna.fa.
- **Genome reference FASTA:** Path to genome.fa.
- **Genome rnnotation GTF:** Path to reference.gtf.
- **BLAST databases:** Ensure BLAST NR and SwissProt databases are downloaded and formatted.
- **PFAM database:** Download and prepare the PFAM database for pfam_scan.pl.
- **rRNA STAR Index:** Ensure you have a STAR index for rRNA references.
- **Genome STAR Index:** Ensure you have a STAR index for the genome.

```bash
Usage: ./rna_seq_pipeline.bash [OPTIONS]

Mandatory Options for Specific Steps:
  --finalGTF PATH               Path to final GTF file (required for Functional Annotation only)
  --alignDir PATH               Path to directory containing BAM files (required for Gene and Transcript Assembly only)
  --SR_RB_gtf_dir PATH          Directory containing SR_RB.gtf files (required for Merging Assemblies)
  --SR_DN_gtf_dir PATH          Directory containing SR_DN.gtf files (optional for Merging Assemblies)
  --MR_RB_gtf_dir PATH          Directory containing MR_RB.gtf files (optional for Merging Assemblies)
  --MR_DN_gtf_dir PATH          Directory containing MR_DN.gtf files (optional for Merging Assemblies)

General Mandatory Options (if applicable based on steps selected):
  --genomeRef PATH              Path to genome reference FASTA file
  --dataDir PATH                Path to input RNA-Seq data directory (must contain short_reads and/or mix_reads folders)

Optional Options:
  --genomeDir PATH              Path to STAR genome directory (will be created if not provided)
  --genomeGTF PATH              Path to genome annotation GTF file (optional, required for Reference-Based assembly)
  --rrnaRef PATH                Path to rRNA reference FASTA file
  --blastDB_SwissProt PATH      Path to BLAST SwissProt database
  --pfamDB PATH                 Path to PFAM database directory
  --outputDir PATH              Path to output directory (default: ./outputDir)
  --threads N                   Number of CPU threads to use (default: 8)
  --minOrfLength N              Minimum ORF length for TransDecoder (default: 100)
  --steps LIST                  Comma-separated list of steps to run (1-8, include 5.1)
  --all                         Run all steps sequentially
  --functionalMethods METHODS    Comma-separated list of functional annotation methods to apply (BLASTp,BLASTx,PFAM; default: all)
  --help                        Display this help message and exit

Steps:
  1 - rRNA Removal
  2 - Read Alignment to Reference Genome
  3 - Gene and Transcript Assembly
  4 - Merge Reference-Based Assemblies
  5 - Merge De Novo Assemblies and Create Pre-Final Annotation
  5.1 - Filter Transcripts with Excessively Long Exons or Genomic Spans
  6 - Isoform Comparison and Annotation
  7 - GTF File Correction and Enhancement
  8 - Functional Annotation and Filtering

Examples:
  Run all steps:
    ./rna_seq_pipeline.bash --genomeRef genome.fa --dataDir ./data --outputDir ./output --threads 4 --all

  Run Steps 1 and 2 only (rRNA Removal and Read Alignment):
    ./rna_seq_pipeline.bash --genomeRef genome.fa --dataDir ./data --outputDir ./output --threads 4 --steps 1,2

  Run Only Step 3 (Gene and Transcript Assembly) with BAM files:
    ./rna_seq_pipeline.bash --alignDir ./bam_files --outputDir ./output --threads 4 --steps 3

  Run Merging Assemblies Steps 4 and 5 with specific GTF directories:
    ./rna_seq_pipeline.bash --SR_RB_gtf_dir ./gtf/SR_RB --SR_DN_gtf_dir ./gtf/SR_DN --outputDir ./output --threads 4 --steps 4,5

  Run Functional Annotation with specific methods:
    ./rna_seq_pipeline.bash --finalGTF final_annotation.gtf --outputDir ./output --threads 4 --steps 8 --functionalMethods BLASTp,PFAM --genomeRef genome.fa
```
