
![SmedAnno Logo](https://github.com/Norreanea/SmedAnno/blob/main/SmedAnno_scheme.png)

# SmedAnno

SmedAnno is a robust and versatile RNA-Seq annotation pipeline designed to streamline RNA-Seq data analysis, from preprocessing to functional annotation.

## Table of Contents

1. [Features](#features)
2. [Prerequisites](#prerequisites)
3. [Installation](#installation)
4. [Usage](#usage)
5. [Output](#output)
6. [License](#license)
7. [Acknowledgements](#acknowledgements)

## Features

- **rRNA Removal:** Efficiently removes ribosomal RNA contaminants from RNA-Seq data.
- **Read Alignment:** Aligns short and long reads to reference genomes using STAR and minimap2.
- **Gene and Transcript Assembly:** Utilizes StringTie for accurate gene and transcript assembly.
- **Assembly Merging:** Merges reference-based and de novo assemblies to create comprehensive annotations.
- **Transcript Filtering:** Filters out transcripts with excessively long exons or genomic spans to ensure data quality.
- **Isoform Comparison:** Compares assembled isoforms against reference annotations for validation.
- **GTF Correction:** Enhances GTF files with corrected and enriched annotations using AGAT.
- **Functional Annotation:** Performs functional annotation using BLASTp, BLASTx, and PFAM to assign biological functions.
- **Integration of Functional Data:** Integrates functional annotations to identify fragmented and chimeric genes.

## Prerequisites

Before installing and running SmedAnno, ensure that the following prerequisites are met:

- **Operating System:** Linux-based system recommended.
- **Conda:** Ensure that Conda is installed on your system.
- **Install System Development Libraries:**
 ```bash
sudo apt-get update
sudo apt-get install zlib1g-dev libcurl4-openssl-dev libxml2-dev
```

## Installation

### Clone the Repository

```bash
git clone https://github.com/Norreanea/SmedAnno.git
cd SmedAnno
```

### Set Up Conda Environments
SmedAnno utilizes two Conda environments for different versions of StringTie:

- **stringtie211:** Contains StringTie version 2.1.1 for steps using short reads.
- **stringtie221:** Contains StringTie version 2.2.1 for steps using mixed reads and merging annotations.

Ensure that you have added the required Conda channels (bioconda and conda-forge) with higher priority than defaults

## Usage

```bash
Usage: smedanno.sh [OPTIONS]
     
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
       6 - Filter Transcripts with Excessively Long Exons or Genomic Spans
       7 - Isoform Comparison and Annotation
       8 - GTF File Correction and Enhancement
       9 - Functional Annotation and Filtering
      10 - Integrate Functional Annotation (including Overlapped Genes and Transcripts, Reversed Duplicates, Fragmmented and Chimeric Genes Identification )
     
     Examples:
       Run all steps:
         smedanno.sh --genomeRef genome.fa --dataDir ./data --outputDir ./output --threads 4 --all
     
       Run Steps 1 and 2 only (rRNA Removal and Read Alignment):
         smedanno.sh --genomeRef genome.fa --dataDir ./data --outputDir ./output --threads 4 --steps 1,2
     
       Run Only Step 3 (Gene and Transcript Assembly) with BAM files:
         smedanno.sh --alignDir ./bam_files --outputDir ./output --threads 4 --steps 3
     
       Run Merging Assemblies Steps 4 and 5 with specific GTF directories:
         smedanno.sh --SR_RB_gtf_dir ./gtf/SR_RB --SR_DN_gtf_dir ./gtf/SR_DN --outputDir ./output --threads 4 --steps 4,5
     
       Run Functional Annotation with specific methods:
         smedanno.sh --finalGTF final_annotation.gtf --outputDir ./output --threads 4 --steps 8 --functionalMethods BLASTp,PFAM --genomeRef genome.fa
```

## Output
Upon successful execution, SmedAnno generates a structured output directory containing:

- Preprocessing Outputs: rRNA-removed reads.
- Alignment Outputs: BAM files aligned to the reference genome.
- Assembly Outputs: Assembled gene and transcript GTF files.
- Merged Annotations: Combined GTF annotations.
- Filtered Annotations: High-quality GTF files after filtering.
- Functional Annotations: BLAST and PFAM results integrated into the final annotation.
- Logs: Detailed logs for each pipeline step, facilitating troubleshooting and verification.

## License
This project is licensed under the MIT License.

## Acknowledgements
- **Bioinformatics Tools:** SmedAnno integrates several powerful tools including StringTie, STAR, minimap2, AGAT, TransDecoder, BLAST, and PFAM.
- **Open-Source Community:** Special thanks to the developers and contributors of the open-source software utilized in this pipeline.
