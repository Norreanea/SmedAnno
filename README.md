
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

- **End-to-end analysis:** Covers all stages from read preprocessing and alignment to transcript assembly and functional annotation.
- **Hybrid read support:** RNA-seq short reads (Illumina), long reads (PacBio/ONT), or a mix of both.
- **Flexible assembly:** Performs both reference-based and de novo transcript assembly with StringTie.
- **Robust annotation:** Integrates results from TransDecoder, BLAST, HMMER (Pfam), and InterProScan to build a high-quality, functionally annotated gene set.
- **Quality control:** The final step analyzes the annotation for structural and functional inconsistencies, flagging potentially fragmented or chimeric genes for manual review.
- **Reproducible environment:** Packaged in a Docker container with all dependencies managed by Conda, ensuring maximum reproducibility.

## Prerequisites

Before installing and running SmedAnno, ensure that the following prerequisites are met:

- Linux-based system recommended.
- Docker and Docker Compose must be installed. Follow the official Docker installation guide for your distribution.

## Installation
SmedAnno is designed to be run via Docker, which simplifies installation by handling all software dependencies automatically.
### Clone the Repository

```bash
git clone https://github.com/Norreanea/SmedAnno.git
cd SmedAnno
```

### Build the Docker Image
The installation is a two-step process. First, a large base image with all Conda dependencies is built. Then, the much smaller final application image is built on top of it:

```bash
# Step 1: Build the base image (can take over 20 min and requires significant disk space)
sudo docker build -t smedanno-base -f base.Dockerfile .

# Step 2: Build the final application image (this is very fast)
sudo docker build -t smedanno .
```
**Important note for Windows/WSL users (disk space)**

The smedanno-base image is very large (~136 GB) due to the extensive Conda environment. By default, Docker on WSL stores its data on your C: drive, which can cause it to fill up unexpectedly.

It is highly recommended to move Docker's data root to a larger drive (e.g., D:) before running the build command.

How to move Docker's data root in WSL:
1. Stop Docker Desktop.
2. Open a WSL terminal and create the Docker configuration directory if it doesn't exist:
```
mkdir -p ~/.docker
```
3. Create or edit the Docker daemon configuration file:
```
nano ~/.docker/daemon.json
```
4. Add the following content to the file, replacing D:\\docker-data with the path to your desired folder on another drive. Use the Windows-style path with double backslashes
```
{"data-root": "D:\\docker-data"}
```
5. Save the file

6. Restart Docker Desktop. Docker will now store all its images, containers, and volumes in the new location. You can now safely run the build commands.

## Usage
All pipeline runs are initiated through the *run_smedanno.sh* wrapper script. This script automatically handles mounting your local data directories into the Docker container:
```bash
		 Usage: ./run_smedanno.sh [OPTIONS]
		 
		 General mandatory options (if applicable based on steps selected):
		   --genomeRef PATH              Absolute path to genome reference FASTA file
		   --dataDirShort PATH           Absolute path to input directory with ONLY short-read FASTQ files.
		   --dataDirMix PATH             Absolute path to input directory with mixed-read data (must contain 'short_reads' and 'long_reads' sub-folders).
	     
		 Mandatory options for specific steps:
		   --finalGTF PATH               Absolute path to final GTF file (optional, required for Functional Annotation only)
		   --alignDirShort PATH          Absolute path to directory with ONLY short-read alignment BAM files (optional, required for Gene and Transcript Assembly only).
	       --alignDirMix PATH            Absolute path to directory with mixed-read alignment BAM files (optional, required for Gene and Transcript Assembly only).
		   --SR_RB_gtf_dir PATH          Directory containing SR_RB.gtf files (optional, required for Merging Assemblies)
		   --SR_DN_gtf_dir PATH          Directory containing SR_DN.gtf files (optional, required for Merging Assemblies)
		   --MR_RB_gtf_dir PATH          Directory containing MR_RB.gtf files (optional, required for Merging Assemblies)
		   --MR_DN_gtf_dir PATH          Directory containing MR_DN.gtf files (optional, required for Merging Assemblies)
		 
		 Long-Read trimming Options (Filtlong for Step 0):
		   --longread_min_len N          Minimum length for long reads (default: 1000)
		   --longread_keep_percent N     Keep the best N percent of reads (default: 95)
		   --longread_target_bases N     Target a total number of bases, keeping the best reads.
		   NOTE: These options provide basic filtering. For best results, use dedicated tools before using SmedAnno.
		 
		 Optional options:
		   --genomeDir PATH              Absolute path to STAR genome directory (will be created if not provided)
		   --genomeGTF PATH              Absolute path to genome annotation GTF file (optional, required for Reference-Based assembly)
		   --rrnaRef PATH                Absolute path to rRNA reference FASTA file
		   --outputDir PATH              Absolute path to output directory (default: ./outputDir)
		   --threads N                   Number of CPU threads to use (default: 8)
		   --Stringtie2minReadCoverage      Minimum read coverage allowed for the predicted transcripts (default: 1)
		   --Stringtie2minIsoformAbundance  Minimum isoform abundance of the predicted transcripts (default: 0.01)
		   --minOrfLength N              Minimum ORF length for TransDecoder (default: 100)
		   --maxExonLength N             Maximum allowed exon length (default: 10000)           
		   --maxTranscriptLength N       Maximum allowed transcript length (default: 100000)      
		   --steps LIST                  Comma-separated list of steps to run (0-10)
		   --all                         Run all steps sequentially
		   --functionalMethods METHODS   Comma-separated list of functional annotation methods to apply (BLASTp,BLASTx,PFAM,INTERPRO; default: BLASTp,BLASTx,PFAM,INTERPRO)
		   --stringtie VERSION               Set stringtie version for both short and mix reads (if not using individual overrides)
		   --stringtie_short VERSION         Set stringtie version for short reads (default: 2.1.1)
		   --stringtie_mix VERSION           Set stringtie version for mixed reads (default: 2.2.1)
		   --trimQual N                 TrimGalore quality cutoff (Phred, default 20)
		   --trimAdapter SEQ            Adapter sequence to remove               (default: auto-detection)
		   --trimGzip                   Gzip-compress trimmed FASTQ              (TRUE/FALSE)
		   --trimLen N                  Minimum read length after trimming       (default: 20 bp.)
		   --genomeType <type>           Specify the type of genome being processed. Options: 'nuclear', 'mito', 'mixed'.
		                                 If not set, the script will auto-detect 'mixed' if headers match --mitoPattern.
		                                 default: 'nuclear'.
		   --mitoPattern <regex>         Regular expression to identify mitochondrial headers for 'mixed' mode.
		                                 default: '[Mm]ito|[Mm]itochondria|mtDNA'.
		   --geneticCodeNucl <code>      Genetic code for nuclear transcripts. Default: 'Universal'.
		   --geneticCodeMito <code>      Genetic code for mitochondrial transcripts. Default: 'Mitochondrial-Vertebrates'.
		 
		_green Available genetic codes:
		   Acetabularia, Candida, Ciliate, Dasycladacean, Euplotid, Hexamita, Mesodinium,
		   Mitochondrial-Ascidian, Mitochondrial-Chlorophycean, Mitochondrial-Echinoderm, 
		   Mitochondrial-Flatworm, Mitochondrial-Invertebrates,  Mitochondrial-Protozoan,
		   Mitochondrial-Pterobranchia, Mitochondrial-Scenedesmus_obliquus, Mitochondrial-Thraustochytrium,
		   Mitochondrial-Trematode, Mitochondrial-Vertebrates, Mitochondrial-Yeast,
		   Pachysolen_tannophilus, Peritrich, SR1_Gracilibacteria, Tetrahymena, Universal
		 
		 Filtering and threshold options (for Step 10):
		   --maxDistance N               Max distance between genes to be considered fragmented (default: 1000)
		   --largeIntronThreshold N      Min size of an intron to be considered unusually large (default: 100000)
		   --blastEvalue FLOAT           E-value threshold for BLAST searches (default: 1e-5)
		   --blastIdentity N             Percent identity threshold for BLAST searches (default: 25)
		   --pfamCEvalue FLOAT           Conditional E-value threshold for Pfam (default: 1e-5)
		   --pfamBitScore N              Bit score threshold for Pfam (default: 10)
		   --pfamCoverage N              Coverage threshold for Pfam domains (default: 50)
		   --interproEvalue FLOAT        E-value threshold for InterProScan hits (default: 1e-5)
		 
		 Steps:
		   0 - Read Trimming / Quality-Filtering
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
		     ./run_smedanno.sh --genomeRef /path/to/genome.fa --dataDirShort /path/to/data --outputDir /path/to/output --threads 4 --all
		 
		   Run Steps 8, 9, 10 only starting from a GTF file:
		     ./run_smedanno.sh --finalGTF /path/to/your.gtf --genomeRef /path/to/genome.fa --outputDir /path/to/output --steps 8,9,10
```
Key points:
- Use absolute paths: Always provide absolute paths for all file and directory arguments (e.g., ```/home/user/data/genome.fa```, not ```../data/genome.fa```)
- Output directory: The ```--outputDir``` argument is always required
- Input data: Depending on the starting step, you must provide either read data (```--dataDirShort```, ```--dataDirMix```) or alignment data (```--alignDirShort```, ```--alignDirMix```)

**Example 1:** Full *de novo* run from short-read FASTQ files.
This is the most common use case. The pipeline will run all steps from 0 to 10.
```
./run_smedanno.sh \
    --all \
    --dataDirShort /path/to/your/short_reads \
    --genomeRef /path/to/your/genome.fa \
    --threads 8 \
    --outputDir /path/to/your/output_directory
```
**Example 2:** Reference-based run from existing BAM files.
Start from Step 3 using pre-aligned BAM files and a reference GTF.
```
./run_smedanno.sh \
    --steps 3,4,5,6,7,8,9,10 \
    --alignDirShort /path/to/your/bams \
    --genomeRef /path/to/your/genome.fa \
    --genomeGTF /path/to/your/reference.gtf \
    --threads 8 \
    --outputDir /path/to/your/output_directory
```
**Example 3:** Annotation-only run on a final GTF.
Run only the functional annotation and quality control steps (9 and 10) on a GTF file you've already generated.
```
./run_smedanno.sh \
    --steps 9,10 \
    --finalGTF /path/to/your/final.gtf \
    --genomeRef /path/to/your/genome.fa \
    --threads 8 \
    --outputDir /path/to/your/output_directory
```
**Full list of options.**
For a complete list of all available options and their descriptions, run:
```
./run_smedanno.sh --help
```

The pipeline is divided into the following modular steps, which can be run all at once (--all) or selectively (--steps 1,2,3...).

- Step 0: Read Trimming / Quality-Filtering 
- Step 1: rRNA Removal 
- Step 2: Read Alignment to Reference Genome 
- Step 3: Gene and Transcript Assembly
- Step 4: Merge Reference-Based Assemblies
- Step 5: Merge *De Novo* Assemblies and Create Pre-Final Annotation 
- Step 6: Filter Transcripts
- Step 7: Isoform Comparison and Annotation 
- Step 8: GTF File Correction and Enhancement 
- Step 9: Functional Annotation 
- Step 10: Integrate Functional Annotation and Quality Control 

## Output
Upon successful execution, SmedAnno generates a structured output directory containing:

- ```preprocessing/```: Cleaned FASTQ files
- ```alignment/```: BAM alignment files
- ```assembly/```: Per-sample GTF files from StringTie
- ```merging/```: Merged GTF files
- ```annotation/```: The main annotation files, including the final corrected GTF
- ```functional_annotation/```: Contains the protein sequences, BLAST/Pfam/InterProScan results, and the final, functionally annotated GTF file (```final_annotation.gtf```)
- ```logs/```: Detailed logs for every step, essential for troubleshooting


## License
This project is licensed under the MIT License.

## Acknowledgements
- **Bioinformatics tools:** SmedAnno integrates several powerful tools including Trim_Galore, Filtlong, STAR, minimap2, StringTie2, AGAT, TransDecoder, BLAST, PFAM,  InterPro and several R libraries
- **Open-source community:** Special thanks to the developers and contributors of the open-source software utilized in this pipeline
  



