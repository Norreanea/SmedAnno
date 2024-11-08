#!/bin/bash

# =====================================================================
# Comprehensive RNA-Seq Bioinformatics Pipeline
# =====================================================================

set -e  # Exit immediately if a command exits with a non-zero status
set -o pipefail  # Pipeline returns the exit status of the last command to fail
#set -x  # Enable debugging
# ---------------------------
# Create Linux-Native Temporary Directory for STAR
# ---------------------------
mkdir -p /tmp/temp_star
chmod 777 /tmp/temp_star

# ---------------------------
# Define Color Variables for Echo Messages
# ---------------------------
RED='\033[0;31m'
GREEN='\033[0;32m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# ---------------------------
# Helper Functions for Colored Echo
# ---------------------------
echo_red() {
    echo -e "${RED}$1${NC}"
}

echo_green() {
    echo -e "${GREEN}$1${NC}"
}

echo_blue() {
    echo -e "${BLUE}$1${NC}"
}
# ---------------------------
# Function to Display Help
# ---------------------------
show_help() {
    echo_blue "Usage: $0 [OPTIONS]"
    echo ""
    echo_blue "Mandatory Options:"
    echo "  --genomeRef PATH          Path to genome reference FASTA file"
    echo "  --dataDir PATH            Path to input RNA-Seq data directory (must contain short_reads and/or mix_reads folders)"
    echo ""
    echo "Optional Options:"
    echo "  --genomeDir PATH          Path to STAR genome directory (will be created if not provided)"
    echo "  --genomeGTF PATH          Path to genome annotation GTF file (optional, required for Reference-Based assembly)"
    echo "  --rrnaRef PATH            Path to rRNA reference FASTA file"
    echo "  --blastDB_SwissProt PATH  Path to BLAST SwissProt database"
    echo "  --pfamDB PATH             Path to PFAM database directory"
    echo "  --outputDir PATH          Path to output directory (default: ./outputDir)"
    echo "  --threads N               Number of CPU threads to use (default: 8)"
    echo "  --minOrfLength N          Minimum ORF length for TransDecoder (default: 100)"
    echo "  --steps LIST              Comma-separated list of steps to run (1-9)"
    echo "  --all                     Run all steps"
    echo "  --help                    Display this help message and exit"
    echo ""
    echo_blue "Steps:"
    echo "  1 - Genome Index Generation"
    echo "  2 - rRNA Removal"
    echo "  3 - Read Alignment to Reference Genome"
    echo "  4 - Gene and Transcript Assembly"
    echo "  5 - Merge Reference-Based Assemblies"
    echo "  6 - Merge De Novo Assemblies and Create Pre-Final Annotation"
    echo "  7 - Isoform Comparison and Annotation"
    echo "  8 - GTF File Correction and Enhancement"
    echo "  9 - Functional Annotation and Filtering"
    echo ""
    echo "Example:"
    echo "  $0 --genomeRef genome.fa --dataDir ./data --outputDir ./output --threads 4 --all"
    exit 1
}

# ---------------------------
# Function to Determine if a File is Gzipped
# ---------------------------
is_gzipped() {
    local file="$1"
    if [[ "${file}" == *.gz ]]; then
        echo "yes"
    else
        echo "no"
    fi
}

# ---------------------------
# Helper Function to Find Read Files (Handles .fq, .fq.gz, .fastq, .fastq.gz)
# ---------------------------
find_read_file() {
    local dir="$1"
    local sample="$2"
    local read_num="$3"
    local extensions=("fq.gz" "fq" "fastq.gz" "fastq")
    
    for ext in "${extensions[@]}"; do
        local file="${dir}/${sample}${read_num}.${ext}"
        if [ -f "${file}" ]; then
            echo "${file}"
            return
        fi
    done
    echo ""
}

# ---------------------------
# Function to Check Global Dependencies
# ---------------------------
check_dependencies() {
    local dependencies=("conda" "wget" "python3")
    for cmd in "${dependencies[@]}"; do
        if ! command -v "$cmd" &> /dev/null; then
            echo_red  "Error: Required tool '$cmd' is not installed or not in PATH."
            echo_red  "Please install '$cmd' before running the pipeline."
            exit 1
        fi
    done
}

# ---------------------------
# Function to Create Conda Environments and Install Tools
# ---------------------------
setup_conda_environments() {
    echo_green  "Setting up conda environments"

    # Initialize Conda
    # shellcheck source=/dev/null
    source "$(conda info --base)/etc/profile.d/conda.sh" || { echo_red  "Error: Failed to source conda.sh"; exit 1; }

    # Create conda environment for StringTie v2.1.1
    if ! conda env list | grep -q 'stringtie211'; then
        echo_green  "Creating conda environment: stringtie211"
        if ! conda create -y -n stringtie211 stringtie=2.1.1 &> "${OUTPUT_DIR}/conda_create_stringtie211.log"; then
            echo_red  "Error: Failed to create stringtie211 environment. Check ${OUTPUT_DIR}/conda_create_stringtie211.log for details."
            exit 1
        fi
    else
        echo_green  "Conda environment 'stringtie211' already exists."
    fi

    # Create conda environment for StringTie v2.2.1 and other tools
    if ! conda env list | grep -q 'stringtie221'; then
        echo_green  "Creating conda environment: stringtie221"
        if ! conda create -y -n stringtie221 stringtie=2.2.1 gffcompare agat TransDecoder blast pfam_scan gffread samtools star minimap2 &> "${OUTPUT_DIR}/conda_create_stringtie221.log"; then
            echo_red  "Error: Failed to create stringtie221 environment. Check ${OUTPUT_DIR}/conda_create_stringtie221.log for details."
            exit 1
        fi
    else
        echo_green  "Conda environment 'stringtie221' already exists."
    fi

    echo_green  "Conda environments set up successfully."
}

# ---------------------------
# Function for Step 1: Genome Index Generation
# ---------------------------
step1_genome_index() {
    echo_blue  "Starting Step 1: Genome Index Generation"
    conda activate stringtie221

    if [ -d "${GENOME_DIR}" ] && [ -f "${GENOME_DIR}/SA" ]; then
        echo_green  "Genome index already exists at ${GENOME_DIR}, skipping index generation."
    else
        echo_green  "Calculating genome length for genomeSAindexNbases parameter"
        
        # Calculate genome length
        GENOME_LENGTH=$(awk '/^>/ {next} {len+=length($0)} END {print len}' "${GENOME_REF}")

        # Calculate genomeSAindexNbases
        SA_INDEX_NBASES=$(python3 -c "import math; l=math.log(${GENOME_LENGTH},2); n=int(l/2 -1); print(min(14, n))")

        echo_green  "Genome length: ${GENOME_LENGTH}"
        echo_green  "Calculated genomeSAindexNbases: ${SA_INDEX_NBASES}"

        echo_green  "Generating genome index with STAR at ${GENOME_DIR}"
        mkdir -p "${GENOME_DIR}"
        STAR --runThreadN "${THREADS}" \
             --runMode genomeGenerate \
             --genomeDir "${GENOME_DIR}" \
             --genomeFastaFiles "${GENOME_REF}" \
             --genomeSAindexNbases "${SA_INDEX_NBASES}" \
             --outTmpDir /tmp/temp_star/temp

        echo_green  "Genome index generation completed."
    fi

    echo_blue  "Step 1 Completed: Genome Index Generation"
    conda deactivate
}

# ---------------------------
# Function for Step 2: Preprocessing - rRNA Removal
# ---------------------------
step2_rrna_removal() {
    echo_blue  "Starting Step 2: rRNA Removal"
    conda activate stringtie221

    if [ -z "${RRNA_REF}" ]; then
        echo_green  "No rRNA reference provided. Skipping rRNA removal."
        conda deactivate
        return
    fi

    # Calculate rRNA reference length for genomeSAindexNbases parameter
    echo_green  "Calculating rRNA reference length for genomeSAindexNbases parameter"
    
    # Calculate rRNA reference length
    RRNA_LENGTH=$(awk '/^>/ {next} {len+=length($0)} END {print len}' "${RRNA_REF}")

    # Calculate genomeSAindexNbases for rRNA reference
    RRNA_SA_INDEX_NBASES=$(python3 -c "import math; l=math.log(${RRNA_LENGTH},2); n=int(l/2 -1); print(min(14, n))")

    echo_green  "rRNA reference length: ${RRNA_LENGTH}"
    echo_green  "Calculated genomeSAindexNbases for rRNA reference: ${RRNA_SA_INDEX_NBASES}"

    # Build STAR index for rRNA reference if not already done
    if [ ! -d "${RRNA_REF_INDEX}" ] || [ ! -f "${RRNA_REF_INDEX}/SA" ]; then
        echo_green  "Building STAR index for rRNA reference at ${RRNA_REF_INDEX}"
        mkdir -p "${RRNA_REF_INDEX}"
        STAR --runThreadN "${THREADS}" \
             --runMode genomeGenerate \
             --genomeDir "${RRNA_REF_INDEX}" \
             --genomeFastaFiles "${RRNA_REF}" \
             --genomeSAindexNbases "${RRNA_SA_INDEX_NBASES}" \
             --outTmpDir /tmp/temp_star/temp
    else
        echo_green  "rRNA STAR index already exists at ${RRNA_REF_INDEX}, skipping index generation."
    fi

    # Process short-only samples
    if [ "${#SHORT_ONLY_SAMPLES[@]}" -gt 0 ]; then
        for SAMPLE in "${SHORT_ONLY_SAMPLES[@]}"; do
            echo_green  "Processing Short-Only Sample: ${SAMPLE}"

            # Find read1 and read2 files
            READ1=$(find_read_file "${SHORT_ONLY_DIR}" "${SAMPLE}" "1")
            READ2=$(find_read_file "${SHORT_ONLY_DIR}" "${SAMPLE}" "2")

            # Check if read files are found
            if [ -z "${READ1}" ] || [ -z "${READ2}" ]; then
                echo_red "Error: One or both short read files for sample ${SAMPLE} not found."
                exit 1
            fi

            # Determine if files are gzipped
            READ1_GZ=$(is_gzipped "${READ1}")
            READ2_GZ=$(is_gzipped "${READ2}")

            # Set readFilesCommand based on compression
            if [ "${READ1_GZ}" == "yes" ] && [ "${READ2_GZ}" == "yes" ]; then
                READ_FILES_CMD="zcat"
            else
                READ_FILES_CMD="cat"
            fi

            # Output files
            NON_RRNA_READ1="${PREPROC_DIR}/${SAMPLE}_non_rrna_R1.fq.gz"
            NON_RRNA_READ2="${PREPROC_DIR}/${SAMPLE}_non_rrna_R2.fq.gz"

            # Remove rRNA using STAR
            echo_green  "Removing rRNA from short reads using STAR for Sample: ${SAMPLE}"
            STAR --runThreadN "${THREADS}" \
                 --genomeDir "${RRNA_REF_INDEX}" \
                 --readFilesIn "${READ1}" "${READ2}" \
                 --readFilesCommand "${READ_FILES_CMD}" \
                 --outFilterType BySJout \
                 --outFilterMultimapNmax 1 \
                 --outFilterMismatchNmax 0 \
                 --outReadsUnmapped Fastx \
                 --outSAMtype BAM SortedByCoordinate \
                 --outFileNamePrefix "${PREPROC_DIR}/${SAMPLE}_STAR_rrna_" \
                 --outTmpDir /tmp/temp_star/temp

            # Extract non-rRNA reads
            gzip "${PREPROC_DIR}/${SAMPLE}_STAR_rrna_Unmapped.out.mate1"
            gzip "${PREPROC_DIR}/${SAMPLE}_STAR_rrna_Unmapped.out.mate2"
            mv "${PREPROC_DIR}/${SAMPLE}_STAR_rrna_Unmapped.out.mate1.gz" "${NON_RRNA_READ1}"
            mv "${PREPROC_DIR}/${SAMPLE}_STAR_rrna_Unmapped.out.mate2.gz" "${NON_RRNA_READ2}"

            echo_green  "rRNA removal completed for Short-Only Sample: ${SAMPLE}"
        done
    fi

    # Process mixed samples
    if [ "${#MIX_SAMPLES[@]}" -gt 0 ]; then
        for SAMPLE in "${MIX_SAMPLES[@]}"; do
            echo "Processing Mixed Sample: ${SAMPLE}"

            # Find mixed short read1 and read2 files
            MIX_READ1=$(find_read_file "${MIX_SHORT_DIR}" "${SAMPLE}" "1")
            MIX_READ2=$(find_read_file "${MIX_SHORT_DIR}" "${SAMPLE}" "2")

            # Check if short read files are found
            if [ -z "${MIX_READ1}" ] || [ -z "${MIX_READ2}" ]; then
                echo_red  "Error: Short read files for mixed sample ${SAMPLE} not found."
                exit 1
            fi

            # Determine if short read files are gzipped
            MIX_READ1_GZ=$(is_gzipped "${MIX_READ1}")
            MIX_READ2_GZ=$(is_gzipped "${MIX_READ2}")

            # Set readFilesCommand for short reads
            if [ "${MIX_READ1_GZ}" == "yes" ] && [ "${MIX_READ2_GZ}" == "yes" ]; then
                MIX_READ_FILES_CMD="zcat"
            else
                MIX_READ_FILES_CMD="cat"
            fi

            # Output files for short reads
            MIX_NON_RRNA_READ1="${PREPROC_DIR}/${SAMPLE}_non_rrna_R1.fq.gz"
            MIX_NON_RRNA_READ2="${PREPROC_DIR}/${SAMPLE}_non_rrna_R2.fq.gz"

            # Remove rRNA from short reads using STAR
            echo_green  "Removing rRNA from mixed sample's short reads using STAR for Sample: ${SAMPLE}"
            STAR --runThreadN "${THREADS}" \
                 --genomeDir "${RRNA_REF_INDEX}" \
                 --readFilesIn "${MIX_READ1}" "${MIX_READ2}" \
                 --readFilesCommand "${MIX_READ_FILES_CMD}" \
                 --outFilterType BySJout \
                 --outFilterMultimapNmax 1 \
                 --outFilterMismatchNmax 0 \
                 --outReadsUnmapped Fastx \
                 --outSAMtype BAM SortedByCoordinate \
                 --outFileNamePrefix "${PREPROC_DIR}/${SAMPLE}_STAR_rrna_" \
                 --outTmpDir /tmp/temp_star/temp

            # Extract non-rRNA reads for short reads
            gzip "${PREPROC_DIR}/${SAMPLE}_STAR_rrna_Unmapped.out.mate1"
            gzip "${PREPROC_DIR}/${SAMPLE}_STAR_rrna_Unmapped.out.mate2"
            mv "${PREPROC_DIR}/${SAMPLE}_STAR_rrna_Unmapped.out.mate1.gz" "${MIX_NON_RRNA_READ1}"
            mv "${PREPROC_DIR}/${SAMPLE}_STAR_rrna_Unmapped.out.mate2.gz" "${MIX_NON_RRNA_READ2}"

            echo_green  "rRNA removal completed for short reads of Mixed Sample: ${SAMPLE}"

            # Process long reads
            # Handle both .fq and .fastq extensions
            MIX_LONG_READ=$(find_read_file "${MIX_LONG_DIR}" "${SAMPLE}" "_long_reads")

            if [ -z "${MIX_LONG_READ}" ]; then
                echo_red  "No long reads found for Mixed Sample: ${SAMPLE}"
                continue
            fi

            MIX_NON_RRNA_LONG_READ="${PREPROC_DIR}/${SAMPLE}_non_rrna_long_reads.fq.gz"
            MIX_LONG_READ_EXTENSION="${MIX_LONG_READ##*.}"  # Get extension

            echo_green  "Removing rRNA from mixed sample's long reads using minimap2 for Sample: ${SAMPLE}"
            
            # Remove rRNA from long reads using minimap2
            if [[ "${MIX_LONG_READ}" == *.gz ]]; then
                minimap2 -t "${THREADS}" -ax map-pb "${RRNA_REF}" <(zcat "${MIX_LONG_READ}") | \
                    samtools view -Sb - > "${ALIGN_DIR}/${SAMPLE}_long_aligned.bam"
            else
                minimap2 -t "${THREADS}" -ax map-pb "${RRNA_REF}" "${MIX_LONG_READ}" | \
                    samtools view -Sb - > "${ALIGN_DIR}/${SAMPLE}_long_aligned.bam"
            fi

            # **New Sorting Step for Long Reads BAM**
            if [ -f "${ALIGN_DIR}/${SAMPLE}_long_aligned.bam" ]; then
                echo_green  "Sorting long reads BAM file for Sample: ${SAMPLE}"
                samtools sort -@ "${THREADS}" -o "${ALIGN_DIR}/${SAMPLE}_long_aligned.sorted.bam" "${ALIGN_DIR}/${SAMPLE}_long_aligned.bam"
                rm "${ALIGN_DIR}/${SAMPLE}_long_aligned.bam"  # Remove unsorted BAM
                mv "${ALIGN_DIR}/${SAMPLE}_long_aligned.sorted.bam" "${ALIGN_DIR}/${SAMPLE}_long_aligned.bam"  # Rename to original name
                echo_green  "Long reads BAM sorted for Sample: ${SAMPLE}"
            fi
			samtools view -b -f 4 "${ALIGN_DIR}/${SAMPLE}_long_aligned.bam" | samtools sort -n | samtools fastq | gzip > "${MIX_NON_RRNA_LONG_READ}"

			#rm "${ALIGN_DIR}/${SAMPLE}_long_aligned.sorted.bam"


            echo_green  "rRNA removal completed for long reads of Mixed Sample: ${SAMPLE}"
            echo_green  "rRNA removal completed for Mixed Sample: ${SAMPLE}"
        done
    fi

    echo_blue  "Step 2 Completed: rRNA Removal"
    conda deactivate 
}

# ---------------------------
# Function for Step 3: Read Alignment to Reference Genome
# ---------------------------
step3_read_alignment() {
    echo_blue  "Starting Step 3: Read Alignment to the Reference Genome"
    conda activate stringtie221

    # Process short-only samples
    if [ "${#SHORT_ONLY_SAMPLES[@]}" -gt 0 ]; then
        for SAMPLE in "${SHORT_ONLY_SAMPLES[@]}"; do
            echo_green  "Aligning Short-Only Sample: ${SAMPLE}"

            # Determine if rRNA removal was performed
            if [ -n "${RRNA_REF}" ]; then
                READ1="${PREPROC_DIR}/${SAMPLE}_non_rrna_R1.fq.gz"
                READ2="${PREPROC_DIR}/${SAMPLE}_non_rrna_R2.fq.gz"
            else
                READ1=$(find_read_file "${SHORT_ONLY_DIR}" "${SAMPLE}" "1")
                READ2=$(find_read_file "${SHORT_ONLY_DIR}" "${SAMPLE}" "2")
            fi

            # Check if read files are found
            if [ -z "${READ1}" ] || [ -z "${READ2}" ]; then
                echo_red  "Error: Read files for sample ${SAMPLE} not found."
                exit 1
            fi

            # Determine if files are gzipped
            READ1_GZ=$(is_gzipped "${READ1}")
            READ2_GZ=$(is_gzipped "${READ2}")

            # Set readFilesCommand based on compression
            if [ "${READ1_GZ}" == "yes" ] && [ "${READ2_GZ}" == "yes" ]; then
                READ_FILES_CMD="zcat"
            else
                READ_FILES_CMD="cat"
            fi

            # Align reads using STAR (sorted BAM)
            echo_green  "Running STAR aligner for Short-Only Sample: ${SAMPLE}"
            STAR --runThreadN "${THREADS}" \
                 --genomeDir "${GENOME_DIR}" \
                 --readFilesIn "${READ1}" "${READ2}" \
                 --readFilesCommand "${READ_FILES_CMD}" \
                 --outSAMtype BAM SortedByCoordinate \
                 --outFileNamePrefix "${ALIGN_DIR}/${SAMPLE}_STAR_" \
                 --outTmpDir /tmp/temp_star/temp

            # Define paths
            ALIGNED_SORTED_BAM="${ALIGN_DIR}/${SAMPLE}_STAR_Aligned.sortedByCoord.out.bam"

            # Check if STAR produced the BAM file
            if [ ! -f "${ALIGNED_SORTED_BAM}" ]; then
                echo_red  "Error: STAR did not produce the expected sorted BAM file for sample ${SAMPLE}."
                exit 1
            fi

            # Add XS attribute
            echo_green  "Adding XS attribute to BAM file for Sample: ${SAMPLE}"
            samtools view -h "${ALIGNED_SORTED_BAM}" | \
                awk -v strType=2 -f ./tagXSstrandedData.awk | \
                samtools view -bSq 10 -F 4 - > "${ALIGN_DIR}/${SAMPLE}_STAR_Aligned.sortedByCoord.outXS.bam"

            echo_green  "Alignment completed for Short-Only Sample: ${SAMPLE}"
        done
    fi

    # Process mixed samples
    if [ "${#MIX_SAMPLES[@]}" -gt 0 ]; then
        for SAMPLE in "${MIX_SAMPLES[@]}"; do
            echo_green  "Aligning Mixed Sample: ${SAMPLE}"

            # Determine if rRNA removal was performed
            if [ -n "${RRNA_REF}" ]; then
                MIX_READ1="${PREPROC_DIR}/${SAMPLE}_non_rrna_R1.fq.gz"
                MIX_READ2="${PREPROC_DIR}/${SAMPLE}_non_rrna_R2.fq.gz"
                MIX_NON_RRNA_LONG_READ="${PREPROC_DIR}/${SAMPLE}_non_rrna_long_reads.fq.gz"
            else
                MIX_READ1=$(find_read_file "${MIX_SHORT_DIR}" "${SAMPLE}" "1")
                MIX_READ2=$(find_read_file "${MIX_SHORT_DIR}" "${SAMPLE}" "2")
                MIX_LONG_READ=$(find_read_file "${MIX_LONG_DIR}" "${SAMPLE}" "_long_reads")
            fi

            # Check if short read files are found
            if [ -z "${MIX_READ1}" ] || [ -z "${MIX_READ2}" ]; then
                echo_red  "Error: Short read files for mixed sample ${SAMPLE} not found."
                exit 1
            fi

            # Determine if short read files are gzipped
            MIX_READ1_GZ=$(is_gzipped "${MIX_READ1}")
            MIX_READ2_GZ=$(is_gzipped "${MIX_READ2}")

            # Set readFilesCommand for short reads
            if [ "${MIX_READ1_GZ}" == "yes" ] && [ "${MIX_READ2_GZ}" == "yes" ]; then
                MIX_READ_FILES_CMD="zcat"
            else
                MIX_READ_FILES_CMD="cat"
            fi

            # Align short reads using STAR (sorted BAM)
            echo_green  "Running STAR aligner for Mixed Sample: ${SAMPLE} (Short Reads)"
            STAR --runThreadN "${THREADS}" \
                 --genomeDir "${GENOME_DIR}" \
                 --readFilesIn "${MIX_READ1}" "${MIX_READ2}" \
                 --readFilesCommand "${MIX_READ_FILES_CMD}" \
                 --outSAMtype BAM SortedByCoordinate \
                 --outFileNamePrefix "${ALIGN_DIR}/${SAMPLE}_STAR_" \
                 --outTmpDir /tmp/temp_star/temp

            # Define paths
            ALIGNED_SHORT_SORTED_BAM="${ALIGN_DIR}/${SAMPLE}_STAR_Aligned.sortedByCoord.out.bam"

            # Check if STAR produced the BAM file
            if [ ! -f "${ALIGNED_SHORT_SORTED_BAM}" ]; then
                echo_red  "Error: STAR did not produce the expected sorted BAM file for sample ${SAMPLE}."
                exit 1
            fi

            # Add XS attribute
            echo_green  "Adding XS attribute to short BAM file for Sample: ${SAMPLE}"
            samtools view -h "${ALIGNED_SHORT_SORTED_BAM}" | \
                awk -v strType=2 -f ./tagXSstrandedData.awk | \
                samtools view -bSq 10 -F 4 - > "${ALIGN_DIR}/${SAMPLE}_STAR_Aligned.sortedByCoord.outXS.bam"

            # Align long reads using minimap2
            if [ -n "${RRNA_REF}" ] && [ -f "${MIX_NON_RRNA_LONG_READ}" ]; then
                echo_green  "Running minimap2 aligner for Mixed Sample: ${SAMPLE} (Long Reads)"
                minimap2 -t "${THREADS}" -ax map-pb "${GENOME_REF}" "${MIX_NON_RRNA_LONG_READ}" | \
                    samtools view -Sb - > "${ALIGN_DIR}/${SAMPLE}_long_aligned.bam"
            elif [ -z "${RRNA_REF}" ] && [ -f "${MIX_LONG_READ}" ]; then
                echo_green  "Running minimap2 aligner for Mixed Sample: ${SAMPLE} (Long Reads)"
                minimap2 -t "${THREADS}" -ax map-pb "${GENOME_REF}" "${MIX_LONG_READ}" | \
                    samtools view -Sb - > "${ALIGN_DIR}/${SAMPLE}_long_aligned.bam"
            else
                echo_green  "No long reads found for Mixed Sample: ${SAMPLE}"
            fi

            # **New Sorting Step for Long Reads BAM**
            if [ -f "${ALIGN_DIR}/${SAMPLE}_long_aligned.bam" ]; then
                echo_green  "Sorting long reads BAM file for Sample: ${SAMPLE}"
                samtools sort -@ "${THREADS}" -o "${ALIGN_DIR}/${SAMPLE}_long_aligned.sorted.bam" "${ALIGN_DIR}/${SAMPLE}_long_aligned.bam"
                rm "${ALIGN_DIR}/${SAMPLE}_long_aligned.bam"  # Remove unsorted BAM
                mv "${ALIGN_DIR}/${SAMPLE}_long_aligned.sorted.bam" "${ALIGN_DIR}/${SAMPLE}_long_aligned.bam"  # Rename to original name
                echo_green  "Long reads BAM sorted for Sample: ${SAMPLE}"
            fi

            echo_green  "Alignment completed for Mixed Sample: ${SAMPLE}"
        done
    fi

    echo_blue  "Step 3 Completed: Read Alignment"
    conda deactivate 
}

# ---------------------------
# Function for Step 4: Gene and Transcript Assembly
# ---------------------------
step4_gene_transcript_assembly() {
    echo_blue  "Starting Step 4: Gene and Transcript Assembly"

    # Process short-only samples
    if [ "${#SHORT_ONLY_SAMPLES[@]}" -gt 0 ]; then
        for SAMPLE in "${SHORT_ONLY_SAMPLES[@]}"; do
            echo_green  "Assembling transcripts for Short-Only Sample: ${SAMPLE}"

            ALIGNED_SORTED_BAM="${ALIGN_DIR}/${SAMPLE}_STAR_Aligned.sortedByCoord.outXS.bam"

            # Reference-Based Annotation (RB) with StringTie2 v2.1.1 (short reads)
            if [ -n "${GENOME_GTF}" ]; then
                echo_green  "Running StringTie2 v2.1.1 for Reference-Based Assembly (RB, SR)"
                conda activate stringtie211
                stringtie -p "${THREADS}" \
                          -G "${GENOME_GTF}" \
                          -c 1.5 \
                          -f 0.02 \
                          -o "${ASSEMBLY_DIR}/${SAMPLE}_SR_RB.gtf" \
                          "${ALIGNED_SORTED_BAM}"
                conda deactivate
            fi

            # De Novo Annotation (DN) with StringTie2 v2.1.1 (short reads)
            echo_green  "Running StringTie2 v2.1.1 for De Novo Assembly (DN, SR)"
            conda activate stringtie211
            stringtie -p "${THREADS}" \
                      -c 1.5 \
                      -f 0.02 \
                      -o "${ASSEMBLY_DIR}/${SAMPLE}_SR_DN.gtf" \
                      "${ALIGNED_SORTED_BAM}"
            conda deactivate

            echo_green  "Transcript assembly completed for Short-Only Sample: ${SAMPLE}"
        done
    fi

    # Process mixed samples
    if [ "${#MIX_SAMPLES[@]}" -gt 0 ]; then
        for SAMPLE in "${MIX_SAMPLES[@]}"; do
            echo_green  "Assembling transcripts for Mixed Sample: ${SAMPLE}"

            ALIGNED_SHORT_SORTED_BAM="${ALIGN_DIR}/${SAMPLE}_STAR_Aligned.sortedByCoord.outXS.bam"
            ALIGNED_LONG_BAM="${ALIGN_DIR}/${SAMPLE}_long_aligned.bam"  # This BAM is now sorted
			
			echo_green  "Short reads BAM: ${ALIGNED_SHORT_SORTED_BAM}"
			echo_green  "Long reads BAM: ${ALIGNED_LONG_BAM}"

			if [ ! -f "${ALIGNED_SHORT_SORTED_BAM}" ]; then
				echo_red  "Error: Short reads BAM file ${ALIGNED_SHORT_SORTED_BAM} does not exist."
				exit 1
			fi

            # Check if long reads alignment exists
            if [ ! -f "${ALIGNED_LONG_BAM}" ]; then
                echo_red  "No aligned long reads BAM found for Sample: ${SAMPLE}, skipping mixed assembly."
                continue
            fi

            # Reference-Based Annotation (RB) with StringTie2 v2.2.1 (Mixed Reads)
            if [ -n "${GENOME_GTF}" ]; then
                echo_green  "Running StringTie2 v2.2.1 for Reference-Based Assembly (RB, MR)"
                conda activate stringtie221
                stringtie --mix -p "${THREADS}" \
                          -G "${GENOME_GTF}" \
                          -c 1.5 \
                          -f 0.02 \
                          -o "${ASSEMBLY_DIR}/${SAMPLE}_MR_RB.gtf" \
                          "${ALIGNED_SHORT_SORTED_BAM}" "${ALIGNED_LONG_BAM}"
                conda deactivate
            fi

            # De Novo Annotation (DN) with StringTie2 v2.2.1 (Mixed Reads)
            echo_green  "Running StringTie2 v2.2.1 for De Novo Assembly (DN, MR)"
            conda activate stringtie221
            stringtie --mix -p "${THREADS}" \
                      -c 1.5 \
                      -f 0.02 \
                      -o "${ASSEMBLY_DIR}/${SAMPLE}_MR_DN.gtf" \
                      "${ALIGNED_SHORT_SORTED_BAM}" "${ALIGNED_LONG_BAM}"
            conda deactivate

            echo_green  "Transcript assembly completed for Mixed Sample: ${SAMPLE}"
        done
    fi

    echo_blue  "Step 4 Completed: Gene and Transcript Assembly"
}

# ---------------------------
# Function for Step 5: Merge Reference-Based Assemblies
# ---------------------------
step5_merging_transcripts_I() {
    echo_blue "Starting Step 5: Merge Reference-Based Assemblies"

    if [ -z "${GENOME_GTF}" ]; then
        echo_red "No reference GTF provided. Skipping Reference-Based merging."
        return
    fi

    # Activate StringTie 2.2.1 environment
    conda activate stringtie221

    # Merge Reference-Based Assemblies for Short Reads (SR)
    if [ "${#SHORT_ONLY_SAMPLES[@]}" -gt 0 ]; then
        echo_green "Merging Reference-Based Assemblies for Short Reads (SR)"
        LIST_SR_RB_GTF="${MERGE_DIR}/SR_RB_gtf_list.txt"
        > "${LIST_SR_RB_GTF}"  # Empty the file

        for SAMPLE in "${SHORT_ONLY_SAMPLES[@]}"; do
            GTF_FILE="${ASSEMBLY_DIR}/${SAMPLE}_SR_RB.gtf"
            if [ -f "${GTF_FILE}" ] && [ -s "${GTF_FILE}" ] && grep -qv "^#" "${GTF_FILE}"; then
                echo "${GTF_FILE}" >> "${LIST_SR_RB_GTF}"
            else
                echo_red "Warning: GTF file ${GTF_FILE} is missing, empty, or contains no transcripts. Skipping."
            fi
        done

        if [ -s "${LIST_SR_RB_GTF}" ]; then
            stringtie --merge -p "${THREADS}" \
                      -G "${GENOME_GTF}" \
                      -o "${MERGE_DIR}/merged_SR_RB.gtf" \
                      "${LIST_SR_RB_GTF}"
            echo_green "Merged SR_RB.gtf created successfully."
        else
            echo_red "No valid SR_RB GTF files available for merging. Skipping."
        fi
    fi

    # Merge Reference-Based Assemblies for Mixed Reads (MR) with SR
    if [ "${#MIX_SAMPLES[@]}" -gt 0 ]; then
        echo_green "Merging Reference-Based Assemblies for Mixed Reads (MR) with SR"
        LIST_MR_RB_GTF="${MERGE_DIR}/MR_RB_gtf_list.txt"
        > "${LIST_MR_RB_GTF}"  # Empty the file

        for SAMPLE in "${MIX_SAMPLES[@]}"; do
            GTF_FILE="${ASSEMBLY_DIR}/${SAMPLE}_MR_RB.gtf"
            if [ -f "${GTF_FILE}" ] && [ -s "${GTF_FILE}" ] && grep -qv "^#" "${GTF_FILE}"; then
                echo "${GTF_FILE}" >> "${LIST_MR_RB_GTF}"
            else
                echo_red "Warning: GTF file ${GTF_FILE} is missing, empty, or contains no transcripts. Skipping."
            fi
        done

        MERGE_INPUT_RB_GTF_LIST="${MERGE_DIR}/merge_input_RB_gtf_list.txt"
        > "${MERGE_INPUT_RB_GTF_LIST}"  # Empty the file

        if [ -f "${MERGE_DIR}/merged_SR_RB.gtf" ] && [ -s "${MERGE_DIR}/merged_SR_RB.gtf" ]; then
            echo "${MERGE_DIR}/merged_SR_RB.gtf" >> "${MERGE_INPUT_RB_GTF_LIST}"
        fi

        if [ -s "${LIST_MR_RB_GTF}" ]; then
            cat "${LIST_MR_RB_GTF}" >> "${MERGE_INPUT_RB_GTF_LIST}"
        fi

        if [ -s "${MERGE_INPUT_RB_GTF_LIST}" ]; then
            stringtie --merge -p "${THREADS}" \
                      -G "${GENOME_GTF}" \
                      -o "${MERGE_DIR}/merged_SR_MR_RB.gtf" \
                      "${MERGE_INPUT_RB_GTF_LIST}"
            echo_green "Merged SR_MR_RB.gtf created successfully."
        else
            echo_red "No valid MR_RB GTF files available for merging. Skipping."
        fi
    fi

    # Deactivate the environment
    conda deactivate

    echo_blue "Step 5 Completed: Merge Reference-Based Assemblies"
}

# ---------------------------
# Function for Step 6: Merge De Novo Assemblies and Create Pre-Final Annotation
# ---------------------------
step6_merging_transcripts_II() {
    echo_blue "Starting Step 6: Merge De Novo Assemblies and Create Pre-Final Annotation"

    conda activate stringtie221

    # Merge De Novo Assemblies for Short Reads (SR)
    if [ "${#SHORT_ONLY_SAMPLES[@]}" -gt 0 ]; then
        echo_green "Merging De Novo Assemblies for Short Reads (SR)"
        LIST_SR_DN_GTF="${MERGE_DIR}/SR_DN_gtf_list.txt"
        > "${LIST_SR_DN_GTF}"  # Empty the file

        for SAMPLE in "${SHORT_ONLY_SAMPLES[@]}"; do
            GTF_FILE="${ASSEMBLY_DIR}/${SAMPLE}_SR_DN.gtf"
            if [ -f "${GTF_FILE}" ] && [ -s "${GTF_FILE}" ] && grep -qv "^#" "${GTF_FILE}"; then
                echo "${GTF_FILE}" >> "${LIST_SR_DN_GTF}"
            else
                echo_red "Warning: GTF file ${GTF_FILE} is missing, empty, or contains no transcripts. Skipping."
            fi
        done

        if [ -s "${LIST_SR_DN_GTF}" ]; then
            stringtie --merge -p "${THREADS}" \
                      -o "${MERGE_DIR}/merged_SR_DN.gtf" \
                      "${LIST_SR_DN_GTF}"
            echo_green "Merged SR_DN.gtf created successfully."
        else
            echo_red "No valid SR_DN GTF files available for merging. Skipping."
        fi
    fi

    # Merge De Novo Assemblies for Mixed Reads (MR) with SR
    if [ "${#MIX_SAMPLES[@]}" -gt 0 ]; then
        echo_green "Merging De Novo Assemblies for Mixed Reads (MR) with SR"
        LIST_MR_DN_GTF="${MERGE_DIR}/MR_DN_gtf_list.txt"
        > "${LIST_MR_DN_GTF}"  # Empty the file

        for SAMPLE in "${MIX_SAMPLES[@]}"; do
            GTF_FILE="${ASSEMBLY_DIR}/${SAMPLE}_MR_DN.gtf"
            if [ -f "${GTF_FILE}" ] && [ -s "${GTF_FILE}" ] && grep -qv "^#" "${GTF_FILE}"; then
                echo "${GTF_FILE}" >> "${LIST_MR_DN_GTF}"
            else
                echo_red "Warning: GTF file ${GTF_FILE} is missing, empty, or contains no transcripts. Skipping."
            fi
        done

        MERGE_INPUT_DN_GTF_LIST="${MERGE_DIR}/merge_input_DN_gtf_list.txt"
        > "${MERGE_INPUT_DN_GTF_LIST}"  # Empty the file

        if [ -f "${MERGE_DIR}/merged_SR_DN.gtf" ] && [ -s "${MERGE_DIR}/merged_SR_DN.gtf" ]; then
            echo "${MERGE_DIR}/merged_SR_DN.gtf" >> "${MERGE_INPUT_DN_GTF_LIST}"
        fi

        if [ -s "${LIST_MR_DN_GTF}" ]; then
            cat "${LIST_MR_DN_GTF}" >> "${MERGE_INPUT_DN_GTF_LIST}"
        fi

        if [ -s "${MERGE_INPUT_DN_GTF_LIST}" ]; then
            stringtie --merge -p "${THREADS}" \
                      -o "${MERGE_DIR}/merged_SR_MR_DN.gtf" \
                      "${MERGE_INPUT_DN_GTF_LIST}"
            echo_green "Merged SR_MR_DN.gtf created successfully."
        else
            echo_red "No valid MR_DN GTF files available for merging. Skipping."
        fi
    fi

    # Merge all assemblies into pre-final annotation
    echo_green "Creating pre-final annotation GTF"

    # Determine which merged files are available
    MERGED_RB_GTF=""
    MERGED_DN_GTF=""
    PREFINAL_GTF="${MERGE_DIR}/prefinal_annotation.gtf"

    if [ -f "${MERGE_DIR}/merged_SR_MR_RB.gtf" ] && [ -s "${MERGE_DIR}/merged_SR_MR_RB.gtf" ]; then
        MERGED_RB_GTF="${MERGE_DIR}/merged_SR_MR_RB.gtf"
    elif [ -f "${MERGE_DIR}/merged_SR_RB.gtf" ] && [ -s "${MERGE_DIR}/merged_SR_RB.gtf" ]; then
        MERGED_RB_GTF="${MERGE_DIR}/merged_SR_RB.gtf"
    fi

    if [ -f "${MERGE_DIR}/merged_SR_MR_DN.gtf" ] && [ -s "${MERGE_DIR}/merged_SR_MR_DN.gtf" ]; then
        MERGED_DN_GTF="${MERGE_DIR}/merged_SR_MR_DN.gtf"
    elif [ -f "${MERGE_DIR}/merged_SR_DN.gtf" ] && [ -s "${MERGE_DIR}/merged_SR_DN.gtf" ]; then
        MERGED_DN_GTF="${MERGE_DIR}/merged_SR_DN.gtf"
    fi

    # Merge based on availability
    if [ -n "${MERGED_RB_GTF}" ] && [ -n "${MERGED_DN_GTF}" ]; then
        echo_green "Merging RB and DN assemblies into pre-final annotation"
        stringtie --merge -p "${THREADS}" \
                  -o "${PREFINAL_GTF}" \
                  "${MERGED_RB_GTF}" "${MERGED_DN_GTF}"
        echo_green "Pre-final annotation GTF created successfully at ${PREFINAL_GTF}"
    elif [ -n "${MERGED_DN_GTF}" ]; then
        echo_green "Using De Novo merged GTF as pre-final annotation"
        cp "${MERGE_DIR}/merged_SR_DN.gtf" "${PREFINAL_GTF}"
        echo_green "Pre-final annotation GTF created successfully at ${PREFINAL_GTF}"
    elif [ -n "${MERGED_RB_GTF}" ]; then
        echo_green "Using Reference-Based merged GTF as pre-final annotation"
        cp "${MERGE_DIR}/merged_SR_RB.gtf" "${PREFINAL_GTF}"
        echo_green "Pre-final annotation GTF created successfully at ${PREFINAL_GTF}"
    else
        echo_red "No merged GTF files available for pre-final annotation. Skipping."
    fi

    conda deactivate

    echo_blue "Step 6 Completed: Merge De Novo Assemblies and Create Pre-Final Annotation"
}

# ---------------------------
# Function for Step 6.1: Filter Out Transcripts with Excessively Long Exons or Genomic Spans
# ---------------------------
step6_filter_transcripts() {
    echo_green "Starting Step 6.1: Filtering Transcripts with Excessively Long Exons or Genomic Spans"

    # Define the maximum exon and transcript lengths
    MAX_EXON_LENGTH=10000
    MAX_TRANSCRIPT_LENGTH=100000

    # Input pre-final GTF
    PREFINAL_GTF="${MERGE_DIR}/prefinal_annotation.gtf"

    # Output filtered GTF
    FILTERED_GTF="${MERGE_DIR}/prefinal_annotation_filtered.gtf"

    # Temporary files
    BAD_TRANSCRIPTS_EXON="${MERGE_DIR}/transcripts_with_long_exons.txt"
    BAD_TRANSCRIPTS_SPAN="${MERGE_DIR}/transcripts_with_long_spans.txt"

    # Check if pre-final GTF exists
    if [ ! -f "${PREFINAL_GTF}" ]; then
        echo_red "Error: Pre-final annotation GTF (${PREFINAL_GTF}) not found. Skipping filtering."
        return
    fi

    # Step 1: Identify transcripts with any exon longer than MAX_EXON_LENGTH
    echo_green "Identifying transcripts with exons longer than ${MAX_EXON_LENGTH} nt in ${PREFINAL_GTF}"
    awk -v max_len="${MAX_EXON_LENGTH}" '
        $3 == "exon" {
            exon_length = $5 - $4 + 1
            if (exon_length > max_len) {
                # Extract transcript_id using regex
                if (match($0, /transcript_id "([^"]+)"/, arr)) {
                    print arr[1]
                }
            }
        }
    ' "${PREFINAL_GTF}" | sort | uniq > "${BAD_TRANSCRIPTS_EXON}"

    # Step 2: Identify transcripts with genomic spans longer than MAX_TRANSCRIPT_LENGTH
    echo_green "Identifying transcripts with genomic spans longer than ${MAX_TRANSCRIPT_LENGTH} nt in ${PREFINAL_GTF}"
    awk '
        $3 == "transcript" {
            transcript_id = ""
            if (match($0, /transcript_id "([^"]+)"/, arr)) {
                transcript_id = arr[1]
                span = $5 - $4 + 1
                if (span > max_span) {
                    print transcript_id
                }
            }
        }
    ' max_span="${MAX_TRANSCRIPT_LENGTH}" "${PREFINAL_GTF}" | sort | uniq > "${BAD_TRANSCRIPTS_SPAN}"

    # Combine bad transcripts
    cat "${BAD_TRANSCRIPTS_EXON}" "${BAD_TRANSCRIPTS_SPAN}" | sort | uniq > "${MERGE_DIR}/transcripts_to_exclude.txt"

    # Check if any bad transcripts were found
    if [ -s "${MERGE_DIR}/transcripts_to_exclude.txt" ]; then
        echo_red "Found $(wc -l < "${MERGE_DIR}/transcripts_to_exclude.txt") transcripts exceeding defined thresholds. These transcripts will be excluded."
    else
        echo_green "No transcripts exceeding defined thresholds were found. No filtering necessary."
        cp "${PREFINAL_GTF}" "${FILTERED_GTF}"
        echo_green "Filtered GTF is identical to the pre-final GTF."
        echo_green "Step 6.1 Completed: Filtering Transcripts with Excessively Long Exons or Genomic Spans"
        return
    fi

    # Step 3: Exclude all entries associated with the identified bad transcripts
    echo_green "Excluding identified transcripts from the GTF file"
    awk -v exclude_file="${MERGE_DIR}/transcripts_to_exclude.txt" '
        BEGIN {
            while ((getline < exclude_file) > 0) {
                exclude[$1] = 1
            }
            close(exclude_file)
        }
        /^#/ { print }
        !/^#/ {
            # Extract transcript_id
            if (match($0, /transcript_id "([^"]+)"/, arr)) {
                transcript_id = arr[1]
                if (!(transcript_id in exclude)) {
                    print
                }
            }
            else {
                # If transcript_id not found, decide to include or exclude
                # Here, we choose to include lines without transcript_id
                print
            }
        }
    ' "${PREFINAL_GTF}" > "${FILTERED_GTF}"

    # Verify that filtering was successful
    if [ -f "${FILTERED_GTF}" ]; then
        echo_green "Filtering completed. Filtered GTF saved at ${FILTERED_GTF}"
        echo_green "Excluded transcripts are listed in ${MERGE_DIR}/transcripts_to_exclude.txt"
    else
        echo_red "Error: Failed to create filtered GTF file."
        exit 1
    fi

    echo_green "Step 6.1 Completed: Filtering Transcripts with Excessively Long Exons or Genomic Spans"
}




# ---------------------------
# Function for Step 7: Isoform Comparison and Annotation
# ---------------------------
step7_isoform_comparison() {
    echo_blue  "Starting Step 7: Isoform Comparison and Annotation"

    if [ ! -f "${MERGE_DIR}/prefinal_annotation_filtered.gtf" ]; then
        echo_red  "Pre-final annotation GTF not found. Skipping Step 7."
        return
    fi

    conda activate stringtie221

    COMP_OUTPUT_PREFIX="${ANNOTATION_DIR}/gffcomp"

    echo_green  "Running gffcompare to compare assembled transcripts with reference (if provided)"
    if [ -n "${GENOME_GTF}" ]; then
        gffcompare -r "${GENOME_GTF}" \
                   -G \
                   --chr-stats \
                   --strict-match \
                   -o "${COMP_OUTPUT_PREFIX}" \
                   "${MERGE_DIR}/prefinal_annotation_filtered.gtf"
    else
        gffcompare \
                   --chr-stats \
                   --strict-match \
                   -o "${COMP_OUTPUT_PREFIX}" \
                   "${MERGE_DIR}/prefinal_annotation_filtered.gtf"
    fi

    conda deactivate

    echo_blue  "Step 7 Completed: Isoform Comparison and Annotation"
}

# ---------------------------
# Function for Step 8: GTF File Correction and Enhancement
# ---------------------------
step8_gtf_correction() {
    echo_blue  "Starting Step 8: GTF Correction and Enhancement"

    conda activate stringtie221

    if [ -f "${ANNOTATION_DIR}/gffcomp.annotated.gtf" ]; then
        ANNOTATED_GTF="${ANNOTATION_DIR}/gffcomp.annotated.gtf"
    elif [ -f "${MERGE_DIR}/prefinal_annotation_filtered.gtf" ]; then
        ANNOTATED_GTF="${MERGE_DIR}/prefinal_annotation_filtered.gtf"
    else
        echo_red  "No GTF file found for correction. Skipping Step 8."
        conda deactivate
        return
    fi

    CORRECTED_GTF="${ANNOTATION_DIR}/corrected.gtf"
    CORRECTED_INTRONS_GTF="${ANNOTATION_DIR}/corrected_with_introns.gtf"

    echo_green  "Converting and correcting GTF file using AGAT"
    # Corrected the option from -i to -g
    agat_convert_sp_gxf2gxf.pl -g "${ANNOTATED_GTF}" -o "${CORRECTED_GTF}"

    echo_green  "Adding intron features to GTF file using AGAT"
    agat_sp_add_introns.pl -g "${CORRECTED_GTF}" -o "${CORRECTED_INTRONS_GTF}"

    conda deactivate

    echo_blue  "Step 8 Completed: GTF Correction and Enhancement"
}

# ---------------------------
# Function for Step 9: Functional Annotation and Filtering
# ---------------------------
step9_functional_annotation() {
    echo_blue  "Starting Step 9: Functional Annotation and Filtering"

    conda activate stringtie221

    if [ -f "${ANNOTATION_DIR}/corrected_with_introns.gtf" ]; then
        ANNOTATED_GTF="${ANNOTATION_DIR}/corrected_with_introns.gtf"
    else
        echo_red  "No GTF file found for functional annotation. Skipping Step 9."
        conda deactivate
        return
    fi

    # Paths to transcript and protein files
    TRANSCRIPTS_FA="${ANNOTATION_DIR}/transcripts.fa"
    LONGEST_ORFS_PEP="${FUNCTIONAL_DIR}/longest_orfs.pep"
    #LONGEST_ORFS_FA="${FUNCTIONAL_DIR}/longest_orfs.fa"

    # Extract transcript sequences (requires gffread)
    echo_green  "Extracting transcript sequences using gffread"
    gffread "${ANNOTATED_GTF}" -g "${GENOME_REF}" -w "${TRANSCRIPTS_FA}"

    # ORF Prediction with TransDecoder
    echo_green  "Running TransDecoder to predict ORFs"
    TransDecoder.LongOrfs -t "${TRANSCRIPTS_FA}" -m "${MIN_ORF_LENGTH}" -O "${FUNCTIONAL_DIR}/TransDecoder"

    # Predict likely coding regions
    echo_green  "Running TransDecoder.Predict to identify likely coding regions"
    TransDecoder.Predict -t "${TRANSCRIPTS_FA}" --single_best_only -O "${FUNCTIONAL_DIR}/TransDecoder"

    # Extract longest ORFs
    echo_green  "Extracting longest ORFs"
    cp "${FUNCTIONAL_DIR}/TransDecoder/transcripts.fa.transdecoder_dir/longest_orfs.pep" "${LONGEST_ORFS_PEP}"
    #cp "${FUNCTIONAL_DIR}/TransDecoder/transcripts.fa.transdecoder_dir/longest_orfs.fa" "${LONGEST_ORFS_FA}"

    # Install SwissProt and PFAM databases if not provided
    if [ -z "${BLAST_DB_SwissProt}" ]; then
        echo_green  "Downloading and formatting SwissProt database"
        mkdir -p "${FUNCTIONAL_DIR}/blast_dbs"
        BLAST_DB_SwissProt="${FUNCTIONAL_DIR}/blast_dbs/swissprot"
        wget -O "${FUNCTIONAL_DIR}/blast_dbs/swissprot.fasta.gz" ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
        gunzip "${FUNCTIONAL_DIR}/blast_dbs/swissprot.fasta.gz"
        makeblastdb -in "${FUNCTIONAL_DIR}/blast_dbs/swissprot.fasta" -dbtype prot -out "${BLAST_DB_SwissProt}"
    fi

    if [ -z "${PFAM_DB}" ]; then
        echo_green  "Downloading PFAM database"
        mkdir -p "${FUNCTIONAL_DIR}/pfam_db"
        PFAM_DB="${FUNCTIONAL_DIR}/pfam_db"
        wget -P "${PFAM_DB}" ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
        gunzip "${PFAM_DB}/Pfam-A.hmm.gz"
        hmmpress "${PFAM_DB}/Pfam-A.hmm"
    fi

    # Homology Search with BLASTp and BLASTx
    echo_green  "Running BLASTp against SwissProt database"
    blastp -query "${LONGEST_ORFS_PEP}" \
           -db "${BLAST_DB_SwissProt}" \
           -outfmt 6 \
           -num_threads "${THREADS}" \
           -out "${FUNCTIONAL_DIR}/blastp_results.out"

    echo_green  "Running BLASTx against SwissProt database"
    blastx -query "${TRANSCRIPTS_FA}" \
           -db "${BLAST_DB_SwissProt}" \
           -outfmt 6 \
           -num_threads "${THREADS}" \
           -out "${FUNCTIONAL_DIR}/blastx_results.out"

    # Domain Identification with PFAM Scan
    echo_green  "Running PFAM Scan to identify protein domains"
    hmmscan --cpu 2 --domtblout "${FUNCTIONAL_DIR}/pfam_results.out" "${PFAM_DB}/Pfam-A.hmm" "${LONGEST_ORFS_PEP}" 
                 
    conda deactivate

    # Filtering Criteria
    echo_green  "Filtering transcripts based on ORFs, PFAM domains, and BLAST hits"

    # Create a list of transcripts with ORFs
    grep ">" "${LONGEST_ORFS_PEP}" | sed 's/>//' > "${FUNCTIONAL_DIR}/transcripts_with_orfs.txt"

    # Filter BLAST results for significant hits (e-value < 1e-5, identity > 30%)
    awk '$11 < 1e-5 && $3 > 30 {print $1}' "${FUNCTIONAL_DIR}/blastp_results.out" | sort | uniq > "${FUNCTIONAL_DIR}/blastp_filtered.txt"
    awk '$11 < 1e-5 && $3 > 30 {print $1}' "${FUNCTIONAL_DIR}/blastx_results.out" | sort | uniq > "${FUNCTIONAL_DIR}/blastx_filtered.txt"

    # Filter PFAM results for significant domains
    awk '{print $1}' "${FUNCTIONAL_DIR}/pfam_results.out" | sort | uniq > "${FUNCTIONAL_DIR}/pfam_filtered.txt"

    # Combine filters
    cat "${FUNCTIONAL_DIR}/blastp_filtered.txt" "${FUNCTIONAL_DIR}/blastx_filtered.txt" "${FUNCTIONAL_DIR}/pfam_filtered.txt" | sort | uniq > "${FUNCTIONAL_DIR}/functional_transcripts.txt"

    # Intersect with transcripts having ORFs
    #comm -12 <(sort "${FUNCTIONAL_DIR}/transcripts_with_orfs.txt") <(sort "${FUNCTIONAL_DIR}/functional_transcripts.txt") > "${FUNCTIONAL_DIR}/final_filtered_transcripts.txt"

    # Generate final annotation GTF
    echo_green  "Generating final filtered GTF file"
	#Here would be R script for filtering.....
    #grep -F -f "${FUNCTIONAL_DIR}/final_filtered_transcripts.txt" "${ANNOTATION_DIR}/corrected_with_introns.gtf" > "${FUNCTIONAL_DIR}/final_annotation.gtf"

    echo_blue  "Step 9 Completed: Functional Annotation and Filtering"
}

# ---------------------------
# Parse Command-Line Arguments with getopt
# ---------------------------
# Initialize variables with default values
GENOME_DIR="./outputDir/genome"
RRNA_REF=""
RRNA_REF_INDEX=""
GENOME_REF=""
GENOME_GTF=""
BLAST_DB_SwissProt=""
PFAM_DB=""
DATA_DIR=""
OUTPUT_DIR="./outputDir"
THREADS=8
MIN_ORF_LENGTH=100
STEPS=()
RUN_ALL=false

# Use getopt for parsing long options
if ! PARSED_OPTIONS=$(getopt -n "$0" -o "" --long genomeDir:,rrnaRef:,genomeRef:,genomeGTF:,blastDB_SwissProt:,pfamDB:,dataDir:,outputDir:,threads:,steps:,all,help,minOrfLength: -- "$@"); then
    show_help
fi

eval set -- "$PARSED_OPTIONS"

# Parse the options
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --genomeDir) GENOME_DIR="$2"; shift 2 ;;
        --rrnaRef) RRNA_REF="$2"; shift 2 ;;
        --genomeRef) GENOME_REF="$2"; shift 2 ;;
        --genomeGTF) GENOME_GTF="$2"; shift 2 ;;
        --blastDB_SwissProt) BLAST_DB_SwissProt="$2"; shift 2 ;;
        --pfamDB) PFAM_DB="$2"; shift 2 ;;
        --dataDir) DATA_DIR="$2"; shift 2 ;;
        --outputDir) OUTPUT_DIR="$2"; shift 2 ;;
        --threads) THREADS="$2"; shift 2 ;;
        --minOrfLength) MIN_ORF_LENGTH="$2"; shift 2 ;;
        --steps) IFS=',' read -ra STEPS <<< "$2"; shift 2 ;;
        --all) RUN_ALL=true; shift ;;
        --help) show_help ;;
        --) shift; break ;;
        *) echo "Unknown parameter passed: $1"; show_help ;;
    esac
done

# ---------------------------
# Validate Required Arguments
# ---------------------------
if [ -z "${GENOME_REF}" ]; then
    echo_red  "Error: --genomeRef is required."
    show_help
fi

if [ -z "${DATA_DIR}" ]; then
    echo_red  "Error: --dataDir is required."
    show_help
fi

if [ ! -d "${DATA_DIR}/short_reads" ] && [ ! -d "${DATA_DIR}/mix_reads" ]; then
    echo_red  "Error: --dataDir must contain 'short_reads' and/or 'mix_reads' directories."
    show_help
fi

if [ "$RUN_ALL" = false ] && [ ${#STEPS[@]} -eq 0 ]; then
    echo_red  "Error: You must specify --all or provide --steps to run."
    show_help
fi

# ---------------------------
# Set Up Output Directories
# ---------------------------
PREPROC_DIR="${OUTPUT_DIR}/preprocessing"
ALIGN_DIR="${OUTPUT_DIR}/alignment"
ASSEMBLY_DIR="${OUTPUT_DIR}/assembly"
MERGE_DIR="${OUTPUT_DIR}/merging"
ANNOTATION_DIR="${OUTPUT_DIR}/annotation"
FUNCTIONAL_DIR="${OUTPUT_DIR}/functional_annotation"
LOGS_DIR="${OUTPUT_DIR}/logs"

mkdir -p "${PREPROC_DIR}" "${ALIGN_DIR}" "${ASSEMBLY_DIR}" "${MERGE_DIR}" "${ANNOTATION_DIR}" "${FUNCTIONAL_DIR}" "${LOGS_DIR}"

# ---------------------------
# Set rRNA Reference Index (STAR requires genomeDir)
# ---------------------------
if [ -n "${RRNA_REF}" ]; then
    RRNA_REF_INDEX="${OUTPUT_DIR}/rRNA_STAR_index"
fi

# ---------------------------
# Define Sample Names
# ---------------------------
# Samples with only short reads
# Define Sample Names
# Samples with only short reads
SHORT_ONLY_DIR="${DATA_DIR}/short_reads"
SHORT_ONLY_SAMPLES=()
if [ -d "${SHORT_ONLY_DIR}" ]; then
    mapfile -t SHORT_ONLY_SAMPLES < <(
        find "${SHORT_ONLY_DIR}" -maxdepth 1 -type f \
            \( -name "*1.fq" -o -name "*2.fq" \
               -o -name "*1.fq.gz" -o -name "*2.fq.gz" \
               -o -name "*1.fastq" -o -name "*2.fastq" \
               -o -name "*1.fastq.gz" -o -name "*2.fastq.gz" \) \
            -printf "%f\n" | \
        sed -E 's/^(.+)[12]\.(fq|fastq)(\.gz)?$/\1/' | \
        sort | uniq
    )
fi

# Samples with mixed reads
MIX_DIR="${DATA_DIR}/mix_reads"
MIX_SHORT_DIR="${MIX_DIR}/short_reads"
MIX_LONG_DIR="${MIX_DIR}/long_reads"
MIX_SAMPLES=()
if [ -d "${MIX_SHORT_DIR}" ]; then
    mapfile -t MIX_SAMPLES < <(
        find "${MIX_SHORT_DIR}" -maxdepth 1 -type f \
            \( -name "*1.fq" -o -name "*2.fq" \
               -o -name "*1.fq.gz" -o -name "*2.fq.gz" \
               -o -name "*1.fastq" -o -name "*2.fastq" \
               -o -name "*1.fastq.gz" -o -name "*2.fastq.gz" \) \
            -printf "%f\n" | \
        sed -E 's/^(.+)[12]\.(fq|fastq)(\.gz)?$/\1/' | \
        sort | uniq
    )
fi

echo "Short-only Samples detected: ${SHORT_ONLY_SAMPLES[*]}"
echo "Mixed Samples detected: ${MIX_SAMPLES[*]}"


# ---------------------------
# Check Dependencies and Set Up Environments
# ---------------------------
check_dependencies
setup_conda_environments

# ---------------------------
# Run Selected Steps
# ---------------------------
if [ "$RUN_ALL" = true ]; then
    STEPS=(1 2 3 4 5 6 6.1 7 8 9)
fi

for STEP in "${STEPS[@]}"; do
    case "$STEP" in
        1)
            step1_genome_index
            ;;
        2)
            step2_rrna_removal
            ;;
        3)
            step3_read_alignment
            ;;
        4)
            step4_gene_transcript_assembly
            ;;
        5)
            step5_merging_transcripts_I
            ;;
        6)
            step6_merging_transcripts_II
            ;;
		6.1)
            step6_filter_transcripts
            ;;
        7)
            step7_isoform_comparison
            ;;
        8)
            step8_gtf_correction
            ;;
        9)
            step9_functional_annotation
            ;;
        *)
            echo_red  "Warning: Step $STEP is not recognized and will be skipped."
            ;;
    esac
done

# =====================================================================
# MOVE ALL LOGS
# =====================================================================
#"${LOGS_DIR}"

# =====================================================================
# Pipeline Completed
# =====================================================================

echo_green  "RNA-Seq Bioinformatics Pipeline Completed Successfully!"
if [ -f "${FUNCTIONAL_DIR}/final_annotation.gtf" ]; then
    echo_green  "Final annotated GTF file is located at ${FUNCTIONAL_DIR}/final_annotation.gtf"
elif [ -f "${ANNOTATION_DIR}/corrected_with_introns.gtf" ]; then
    echo_green  "Final annotated GTF file is located at ${ANNOTATION_DIR}/corrected_with_introns.gtf"
fi
