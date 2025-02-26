#!/bin/bash

# =====================================================================
# Comprehensive RNA-Seq annotation pipeline 
# =====================================================================
# Master script 

set -e  # Exit immediately if a command exits with a non-zero status
set -o pipefail  # Pipeline returns the exit status of the last command to fail
# set -x  # Enable debugging

# ---------------------------
# Create Linux-Native Temporary Directory for STAR
# ---------------------------
#mkdir -p /tmp/temp_star
#chmod 777 /tmp/temp_star

# ---------------------------
# Define Color Variables for Echo Messages
# ---------------------------
RED='\033[0;31m'
GREEN='\033[0;32m'
BLUE='\033[1;34m'
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
    echo "Usage: SmedAnno [OPTIONS]"
    echo ""
	echo "Steps:"
    echo "  1  - rRNA Removal"
    echo "  2  - Read Alignment to Reference Genome"
    echo "  3  - Gene and Transcript Assembly"
    echo "  4  - Merge Reference-Based Assemblies"
    echo "  5  - Merge De Novo Assemblies and Create Pre-Final Annotation"
    echo "  6  - Filter Transcripts with Excessively Long Exons or Genomic Spans"
    echo "  7  - Isoform Comparison and Annotation"
    echo "  8  - GTF File Correction and Enhancement"
    echo "  9  - Functional Annotation and Filtering"
	echo "  10 - Integrate Functional Annotation (including Overlapped Genes and Transcripts, Reversed Duplicates, Fragmmented and Chimeric Genes Identification )"
    echo ""
	echo "General Mandatory Options (if applicable based on steps selected):"
    echo "  --genomeRef PATH              Path to genome reference FASTA file"
    echo "  --dataDir PATH                Path to input RNA-Seq data directory (must contain short_reads and/or mix_reads folders)"
    echo ""
    echo "Mandatory Options for Specific Steps:"
    echo "  --alignDir PATH               Path to directory containing BAM files (required for Step 3)"
    echo "  --SR_RB_gtf_dir PATH          Directory containing SR_RB.gtf files (required for Step 4 and 5)"
    echo "  --SR_DN_gtf_dir PATH          Directory containing SR_DN.gtf files (required for Step 4 and 5)"
    echo "  --MR_RB_gtf_dir PATH          Directory containing MR_RB.gtf files (required for Step 4 and 5)"
    echo "  --MR_DN_gtf_dir PATH          Directory containing MR_DN.gtf files (required for Step 4 and 5)"
	echo "  --finalGTF PATH               Path to final GTF file (required for Step 9 and 10)"
	echo "  --outputDir PATH              Path to output directory with functional annotation folder (where PFAM, BLAST, and/or ORF outputs are stored; required for Step 10)"
    echo ""
    echo "Optional Options:"
    echo "  --genomeDir PATH              Path to STAR genome directory (will be created if not provided)"
    echo "  --genomeGTF PATH              Path to genome annotation GTF file (optional, required for Reference-Based assembly)"
    echo "  --rrnaRef PATH                Path to rRNA reference FASTA file"
    echo "  --blastDB_SwissProt PATH      Path to BLAST SwissProt database"
    echo "  --pfamDB PATH                 Path to PFAM database directory"
    echo "  --outputDir PATH              Path to output directory (default: ./outputDir)"
    echo "  --threads N                   Number of CPU threads to use (default: 8)"
    echo "  --minOrfLength N              Minimum ORF length for TransDecoder (default: 100)"
	echo "  --maxExonLength N             Maximum allowed exon length (default: 10000)"           
    echo "  --maxTranscriptLength N       Maximum allowed transcript length (default: 100000)"      
    echo "  --steps LIST                  Comma-separated list of steps to run (1-10)"
    echo "  --all                         Run all steps sequentially"
    echo "  --functionalMethods METHODS   Comma-separated list of functional annotation methods to apply (BLASTp,BLASTx,PFAM; default: all)"
	echo "  --conda VERSION               Set stringtie version for both short and mix reads (if not using individual overrides)"
    echo "  --conda_short VERSION         Set stringtie version for short reads (default: 2.1.1)"
    echo "  --conda_mix VERSION           Set stringtie version for mixed reads (default: 2.2.1)"
    echo "  --help                        Display this help message and exit"
    echo ""
    echo "Examples:"
    echo "  Run all steps:"
    echo "    ./smedanno.sh --genomeRef genome.fa --dataDir ./data --outputDir ./output --threads 4 --all"
    echo ""
	echo "  Run all steps with custom stringtie version 3.0.0 for both environments:"
    echo "    ./smedanno.sh --genomeRef genome.fa --dataDir ./data --outputDir ./output --threads 4 --all --conda 3.0.0"
    echo ""
    echo "  Run Steps 1 and 2 only (rRNA Removal and Read Alignment):"
    echo "    ./smedanno.sh --genomeRef genome.fa --dataDir ./data --outputDir ./output --threads 4 --steps 1,2"
    echo ""
    echo "  Run Only Step 3 (Gene and Transcript Assembly) with BAM files:"
    echo "    ./smedanno.sh --alignDir ./bam_files --outputDir ./output --threads 4 --steps 3"
    echo ""
    echo "  Run Merging Assemblies Steps 4 and 5 with specific GTF directories:"
    echo "    ./smedanno.sh --SR_RB_gtf_dir ./gtf/SR_RB --SR_DN_gtf_dir ./gtf/SR_DN --outputDir ./output --threads 4 --steps 4,5"
    echo ""
    echo "  Run Functional Annotation with specific methods:"
    echo "    ./smedanno.sh --finalGTF corrected_with_introns.gtf --outputDir ./output --threads 4 --steps 8 --functionalMethods BLASTp,PFAM --genomeRef genome.fa"
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
# Define a helper function to check if a step is requested
# ---------------------------

contains_step() {
    local target="$1"
    shift
    for elem in "$@"; do
        if [ "$elem" = "$target" ]; then
            return 0
        fi
    done
    return 1
}

# ---------------------------
# Global Default Versions for Conda Environments
# ---------------------------
DEFAULT_CONDA_SHORT="2.1.1"
DEFAULT_CONDA_MIX="2.2.1"
CONDA_SHORT_VERSION="${DEFAULT_CONDA_SHORT}"
CONDA_MIX_VERSION="${DEFAULT_CONDA_MIX}"
CONDA_VERSION=""
ENV_APPS="bioapps_env"
ENV_R="R_env"


# ---------------------------
# Default maximum lengths (can be overridden by user)
# ---------------------------
DEFAULT_MAX_EXON_LENGTH=10000
DEFAULT_MAX_TRANSCRIPT_LENGTH=100000
MAX_EXON_LENGTH="${DEFAULT_MAX_EXON_LENGTH}"
MAX_TRANSCRIPT_LENGTH="${DEFAULT_MAX_TRANSCRIPT_LENGTH}"

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
    echo_green "Setting up conda environments"
    source "$(conda info --base)/etc/profile.d/conda.sh" || { echo_red "Error: Failed to source conda.sh"; exit 1; }
	conda config --set channel_priority flexible
	
	# Create the logs directory if it doesn't exist
    mkdir -p "${OUTPUT_DIR}/logs"

    if ! conda env list | grep -q "${ENV_SHORT}"; then
        echo_green "Creating conda environment: ${ENV_SHORT} with stringtie version ${CONDA_SHORT_VERSION}"
        if ! conda create -y -n "${ENV_SHORT}" stringtie=${CONDA_SHORT_VERSION} &> "${OUTPUT_DIR}/logs/conda_create_${ENV_SHORT}.log"; then
            echo_red "Error: Failed to create ${ENV_SHORT} environment. Check ${OUTPUT_DIR}/logs/conda_create_${ENV_SHORT}.log for details."
            exit 1
        fi
    else
        echo_green "Conda environment '${ENV_SHORT}' already exists."
    fi

    if ! conda env list | grep -q "${ENV_MIX}"; then
        echo_green "Creating conda environment: ${ENV_MIX} with stringtie version ${CONDA_MIX_VERSION} and additional tools"
        if ! conda create -y -n "${ENV_MIX}" stringtie=${CONDA_MIX_VERSION} &> "${OUTPUT_DIR}/logs/conda_create_${ENV_MIX}.log"; then
            echo_red "Error: Failed to create ${ENV_MIX} environment. Check ${OUTPUT_DIR}/logs/conda_create_${ENV_MIX}.log for details."
            exit 1
        fi
    else
        echo_green "Conda environment '${ENV_MIX}' already exists."
    fi
	
	if ! conda env list | grep -q "${ENV_APPS}"; then
		echo_green "Creating conda environment: ${ENV_APPS} with gffcompare, agat, TransDecoder, blast, hmmer, pfam_scan, gffread, samtools, star, and minimap2"
		if ! conda create -y -n "${ENV_APPS}" gffcompare agat TransDecoder blast hmmer pfam_scan gffread samtools star minimap2&> "${OUTPUT_DIR}/logs/conda_create_${ENV_APPS}.log"; then
			echo_red "Error: Failed to create ${ENV_APPS} environment. Check ${OUTPUT_DIR}/logs/conda_create_${ENV_APPS}.log for details."
			exit 1
		fi
	else
    echo_green "Conda environment '${ENV_APPS}' already exists."
	fi 
	
	if ! conda env list | grep -q "${ENV_R}"; then
		echo_green "Creating conda environment: ${ENV_R} with r-base, rtracklayer, genomicfeatures, genomicranges, zlib, libcurl, and libxml2"
		if ! conda create -y -n "${ENV_R}" r-base=4.4.2 zlib libcurl libxml2 bioconductor-rtracklayer bioconductor-genomicfeatures bioconductor-genomicranges&> "${OUTPUT_DIR}/logs/conda_create_${ENV_R}.log"; then
			echo_red "Error: Failed to create ${ENV_R} environment. Check ${OUTPUT_DIR}/logs/conda_create_${ENV_R}.log for details."
			exit 1
		fi
	else
    echo_green "Conda environment '${ENV_R}' already exists."
    fi


    echo_green "Conda environments set up successfully."
	}

# ---------------------------
# Helper Function: contains
# Checks if a value exists in a list
# ---------------------------
contains() {
    local e
    for e in "$@"; do
        if [[ "$e" == "$1" ]]; then
            return 0
        fi
    done
    return 1
}

# ---------------------------
# Function for Step 1: Preprocessing - rRNA Removal
# ---------------------------
step1_rrna_removal() {
    echo_blue  "Starting Step 1: rRNA Removal"
    conda activate ${ENV_APPS}

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
		TMP_DIR="${OUTPUT_DIR}/tmp_star"
		rm -rf "${TMP_DIR}"
		mkdir -p "${TMP_DIR}"
		chmod 777 "${TMP_DIR}"
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
			TMP_DIR="${OUTPUT_DIR}/tmp_star"
		    rm -rf "${TMP_DIR}"
		    mkdir -p "${TMP_DIR}"
		    chmod 777 "${TMP_DIR}"
			
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
            echo_green "Processing Mixed Sample: ${SAMPLE}"

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

            # Align short reads using STAR (sorted BAM)
            echo_green  "Running STAR aligner for Mixed Sample: ${SAMPLE} (Short Reads)"
			TMP_DIR="${OUTPUT_DIR}/tmp_star"
			rm -rf "${TMP_DIR}"
			mkdir -p "${TMP_DIR}"
			chmod 777 "${TMP_DIR}"
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
			# Output files
            NON_RRNA_READ1="${PREPROC_DIR}/${SAMPLE}_non_rrna_R1.fq.gz"
            NON_RRNA_READ2="${PREPROC_DIR}/${SAMPLE}_non_rrna_R2.fq.gz"
			# Extract non-rRNA reads
            gzip "${PREPROC_DIR}/${SAMPLE}_STAR_rrna_Unmapped.out.mate1"
            gzip "${PREPROC_DIR}/${SAMPLE}_STAR_rrna_Unmapped.out.mate2"
            mv "${PREPROC_DIR}/${SAMPLE}_STAR_rrna_Unmapped.out.mate1.gz" "${NON_RRNA_READ1}"
            mv "${PREPROC_DIR}/${SAMPLE}_STAR_rrna_Unmapped.out.mate2.gz" "${NON_RRNA_READ2}"

            # Define paths
            ALIGNED_SORTED_BAM="${PREPROC_DIR}/${SAMPLE}_STAR_rrna_Unmapped.out.bam"

            # Check if STAR produced the BAM file
            if [ ! -f "${ALIGNED_SORTED_BAM}" ]; then
                echo_red  "Error: STAR did not produce the expected BAM file for mixed sample ${SAMPLE}."
                exit 1
            fi

            # **Remove XS Attribute Addition:**
            # Commenting out the XS attribute addition as it's not needed here
            # echo_green  "Adding XS attribute to short BAM file for Sample: ${SAMPLE}"
            # samtools view -h "${ALIGNED_SORTED_BAM}" | \
            #     awk -v strType=2 -f ./tagXSstrandedData.awk | \
            #     samtools view -bSq 10 -F 4 - > "${ALIGN_DIR}/${SAMPLE}_STAR_Aligned.sortedByCoord.outXS.bam"

            # **Extract Unmapped Reads from Long Reads:**

            # Process long reads
            # Handle both .fq and .fastq extensions
            MIX_LONG_READ=$(find_read_file "${MIX_LONG_DIR}" "${SAMPLE}" "_long_reads")

            if [ -z "${MIX_LONG_READ}" ]; then
                echo_red  "No long reads found for Mixed Sample: ${SAMPLE}"
                continue
            fi

            MIX_NON_RRNA_LONG_READ="${PREPROC_DIR}/${SAMPLE}_non_rrna_long_reads.fq.gz"

            echo_green  "Removing rRNA from mixed sample's long reads using minimap2 for Sample: ${SAMPLE}"
            
            # Remove rRNA from long reads using minimap2
            if [[ "${MIX_LONG_READ}" == *.gz ]]; then
                minimap2 -t "${THREADS}" -ax map-pb "${RRNA_REF}" <(zcat "${MIX_LONG_READ}") | \
                    samtools view -Sb - > "${ALIGN_DIR}/${SAMPLE}_long_aligned.bam"
            else
                minimap2 -t "${THREADS}" -ax map-pb "${RRNA_REF}" "${MIX_LONG_READ}" | \
                    samtools view -Sb - > "${ALIGN_DIR}/${SAMPLE}_long_aligned.bam"
            fi

            # **Sort Long Reads BAM:**
            if [ -f "${ALIGN_DIR}/${SAMPLE}_long_aligned.bam" ]; then
                echo_green  "Sorting long reads BAM file for Sample: ${SAMPLE}"
                samtools sort -@ "${THREADS}" -o "${ALIGN_DIR}/${SAMPLE}_long_aligned.sorted.bam" "${ALIGN_DIR}/${SAMPLE}_long_aligned.bam"
                rm "${ALIGN_DIR}/${SAMPLE}_long_aligned.bam"  # Remove unsorted BAM
                mv "${ALIGN_DIR}/${SAMPLE}_long_aligned.sorted.bam" "${ALIGN_DIR}/${SAMPLE}_long_aligned.bam"  # Rename to original name
                echo_green  "Long reads BAM sorted for Sample: ${SAMPLE}"
            fi

            # **Extract Unmapped (Non-rRNA) Long Reads and Convert to FASTQ:**
            echo_green  "Extracting unmapped (non-rRNA) long reads for Sample: ${SAMPLE}"
            samtools view -b -f 4 "${ALIGN_DIR}/${SAMPLE}_long_aligned.bam" > "${ALIGN_DIR}/${SAMPLE}_long_unmapped.bam"

            echo_green  "Converting unmapped long reads BAM to paired FASTQ for Sample: ${SAMPLE}"
            samtools fastq "${MIX_NON_RRNA_LONG_READ}" \
                          "${ALIGN_DIR}/${SAMPLE}_long_unmapped.bam"

            # **Clean Up Intermediate Files:**
            rm "${ALIGN_DIR}/${SAMPLE}_long_aligned.bam" "${ALIGN_DIR}/${SAMPLE}_long_unmapped.bam"

            echo_green  "Non-rRNA long reads extracted and converted to FASTQ for Sample: ${SAMPLE}"
        done
    fi
}


    # ---------------------------
    # Function for Step 2: Read Alignment to Reference Genome
    # ---------------------------
    step2_read_alignment() {
        echo_blue  "Starting Step 2: Read Alignment to the Reference Genome"
        conda activate ${ENV_APPS}
		
		if [ -z "${GENOME_DIR}" ]; then
			GENOME_DIR="${OUTPUT_DIR}/genome"
			echo_green "Step 2 requires a genome index. GENOME_DIR not provided; auto-setting GENOME_DIR to ${GENOME_DIR}"
		fi

        # Integrate Genome Index Generation within Step 2
        if [ ! -d "${GENOME_DIR}" ] || [ ! -f "${GENOME_DIR}/SA" ]; then
            echo_green  "Genome index not found. Generating genome index as part of Step 2."
            # Calculate genome length
            GENOME_LENGTH=$(awk '/^>/ {next} {len+=length($0)} END {print len}' "${GENOME_REF}")

            # Calculate genomeSAindexNbases
            SA_INDEX_NBASES=$(python3 -c "import math; l=math.log(${GENOME_LENGTH},2); n=int(l/2 -1); print(min(14, n))")

            echo_green  "Genome length: ${GENOME_LENGTH}"
            echo_green  "Calculated genomeSAindexNbases: ${SA_INDEX_NBASES}"

            echo_green  "Generating genome index with STAR at ${GENOME_DIR}"
			TMP_DIR="${OUTPUT_DIR}/tmp_star"
		    rm -rf "${TMP_DIR}"
		    mkdir -p "${TMP_DIR}"
		    chmod 777 "${TMP_DIR}"
            mkdir -p "${GENOME_DIR}"
            STAR --runThreadN "${THREADS}" \
                 --runMode genomeGenerate \
                 --genomeDir "${GENOME_DIR}" \
                 --genomeFastaFiles "${GENOME_REF}" \
                 --genomeSAindexNbases "${SA_INDEX_NBASES}" \
                 --outTmpDir /tmp/temp_star/temp

            echo_green  "Genome index generation completed."
        else
            echo_green  "Genome index already exists at ${GENOME_DIR}, skipping index generation."
        fi

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
				TMP_DIR="${OUTPUT_DIR}/tmp_star"
		        rm -rf "${TMP_DIR}"
		        mkdir -p "${TMP_DIR}"
		        chmod 777 "${TMP_DIR}"
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

                # Remove intermediate BAM file
                rm "${ALIGNED_SORTED_BAM}"

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
				TMP_DIR="${OUTPUT_DIR}/tmp_star"
		        rm -rf "${TMP_DIR}"
		        mkdir -p "${TMP_DIR}"
		        chmod 777 "${TMP_DIR}"
                STAR --runThreadN "${THREADS}" \
                     --genomeDir "${GENOME_DIR}" \
                     --readFilesIn "${MIX_READ1}" "${MIX_READ2}" \
                     --readFilesCommand "${MIX_READ_FILES_CMD}" \
                     --outFilterType BySJout \
                     --outFilterMultimapNmax 1 \
                     --outFilterMismatchNmax 0 \
                     --outReadsUnmapped Fastx \
                     --outSAMtype BAM SortedByCoordinate \
                     --outFileNamePrefix "${ALIGN_DIR}/${SAMPLE}_STAR_" \
                     --outTmpDir /tmp/temp_star/temp

                # Define paths
                ALIGNED_SHORT_SORTED_BAM="${ALIGN_DIR}/${SAMPLE}_STAR_Aligned.sortedByCoord.out.bam"

                # Check if STAR produced the BAM file
                if [ ! -f "${ALIGNED_SHORT_SORTED_BAM}" ]; then
                    echo_red  "Error: STAR did not produce the expected sorted BAM file for mixed sample ${SAMPLE}."
                    exit 1
                fi

                # Add XS attribute
                echo_green  "Adding XS attribute to short BAM file for Sample: ${SAMPLE}"
                samtools view -h "${ALIGNED_SHORT_SORTED_BAM}" | \
                    awk -v strType=2 -f ./tagXSstrandedData.awk | \
                    samtools view -bSq 10 -F 4 - > "${ALIGN_DIR}/${SAMPLE}_STAR_Aligned.sortedByCoord.outXS.bam"

                # Remove intermediate BAM file
                rm "${ALIGNED_SHORT_SORTED_BAM}"

                # Align long reads using minimap2
                if [ -n "${RRNA_REF}" ]; then
                    MIX_LONG_READ="${MIX_NON_RRNA_LONG_READ}"
                else
                    MIX_LONG_READ=$(find_read_file "${MIX_LONG_DIR}" "${SAMPLE}" "_long_reads")
                fi

                if [ -n "${MIX_LONG_READ}" ] && [ -f "${MIX_LONG_READ}" ]; then
                    echo_green  "Running minimap2 aligner for Mixed Sample: ${SAMPLE} (Long Reads)"
                    if [[ "${MIX_LONG_READ}" == *.gz ]]; then
                        minimap2 -t "${THREADS}" -ax map-pb "${GENOME_REF}" <(zcat "${MIX_LONG_READ}") | \
                            samtools view -Sb - > "${ALIGN_DIR}/${SAMPLE}_long_aligned.bam"
                    else
                        minimap2 -t "${THREADS}" -ax map-pb "${GENOME_REF}" "${MIX_LONG_READ}" | \
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
                else
                    echo_green  "No long reads found or file does not exist for Mixed Sample: ${SAMPLE}"
                fi

                echo_green  "Alignment completed for Mixed Sample: ${SAMPLE}"
            done
        fi

        echo_blue  "Step 2 Completed: Read Alignment"
        conda deactivate 
    }

    # ---------------------------
    # Function for Step 3: Gene and Transcript Assembly
    # ---------------------------
    step3_gene_transcript_assembly() {
        echo_blue  "Starting Step 3: Gene and Transcript Assembly"

        # Ensure ALIGN_DIR is provided and exists
        if [ -z "${ALIGN_DIR}" ] || [ ! -d "${ALIGN_DIR}" ]; then
            echo_red  "Error: --alignDir must be provided and must exist when running Step 3 (Gene and Transcript Assembly)."
            exit 1
        fi

        # Process short-only samples
        if [ "${#SHORT_ONLY_SAMPLES[@]}" -gt 0 ]; then
            for SAMPLE in "${SHORT_ONLY_SAMPLES[@]}"; do
                echo_green  "Assembling transcripts for Short-Only Sample: ${SAMPLE}"

                ALIGNED_SORTED_BAM="${ALIGN_DIR}/${SAMPLE}_STAR_Aligned.sortedByCoord.outXS.bam"

                # Reference-Based Annotation (RB) with StringTie2 (short reads)
                if [ -n "${GENOME_GTF}" ]; then
                    echo_green  "Running StringTie2  for Reference-Based Assembly (RB, SR)"
                    conda activate ${ENV_SHORT}
                    stringtie -p "${THREADS}" \
                              -G "${GENOME_GTF}" \
                              -c 1.5 \
                              -f 0.02 \
                              -o "${SR_RB_GTF_DIR}/${SAMPLE}_SR_RB.gtf" \
                              "${ALIGNED_SORTED_BAM}"
                    conda deactivate
                fi

                # De Novo Annotation (DN) with StringTie2 (short reads)
                echo_green  "Running StringTie2 for De Novo Assembly (DN, SR)"
                conda activate ${ENV_SHORT}
                stringtie -p "${THREADS}" \
                          -c 1.5 \
                          -f 0.02 \
                          -o "${SR_DN_GTF_DIR}/${SAMPLE}_SR_DN.gtf" \
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

                # Reference-Based Annotation (RB) with StringTie2 (Mixed Reads)
                if [ -n "${GENOME_GTF}" ]; then
                    echo_green  "Running StringTie2 for Reference-Based Assembly (RB, MR)"
                    conda activate ${ENV_MIX}
                    stringtie --mix -p "${THREADS}" \
                              -G "${GENOME_GTF}" \
                              -c 1.5 \
                              -f 0.02 \
                              -o "${MR_RB_GTF_DIR}/${SAMPLE}_MR_RB.gtf" \
                              "${ALIGNED_SHORT_SORTED_BAM}" "${ALIGNED_LONG_BAM}"
                    conda deactivate
                fi

                # De Novo Annotation (DN) with StringTie2 (Mixed Reads)
                echo_green  "Running StringTie2 for De Novo Assembly (DN, MR)"
                conda activate ${ENV_MIX}
                stringtie --mix -p "${THREADS}" \
                          -c 1.5 \
                          -f 0.02 \
                          -o "${MR_DN_GTF_DIR}/${SAMPLE}_MR_DN.gtf" \
                          "${ALIGNED_SHORT_SORTED_BAM}" "${ALIGNED_LONG_BAM}"
                conda deactivate

                echo_green  "Transcript assembly completed for Mixed Sample: ${SAMPLE}"
            done
        fi

        echo_blue  "Step 3 Completed: Gene and Transcript Assembly"
    }

    # ---------------------------
    # Function for Step 4: Merge Reference-Based Assemblies
    # ---------------------------
    step4_merging_transcripts_I() {
        echo_blue "Starting Step 4: Merge Reference-Based Assemblies"

        if [ -z "${GENOME_GTF}" ]; then
            echo_red "No reference GTF provided. Skipping Reference-Based merging."
            return
        fi

        # Activate StringTie environment
        conda activate ${ENV_MIX}

        # Merge Reference-Based Assemblies for Short Reads (SR)
        if [ -d "${SR_RB_GTF_DIR}" ]; then
            echo_green "Merging Reference-Based Assemblies for Short Reads (SR)"
            LIST_SR_RB_GTF="${MERGE_DIR}/SR_RB_gtf_list.txt"
            > "${LIST_SR_RB_GTF}"  # Empty the file

            for GTF_FILE in "${SR_RB_GTF_DIR}"/*.gtf; do
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
        else
            echo_red "Directory for SR_RB.gtf files (${SR_RB_GTF_DIR}) not found or empty. Skipping."
        fi

        # Merge Reference-Based Assemblies for Mixed Reads (MR)
        if [ -d "${MR_RB_GTF_DIR}" ]; then
            echo_green "Merging Reference-Based Assemblies for Mixed Reads (MR)"
            LIST_MR_RB_GTF="${MERGE_DIR}/MR_RB_gtf_list.txt"
            > "${LIST_MR_RB_GTF}"  # Empty the file

            for GTF_FILE in "${MR_RB_GTF_DIR}"/*.gtf; do
                if [ -f "${GTF_FILE}" ] && [ -s "${GTF_FILE}" ] && grep -qv "^#" "${GTF_FILE}"; then
                    echo "${GTF_FILE}" >> "${LIST_MR_RB_GTF}"
                else
                    echo_red "Warning: GTF file ${GTF_FILE} is missing, empty, or contains no transcripts. Skipping."
                fi
            done

            if [ -s "${LIST_MR_RB_GTF}" ]; then
                stringtie --merge -p "${THREADS}" \
                          -G "${GENOME_GTF}" \
                          -o "${MERGE_DIR}/merged_MR_RB.gtf" \
                          "${LIST_MR_RB_GTF}"
                echo_green "Merged MR_RB.gtf created successfully."
            else
                echo_red "No valid MR_RB GTF files available for merging. Skipping."
            fi
        else
            echo_red "Directory for MR_RB.gtf files (${MR_RB_GTF_DIR}) not found or empty. Skipping."
        fi

        # Deactivate the environment
        conda deactivate

        echo_blue "Step 4 Completed: Merge Reference-Based Assemblies"
    }

    # ---------------------------
    # Function for Step 5: Merge De Novo Assemblies and Create Pre-Final Annotation
    # ---------------------------
    step5_merging_transcripts_II() {
        echo_blue "Starting Step 5: Merge De Novo Assemblies and Create Pre-Final Annotation"

        conda activate ${ENV_MIX}

        # Merge De Novo Assemblies for Short Reads (SR)
        if [ -d "${SR_DN_GTF_DIR}" ]; then
            echo_green "Merging De Novo Assemblies for Short Reads (SR)"
            LIST_SR_DN_GTF="${MERGE_DIR}/SR_DN_gtf_list.txt"
            > "${LIST_SR_DN_GTF}"  # Empty the file

            for GTF_FILE in "${SR_DN_GTF_DIR}"/*.gtf; do
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
        else
            echo_red "Directory for SR_DN.gtf files (${SR_DN_GTF_DIR}) not found or empty. Skipping."
        fi

        # Merge De Novo Assemblies for Mixed Reads (MR)
        if [ -d "${MR_DN_GTF_DIR}" ]; then
            echo_green "Merging De Novo Assemblies for Mixed Reads (MR)"
            LIST_MR_DN_GTF="${MERGE_DIR}/MR_DN_gtf_list.txt"
            > "${LIST_MR_DN_GTF}"  # Empty the file

            for GTF_FILE in "${MR_DN_GTF_DIR}"/*.gtf; do
                if [ -f "${GTF_FILE}" ] && [ -s "${GTF_FILE}" ] && grep -qv "^#" "${GTF_FILE}"; then
                    echo "${GTF_FILE}" >> "${LIST_MR_DN_GTF}"
                else
                    echo_red "Warning: GTF file ${GTF_FILE} is missing, empty, or contains no transcripts. Skipping."
                fi
            done

            if [ -s "${LIST_MR_DN_GTF}" ]; then
                stringtie --merge -p "${THREADS}" \
                          -o "${MERGE_DIR}/merged_MR_DN.gtf" \
                          "${LIST_MR_DN_GTF}"
                echo_green "Merged MR_DN.gtf created successfully."
            else
                echo_red "No valid MR_DN GTF files available for merging. Skipping."
            fi
        else
            echo_red "Directory for MR_DN.gtf files (${MR_DN_GTF_DIR}) not found or empty. Skipping."
        fi

        # Merge all assemblies into pre-final annotation
        echo_green "Creating pre-final annotation GTF"

        # Determine which merged files are available
        MERGED_RB_GTF=""
        MERGED_DN_GTF=""
        PREFINAL_GTF="${MERGE_DIR}/prefinal_annotation.gtf"

        if [ -f "${MERGE_DIR}/merged_SR_RB.gtf" ] && [ -s "${MERGE_DIR}/merged_SR_RB.gtf" ]; then
            MERGED_RB_GTF="${MERGE_DIR}/merged_SR_RB.gtf"
        elif [ -f "${MERGE_DIR}/merged_MR_RB.gtf" ] && [ -s "${MERGE_DIR}/merged_MR_RB.gtf" ]; then
            MERGED_RB_GTF="${MERGE_DIR}/merged_MR_RB.gtf"
        fi

        if [ -f "${MERGE_DIR}/merged_SR_DN.gtf" ] && [ -s "${MERGE_DIR}/merged_SR_DN.gtf" ]; then
            MERGED_DN_GTF="${MERGE_DIR}/merged_SR_DN.gtf"
        elif [ -f "${MERGE_DIR}/merged_MR_DN.gtf" ] && [ -s "${MERGE_DIR}/merged_MR_DN.gtf" ]; then
            MERGED_DN_GTF="${MERGE_DIR}/merged_MR_DN.gtf"
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

        echo_blue "Step 5 Completed: Merge De Novo Assemblies and Create Pre-Final Annotation"
    }

    # ---------------------------
    # Function for Step 6: Filter Out Transcripts with Excessively Long Exons or Genomic Spans
    # ---------------------------
    step6_filter_transcripts() {
        echo_green "Starting Step 6: Filtering Transcripts with Excessively Long Exons or Genomic Spans"
        conda activate ${ENV_APPS}
        # Define the maximum exon and transcript lengths
        #MAX_EXON_LENGTH=10000
        #MAX_TRANSCRIPT_LENGTH=100000

        # Input pre-final GTF
        PREFINAL_GTF="${MERGE_DIR}/prefinal_annotation.gtf"

        # Output filtered GTF
        FILTERED_GTF="${MERGE_DIR}/filtered_annotation.gtf"

        # Temporary files
        BAD_TRANSCRIPTS_EXON="${MERGE_DIR}/transcripts_with_long_exons.txt"
        BAD_TRANSCRIPTS_SPAN="${MERGE_DIR}/transcripts_with_long_spans.txt"

        # Check if pre-final GTF exists
        if [ ! -f "${PREFINAL_GTF}" ]; then
            echo_red "Error: Pre-final annotation GTF (${PREFINAL_GTF}) not found. Skipping filtering."
            return
        fi

        # Identify transcripts with any exon longer than MAX_EXON_LENGTH
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

        # Identify transcripts with genomic spans longer than MAX_TRANSCRIPT_LENGTH
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
            echo_green "Step 6 Completed: Filtering Transcripts with Excessively Long Exons or Genomic Spans"
            return
        fi

        # Exclude all entries associated with the identified bad transcripts
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
		conda deactivate

        echo_green "Step 6 Completed: Filtering Transcripts with Excessively Long Exons or Genomic Spans"
    }

    # ---------------------------
    # Function for Step 7: Isoform Comparison and Annotation
    # ---------------------------
    step7_isoform_comparison() {
        echo_blue  "Starting Step 7: Isoform Comparison and Annotation"

        # Determine which GTF to use
        if [ -f "${MERGE_DIR}/filtered_annotation.gtf" ]; then
            ANNOTATED_GTF="${MERGE_DIR}/filtered_annotation.gtf"
        elif [ -n "${finalGTF}" ] && [ -f "${finalGTF}" ]; then
            ANNOTATED_GTF="${finalGTF}"
        else
            echo_red  "Pre-final annotation GTF not found. Skipping Step 7."
            return
        fi

        conda activate ${ENV_APPS}

        COMP_OUTPUT_PREFIX="${ANNOTATION_DIR}/gffcomp"

        echo_green  "Running gffcompare to compare assembled transcripts with reference (if provided)"
        if [ -n "${GENOME_GTF}" ]; then
            gffcompare -r "${GENOME_GTF}" \
                       -G \
                       --chr-stats \
                       --strict-match \
                       -o "${COMP_OUTPUT_PREFIX}" \
                       "${ANNOTATED_GTF}"
        else
            gffcompare \
                       --chr-stats \
                       --strict-match \
                       -o "${COMP_OUTPUT_PREFIX}" \
                       "${ANNOTATED_GTF}"
        fi

        conda deactivate

        echo_blue  "Step 7 Completed: Isoform Comparison and Annotation"
    }

    # ---------------------------
    # Function for Step 8: GTF File Correction and Enhancement
    # ---------------------------
    step8_gtf_correction() {
        echo_blue  "Starting Step 8: GTF Correction and Enhancement"

        conda activate ${ENV_APPS}

        # Determine which GTF to use for correction
        if [ -f "${ANNOTATION_DIR}/gffcomp.annotated.gtf" ]; then
            ANNOTATED_GTF="${ANNOTATION_DIR}/gffcomp.annotated.gtf"
        elif [ -f "${MERGE_DIR}/filtered_annotation.gtf" ]; then
            ANNOTATED_GTF="${MERGE_DIR}/filtered_annotation.gtf"
        elif [ -n "${finalGTF}" ] && [ -f "${finalGTF}" ]; then
            ANNOTATED_GTF="${finalGTF}"
        else
            echo_red  "No GTF file found for correction. Skipping Step 8."
            conda deactivate
            return
        fi

        CORRECTED_GTF="${ANNOTATION_DIR}/corrected.gtf"
        CORRECTED_INTRONS_GTF="${ANNOTATION_DIR}/corrected_with_introns.gtf"

        echo_green  "Converting and correcting GTF file using AGAT"
        agat_convert_sp_gxf2gxf.pl -g "${ANNOTATED_GTF}" -o "${CORRECTED_GTF}"

        echo_green  "Adding intron features to GTF file using AGAT"
        agat_sp_add_introns.pl -g "${CORRECTED_GTF}" -o "${CORRECTED_INTRONS_GTF}"
		
		# Move all AGAT log files (*.agat.log) from the current directory to the logs directory in the output folder
		if compgen -G "./*.agat.log" > /dev/null; then
			echo_green "Moving AGAT log files to ${OUTPUT_DIR}/logs/"
			mv ./*.agat.log "${OUTPUT_DIR}/logs/"
		else
			echo_green "No AGAT log files (*.agat.log) found in the current directory."
		fi

        conda deactivate

        echo_blue  "Step 8 Completed: GTF Correction and Enhancement"
    }

    # ---------------------------
    # Function for Step 9: Functional Annotation and Filtering
    # ---------------------------
    step9_functional_annotation() {
        echo_blue  "Starting Step 9: Functional Annotation and Filtering"

        if [ -f "${ANNOTATION_DIR}/corrected_with_introns.gtf" ]; then
            ANNOTATED_GTF="${ANNOTATION_DIR}/corrected_with_introns.gtf"
        elif [ -n "${finalGTF}" ] && [ -f "${finalGTF}" ]; then
            ANNOTATED_GTF="${finalGTF}"
        else
            echo_red  "No GTF file found for functional annotation. Skipping Step 9."
            return
        fi

        conda activate ${ENV_APPS}

        # Paths to transcript and protein files
        TRANSCRIPTS_FA="${ANNOTATION_DIR}/transcripts.fa"
        LONGEST_ORFS_PEP="${FUNCTIONAL_DIR}/longest_orfs.pep"

        # Extract transcript sequences (requires gffread)
        echo_green  "Extracting transcript sequences using gffread"
        gffread "${ANNOTATED_GTF}" -g "${GENOME_REF}" -w "${TRANSCRIPTS_FA}"

        # ORF Prediction with TransDecoder
        echo_green  "Running TransDecoder to predict ORFs"
        TransDecoder.LongOrfs -t "${TRANSCRIPTS_FA}" -m "${MIN_ORF_LENGTH}" -O "${FUNCTIONAL_DIR}/TransDecoder"
        #TransDecoder.LongOrfs -t "${TRANSCRIPTS_FA}" --min_prot_len "${MIN_ORF_LENGTH}"

        # Predict likely coding regions
        echo_green  "Running TransDecoder.Predict to identify likely coding regions"
        TransDecoder.Predict -t "${TRANSCRIPTS_FA}" --single_best_only -O "${FUNCTIONAL_DIR}/TransDecoder"

        # Extract longest ORFs
        echo_green  "Extracting longest ORFs"
        cp "${FUNCTIONAL_DIR}/TransDecoder/transcripts.fa.transdecoder_dir/longest_orfs.pep" "${LONGEST_ORFS_PEP}"

        # Install SwissProt and PFAM databases if not provided
        if [ -z "${BLAST_DB_SwissProt}" ]; then
            echo_green  "Downloading and formatting SwissProt database"
            mkdir -p "${FUNCTIONAL_DIR}/blast_dbs"
            BLAST_DB_SwissProt="${FUNCTIONAL_DIR}/blast_dbs/swissprot"
            wget -O "${FUNCTIONAL_DIR}/blast_dbs/swissprot.fasta.gz" ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
            gunzip "${FUNCTIONAL_DIR}/blast_dbs/swissprot.fasta.gz"
            makeblastdb -in "${FUNCTIONAL_DIR}/blast_dbs/swissprot.fasta" -dbtype prot -out "${BLAST_DB_SwissProt}"
        fi

        

        # Homology Search with BLASTp and BLASTx
        if [[ "${FUNCTIONAL_METHODS}" == *"BLASTp"* ]]; then
            echo_green  "Running BLASTp against SwissProt database"
            blastp -query "${LONGEST_ORFS_PEP}" \
                   -db "${BLAST_DB_SwissProt}" \
                   -outfmt 6 \
                   -num_threads "${THREADS}" \
                   -out "${FUNCTIONAL_DIR}/blastp_results.out"
        else
            echo_green  "Skipping BLASTp as per user selection."
        fi

        if [[ "${FUNCTIONAL_METHODS}" == *"BLASTx"* ]]; then
            echo_green  "Running BLASTx against SwissProt database"
            blastx -query "${TRANSCRIPTS_FA}" \
                   -db "${BLAST_DB_SwissProt}" \
                   -outfmt 6 \
                   -num_threads "${THREADS}" \
                   -out "${FUNCTIONAL_DIR}/blastx_results.out"
        else
            echo_green  "Skipping BLASTx as per user selection."
        fi

        if [[ "${FUNCTIONAL_METHODS}" == *"PFAM"* ]]; then
		    if [ -z "${PFAM_DB}" ]; then
				echo_green  "Downloading PFAM database"
				mkdir -p "${FUNCTIONAL_DIR}/pfam_db"
				PFAM_DB="${FUNCTIONAL_DIR}/pfam_db"
				wget -P "${PFAM_DB}" ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
				gunzip "${PFAM_DB}/Pfam-A.hmm.gz"
				hmmpress "${PFAM_DB}/Pfam-A.hmm"
			fi
            echo_green  "Running PFAM Scan to identify protein domains"
            hmmscan --cpu "${THREADS}" --domtblout "${FUNCTIONAL_DIR}/pfam_results.out" "${PFAM_DB}/Pfam-A.hmm" "${LONGEST_ORFS_PEP}" 
        else
            echo_green  "Skipping PFAM Scan as per user selection."
        fi

        conda deactivate

        # Filtering Criteria
        echo_green  "Filtering transcripts based on ORFs, PFAM domains, and BLAST hits"

        # Create a list of transcripts with ORFs
        grep ">" "${LONGEST_ORFS_PEP}" | sed 's/>//' > "${FUNCTIONAL_DIR}/transcripts_with_orfs.txt"

        # Initialize functional transcripts file
        FUNCTIONAL_TRANSCRIPTS="${FUNCTIONAL_DIR}/functional_transcripts.txt"
        > "${FUNCTIONAL_TRANSCRIPTS}"

        # Filter BLAST results for significant hits (e-value < 1e-5, identity > 30%)
        if [[ "${FUNCTIONAL_METHODS}" == *"BLASTp"* ]]; then
            awk '$11 < 1e-5 && $3 > 30 {print $1}' "${FUNCTIONAL_DIR}/blastp_results.out" | sort | uniq >> "${FUNCTIONAL_TRANSCRIPTS}"
        fi
        if [[ "${FUNCTIONAL_METHODS}" == *"BLASTx"* ]]; then
            awk '$11 < 1e-5 && $3 > 30 {print $1}' "${FUNCTIONAL_DIR}/blastx_results.out" | sort | uniq >> "${FUNCTIONAL_TRANSCRIPTS}"
        fi

        # Filter PFAM results for significant domains
        if [[ "${FUNCTIONAL_METHODS}" == *"PFAM"* ]]; then
            awk '{print $1}' "${FUNCTIONAL_DIR}/pfam_results.out" | sort | uniq >> "${FUNCTIONAL_TRANSCRIPTS}"
        fi

        # Remove duplicates
        sort "${FUNCTIONAL_TRANSCRIPTS}" | uniq > "${FUNCTIONAL_DIR}/functional_transcripts_sorted.txt"

        # Intersect with transcripts having ORFs
        comm -12 <(sort "${FUNCTIONAL_DIR}/transcripts_with_orfs.txt") <(sort "${FUNCTIONAL_DIR}/functional_transcripts_sorted.txt") > "${FUNCTIONAL_DIR}/final_filtered_transcripts.txt"

        # Generate final annotation GTF
        #echo_green  "Generating final filtered GTF file"
        #grep -F -f "${FUNCTIONAL_DIR}/final_filtered_transcripts.txt" "${ANNOTATED_GTF}" > "${FUNCTIONAL_DIR}/annotation_part_I.gtf"

        #echo_green  "Final filtered GTF file created at ${FUNCTIONAL_DIR}/annotation_part_I.gtf"

        echo_blue  "Step 9 Completed: Functional Annotation and Filtering"
    }
	# ---------------------------
    # Function for Step 10: Integrate Functional Annotation (R Script)
    # ---------------------------

	step10_integrate_functional_annotation() {
		echo_blue "Starting Step 10: Integrate Functional Annotation (R Script)"
		conda activate ${ENV_R}
		# Define paths
		FUNCTIONAL_SCRIPT="functional_annotation.R"  
		ANNOTATION_DIR="${OUTPUT_DIR}/annotation"
		FUNCTIONAL_DIR="${OUTPUT_DIR}/functional_annotation"
		GENOME_GTF="${GENOME_GTF}"  # Already defined earlier
		OUTPUT_DIR="${OUTPUT_DIR}"  # Already defined earlier

		# Check if R script exists
		if [ ! -f "${FUNCTIONAL_SCRIPT}" ]; then
			echo_red "Error: R script not found at ${FUNCTIONAL_SCRIPT}"
			exit 1
		fi
		
		# Set up a personal R library in your output directory
		export R_LIBS_USER="${OUTPUT_DIR}/R_lib"
		mkdir -p "${R_LIBS_USER}"


		# Run the R script with required arguments
		if [ -n "$GENOME_GTF" ]; then
		  Rscript functional_annotation.R --annotation_dir "$ANNOTATION_DIR" \
									   --functional_dir "$FUNCTIONAL_DIR" \
									   --genome_gtf "$GENOME_GTF" \
									   --output_dir "$OUTPUT_DIR" \
									   >> "$LOGS_DIR/functional_annotation.log" 2>&1
		else
		  Rscript functional_annotation.R --annotation_dir "$ANNOTATION_DIR" \
									   --functional_dir "$FUNCTIONAL_DIR" \
									   --output_dir "$OUTPUT_DIR" \
									   >> "$LOGS_DIR/functional_annotation.log" 2>&1
        fi
		EXIT_STATUS=$?
		conda deactivate
		
		FINAL_GTF="${FUNCTIONAL_DIR}/annotation_part_I.gtf"
		if [ ! -f "${FINAL_GTF}" ]; then
			echo_red "Error: Final annotated GTF file not found. Pipeline may not have completed successfully. Please check ${LOGS_DIR}/functional_annotation.log for details."
			exit 1
		fi
		# Check log file for any error keywords (case-insensitive)
		if grep -qi "error" "$LOGS_DIR/functional_annotation.log"; then
			echo_red "Error: Functional annotation R script log indicates errors. Check ${LOGS_DIR}/functional_annotation.log for details."
			exit 1
		fi
		
		# Check if R script ran successfully
		if [ $EXIT_STATUS -ne 0 ]; then
			echo_red "Error: Functional annotation R script failed with exit code $EXIT_STATUS. Check ${LOGS_DIR}/functional_annotation.log for details."
			exit 1
		else
			echo_green "Functional annotation completed successfully."
			echo_green  "Final filtered GTF file created at ${FUNCTIONAL_DIR}/annotation_part_I.gtf"
		fi

		echo_blue "Step 10 Completed: Integrate Functional Annotation (R Script)"
	}

	# Add Step 10 to the pipeline
	if [[ " ${STEPS_TO_RUN[*]} " =~ " 10 " ]]; then
		step10_integrate_functional_annotation
	fi

    # ---------------------------
    # Parse Command-Line Arguments with getopt
    # ---------------------------
    # Initialize variables with default values
    GENOME_DIR="./outputDir/genome"
    RRNA_REF=""
    RRNA_REF_INDEX=""
    GENOME_REF=""
    GENOME_GTF=""
    ALIGN_DIR=""
    finalGTF=""
    SR_RB_GTF_DIR=""
    SR_DN_GTF_DIR=""
    MR_RB_GTF_DIR=""
    MR_DN_GTF_DIR=""
    BLAST_DB_SwissProt=""
    PFAM_DB=""
    DATA_DIR=""
	functional_dir=""
    OUTPUT_DIR="./outputDir"
    THREADS=8
    MIN_ORF_LENGTH=100
    STEPS=()
    RUN_ALL=false
    FUNCTIONAL_METHODS="BLASTp,BLASTx,PFAM"  # Default: all methods

    # Use getopt for parsing long options
    if ! PARSED_OPTIONS=$(getopt -n "$0" -o "" --long genomeDir:,rrnaRef:,genomeRef:,genomeGTF:,ALIGN_DIR:,finalGTF:,SR_RB_gtf_dir:,SR_DN_gtf_dir:,MR_RB_gtf_dir:,MR_DN_gtf_dir:,blastDB_SwissProt:,functional_dir:,pfamDB:,dataDir:,outputDir:,threads:,minOrfLength:,maxExonLength:,maxTranscriptLength:,steps:,all,help,functionalMethods:,conda:,conda_short:,conda_mix: -- "$@"); then
		show_help
	fi
	eval set -- "$PARSED_OPTIONS"

	while [[ "$#" -gt 0 ]]; do
		case $1 in
			--genomeDir) GENOME_DIR="$2"; shift 2 ;;
			--rrnaRef) RRNA_REF="$2"; shift 2 ;;
			--genomeRef) GENOME_REF="$2"; shift 2 ;;
			--genomeGTF) GENOME_GTF="$2"; shift 2 ;;
			--ALIGN_DIR) ALIGN_DIR="$2"; shift 2 ;;
			--finalGTF) finalGTF="$2"; shift 2 ;;
			--SR_RB_gtf_dir) SR_RB_GTF_DIR="$2"; shift 2 ;;
			--SR_DN_gtf_dir) SR_DN_GTF_DIR="$2"; shift 2 ;;
			--MR_RB_gtf_dir) MR_RB_GTF_DIR="$2"; shift 2 ;;
			--MR_DN_gtf_dir) MR_DN_GTF_DIR="$2"; shift 2 ;;
			--blastDB_SwissProt) BLAST_DB_SwissProt="$2"; shift 2 ;;
			--pfamDB) PFAM_DB="$2"; shift 2 ;;
			--dataDir) DATA_DIR="$2"; shift 2 ;;
			--functional_dir) functional_dir="$2"; shift 2 ;;
			--outputDir) OUTPUT_DIR="$2"; shift 2 ;;
			--threads) THREADS="$2"; shift 2 ;;
			--minOrfLength) MIN_ORF_LENGTH="$2"; shift 2 ;;
			--maxExonLength) MAX_EXON_LENGTH="$2"; shift 2 ;;
			--maxTranscriptLength) MAX_TRANSCRIPT_LENGTH="$2"; shift 2 ;;
			--steps) IFS=',' read -ra STEPS <<< "$2"; shift 2 ;;
			--all) RUN_ALL=true; shift ;;
			--functionalMethods) FUNCTIONAL_METHODS="$2"; shift 2 ;;
			--conda) CONDA_VERSION="$2"; shift 2 ;;
			--conda_short) CONDA_SHORT_VERSION="$2"; shift 2 ;;
			--conda_mix) CONDA_MIX_VERSION="$2"; shift 2 ;;
			-h|--help) show_help ;;
			--) shift; break ;;
			*) echo "Unknown parameter passed: $1"; show_help ;;
		esac
	done
	
	# If step 2 is requested
	if contains_step "2" "${STEPS_TO_RUN[@]}"; then
		if [ -z "${GENOME_DIR}" ]; then
			GENOME_DIR="${OUTPUT_DIR}/genome"
			echo_green "Step 2 requires a genome index. GENOME_DIR not provided; auto-setting GENOME_DIR to ${GENOME_DIR}"
		fi
	fi

	# If --conda was provided and the specific versions are still at their defaults, override them
	if [ -n "$CONDA_VERSION" ]; then
		if [ "$CONDA_SHORT_VERSION" = "$DEFAULT_CONDA_SHORT" ]; then
			CONDA_SHORT_VERSION="$CONDA_VERSION"
		fi
		if [ "$CONDA_MIX_VERSION" = "$DEFAULT_CONDA_MIX" ]; then
			CONDA_MIX_VERSION="$CONDA_VERSION"
		fi
	fi
	
	if [ "$CONDA_SHORT_VERSION" = "$CONDA_MIX_VERSION" ]; then
		ENV_COMMON="stringtie_${CONDA_SHORT_VERSION//./}"
		ENV_SHORT="${ENV_COMMON}"
		ENV_MIX="${ENV_COMMON}"
	else
		ENV_SHORT="stringtie_short_${CONDA_SHORT_VERSION//./}"
		ENV_MIX="stringtie_mix_${CONDA_MIX_VERSION//./}"
	fi
	
	# Validate that MAX_EXON_LENGTH and MAX_TRANSCRIPT_LENGTH are positive integers
	if ! [[ "$MAX_EXON_LENGTH" =~ ^[1-9][0-9]*$ ]]; then
		echo_red "Error: --maxExonLength must be a positive integer."
		exit 1
	fi
	if ! [[ "$MAX_TRANSCRIPT_LENGTH" =~ ^[1-9][0-9]*$ ]]; then
		echo_red "Error: --maxTranscriptLength must be a positive integer."
		exit 1
	fi

    # ---------------------------
    # Validate Mandatory Arguments Based on Selected Steps
    # ---------------------------
    if [ "$RUN_ALL" = false ] && [ ${#STEPS[@]} -eq 0 ]; then
        echo_red  "Error: You must specify --all or provide --steps to run."
        show_help
    fi

    # Validate Functional Methods
    IFS=',' read -ra VALID_METHODS <<< "BLASTp,BLASTx,PFAM"
    for method in $(echo "${FUNCTIONAL_METHODS}" | tr ',' ' '); do
        if [[ ! " ${VALID_METHODS[*]} " =~ " ${method} " ]]; then
            echo_red "Error: Invalid functional annotation method '${method}'. Valid options are BLASTp, BLASTx, PFAM."
            exit 1
        fi
    done

    # ---------------------------
    # Set Up Output Directories
    # ---------------------------
	if [ "$RUN_ALL" = true ]; then
		GENOME_DIR="${OUTPUT_DIR}/genome"
		PREPROC_DIR="${OUTPUT_DIR}/preprocessing"
		ALIGN_DIR="${OUTPUT_DIR}/alignment"
		ASSEMBLY_DIR="${OUTPUT_DIR}/assembly"
		MERGE_DIR="${OUTPUT_DIR}/merging"
		ANNOTATION_DIR="${OUTPUT_DIR}/annotation"
		FUNCTIONAL_DIR="${OUTPUT_DIR}/functional_annotation"
		LOGS_DIR="${OUTPUT_DIR}/logs"
		SR_RB_GTF_DIR="${ASSEMBLY_DIR}/SR_RB"
		SR_DN_GTF_DIR="${ASSEMBLY_DIR}/SR_DN"
		MR_RB_GTF_DIR="${ASSEMBLY_DIR}/MR_RB"
		MR_DN_GTF_DIR="${ASSEMBLY_DIR}/MR_DN"
		
		#mkdir -p "${PREPROC_DIR}" "${ALIGN_DIR}" "${ASSEMBLY_DIR}" "${MERGE_DIR}" "${ANNOTATION_DIR}" "${FUNCTIONAL_DIR}" "${LOGS_DIR}" "${SR_RB_GTF_DIR}" "${SR_DN_GTF_DIR}" "${MR_RB_GTF_DIR}" "${MR_DN_GTF_DIR}"
	fi

    # ---------------------------
	# Set up default output directories if not provided
	# ---------------------------
	if [ -z "${OUTPUT_DIR}" ]; then
		OUTPUT_DIR="./outputDir"
	fi
	if [ -z "${GENOME_DIR}" ]; then
		GENOME_DIR="${OUTPUT_DIR}/genome"
	fi
	if [ -z "${PREPROC_DIR}" ]; then
		PREPROC_DIR="${OUTPUT_DIR}/preprocessing"
	fi
	if [ -z "${ALIGN_DIR}" ]; then
		ALIGN_DIR="${OUTPUT_DIR}/alignment"
	fi
	if [ -z "${ASSEMBLY_DIR}" ]; then
		ASSEMBLY_DIR="${OUTPUT_DIR}/assembly"
	fi
	if [ -z "${MERGE_DIR}" ]; then
		MERGE_DIR="${OUTPUT_DIR}/merging"
	fi
	if [ -z "${ANNOTATION_DIR}" ]; then
		ANNOTATION_DIR="${OUTPUT_DIR}/annotation"
	fi
	if [ -z "${FUNCTIONAL_DIR}" ]; then
		FUNCTIONAL_DIR="${OUTPUT_DIR}/functional_annotation"
	fi
	if [ -z "${LOGS_DIR}" ]; then
		LOGS_DIR="${OUTPUT_DIR}/logs"
	fi
	if [ -z "${SR_RB_GTF_DIR}" ]; then
		SR_RB_GTF_DIR="${ASSEMBLY_DIR}/SR_RB"
	fi
	if [ -z "${SR_DN_GTF_DIR}" ]; then
		SR_DN_GTF_DIR="${ASSEMBLY_DIR}/SR_DN"
	fi
	if [ -z "${MR_RB_GTF_DIR}" ]; then
		MR_RB_GTF_DIR="${ASSEMBLY_DIR}/MR_RB"
	fi
	if [ -z "${MR_DN_GTF_DIR}" ]; then
		MR_DN_GTF_DIR="${ASSEMBLY_DIR}/MR_DN"
	fi

    # Create directories
    mkdir -p "${OUTPUT_DIR}" "${PREPROC_DIR}" "${ALIGN_DIR}" "${ASSEMBLY_DIR}" "${MERGE_DIR}" "${ANNOTATION_DIR}" "${FUNCTIONAL_DIR}" "${LOGS_DIR}" "${SR_RB_GTF_DIR}" "${SR_DN_GTF_DIR}" "${MR_RB_GTF_DIR}" "${MR_DN_GTF_DIR}"


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

    if [ "${#SHORT_ONLY_SAMPLES[@]}" -gt 0 ]; then
    echo_green "Short-only Samples detected: ${SHORT_ONLY_SAMPLES[*]}"
	fi
	if [ "${#MIX_SAMPLES[@]}" -gt 0 ]; then
    echo_green "Mixed Samples detected: ${MIX_SAMPLES[*]}"
	fi
	


    # ---------------------------
    # Check Dependencies and Set Up Environments
    # ---------------------------
    check_dependencies
    setup_conda_environments
	# ---------------------------
	# Define Step Order
	# ---------------------------
	STEP_ORDER=(1 2 3 4 5 6 7 8 9 10)
	
    # ---------------------------
	# Determine Steps to Run
	# ---------------------------
	if [ "$RUN_ALL" = true ]; then
		STEPS_TO_RUN=("${STEP_ORDER[@]}")
	else
		if [ ${#STEPS[@]} -eq 0 ]; then
			echo_red "Error: You must specify --all or provide --steps to run."
			show_help
		fi
		for STEP in "${STEPS[@]}"; do
			if [[ ! " ${STEP_ORDER[*]} " =~ " ${STEP} " ]]; then
				echo_red "Error: Invalid step number '${STEP}'. Must be between 1 and 10."
				exit 1
			fi
		done
		STEPS_TO_RUN=("${STEPS[@]}")
	fi

	# ---------------------------
	# Define Step Functions Mapping
	# ---------------------------
	declare -A STEP_FUNCTIONS=(
		[1]=step1_rrna_removal
		[2]=step2_read_alignment
		[3]=step3_gene_transcript_assembly
		[4]=step4_merging_transcripts_I
		[5]=step5_merging_transcripts_II
		[6]=step6_filter_transcripts
		[7]=step7_isoform_comparison
		[8]=step8_gtf_correction
		[9]=step9_functional_annotation
		[10]=step10_integrate_functional_annotation
	)

	

	# ---------------------------
	# Run Selected Steps in Order
	# ---------------------------
	for STEP in "${STEP_ORDER[@]}"; do
		if contains_step "$STEP" "${STEPS_TO_RUN[@]}"; then
			STEP_FUNCTION="${STEP_FUNCTIONS[$STEP]}"
			if [ -z "${STEP_FUNCTION}" ]; then
				echo_red "Error: No function defined for Step $STEP."
				exit 1
			fi
			echo_green "Executing Step $STEP: ${STEP_FUNCTION}"
			"${STEP_FUNCTION}"
		fi
	done
	



    # ---------------------------
    # Validate Conditional Mandatory Options Based on Selected Steps
    # ---------------------------
    for STEP in "${STEPS_TO_RUN[@]}"; do
        case "$STEP" in
            1|2)
                #if [ -z "${genomeRef}" ]; then
                #    echo_red "Error: --genomeRef must be provided when running Step 8 (Functional Annotation and Filtering)."
                #    exit 1
                #fi
				if [ -z "${DATA_DIR}" ]; then
					echo_red  "Error: --dataDir is required."
					show_help
				fi


				if [ ! -d "${DATA_DIR}/short_reads" ] && [ ! -d "${DATA_DIR}/mix_reads" ]; then
					echo_red  "Error: --dataDir must contain 'short_reads' and/or 'mix_reads' directories."
					show_help
				fi
                ;;
            3)
			    if [ "$RUN_ALL" = false ]; then
					if [ -z "${ALIGN_DIR}" ] || [ ! -d "${ALIGN_DIR}" ]; then
						echo_red "Error: --alignDir must be provided and must exist when running Step 3 (Gene and Transcript Assembly)."
						exit 1
					fi
				fi
                ;;
            4|5|6)
			    if [ "$RUN_ALL" = false ]; then
					# For merging assemblies, require specific GTF directories
					if [ -z "${SR_RB_GTF_DIR}" ] && [ -z "${SR_DN_GTF_DIR}" ] && [ -z "${MR_RB_GTF_DIR}" ] && [ -z "${MR_DN_GTF_DIR}" ]; then
						echo_red "Error: At least one of --SR_RB_gtf_dir, --SR_DN_gtf_dir, --MR_RB_gtf_dir, or --MR_DN_gtf_dir must be provided when running Steps 4, 5, or 6 (Merging Assemblies)."
						exit 1
					fi
				fi
                ;;
            7|8)
			    if [ "$RUN_ALL" = false ]; then
					if [ -z "${finalGTF}" ] && [ ! -f "${MERGE_DIR}/prefinal_annotation.gtf" ]; then
						echo_red "Error: --finalGTF must be provided or prefinal_annotation.gtf must exist when running Steps 7 or 8 (Isoform Comparison and Annotation, GTF Correction)."
						exit 1
					fi
				fi
                ;;
            9)  
			    if [ "$RUN_ALL" = false ]; then
					if [ -z "${finalGTF}" ] && [ ! -f "${ANNOTATION_DIR}/corrected_with_introns.gtf" ] && [ ! -f "${MERGE_DIR}/filtered_annotation.gtf" ]; then
						echo_red "Error: --finalGTF must be provided when running Step 9 (Functional Annotation and Filtering)."
						exit 1
					fi
				fi
				#if [ -z "${genomeRef}" ] && [ ! -f "${GENOME_REF}" ]; then
                #    echo_red "Error: --genomeRef must be provided when running Step 9 (Functional Annotation and Filtering)."
                #    exit 1
                #fi
				
                ;;
			10)
			    if [ "$RUN_ALL" = false ]; then
                    if [ -z "${outputDir}" ]; then
					echo_red "Error: --outputDir is required for Step 10 (Integrate Functional Annotation)."
					exit 1
				fi
				# Then check that at least one pre-final annotation GTF exists.
				if [ ! -f "${MERGE_DIR}/filtered_annotation.gtf" ] && [ ! -f "${ANNOTATION_DIR}/corrected_with_introns.gtf" ]; then
					echo_red "Error: Pre-final annotation GTF is missing. Ensure Steps 1-9 have been successfully completed before running Step 10."
					exit 1
				fi
                fi
                    ;;
            *)
                ;;
        esac
    done

    # ---------------------------
    # Run Selected Steps in Order
    # ---------------------------
    for STEP in "${STEP_ORDER[@]}"; do
        if [[ " ${STEPS_TO_RUN[*]} " =~ " ${STEP} " ]]; then
            STEP_FUNCTION="${STEP_FUNCTIONS[$STEP]}"
            if [ -z "${STEP_FUNCTION}" ]; then
                echo_red  "Error: No function defined for Step $STEP."
                exit 1
            fi
            echo_green  "Executing Step $STEP: ${STEP_FUNCTION}"
            "${STEP_FUNCTION}"
        fi
    done

    # =====================================================================
    # Pipeline Completed
    # =====================================================================

    echo_green  "RNA-Seq Bioinformatics Pipeline Completed Successfully!"
    if [ -f "${FUNCTIONAL_DIR}/annotation_part_I.gtf" ]; then
        echo_green  "Final annotated GTF file is located at ${FUNCTIONAL_DIR}/annotation_part_I.gtf"
    elif [ -f "${ANNOTATION_DIR}/corrected_with_introns.gtf" ]; then
        echo_green  "Final annotated GTF file is located at ${ANNOTATION_DIR}/corrected_with_introns.gtf"
    fi
