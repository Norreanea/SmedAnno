#!/bin/bash

# =====================================================================
# Enhanced RNA-Seq Bioinformatics Pipeline
# =====================================================================

set -e  # Exit immediately if a command exits with a non-zero status
set -o pipefail  # Pipeline returns the exit status of the last command to fail

# ---------------------------
# Function to Display Help
# ---------------------------
show_help() {
    echo "Usage: $0 [OPTIONS]"
    echo ""
    echo "Options:"
    echo "  --genomeDir PATH          Path to STAR genome directory"
    echo "  --rrnaRef PATH            Path to rRNA reference FASTA file"
    echo "  --genomeRef PATH          Path to genome reference FASTA file"
    echo "  --genomeGTF PATH          Path to genome annotation GTF file"
    echo "  --blastDB_NR PATH         Path to BLAST NR database"
    echo "  --blastDB_SwissProt PATH  Path to BLAST SwissProt database"
    echo "  --pfamDB PATH             Path to PFAM database directory"
    echo "  --dataDir PATH            Path to input RNA-Seq data directory"
    echo "  --outputDir PATH          Path to output directory"
    echo "  --threads N               Number of CPU threads to use (default: 8)"
    echo "  --steps LIST              Comma-separated list of steps to run (2-9)"
    echo "  --all                     Run all steps"
    echo "  --help                    Display this help message and exit"
    echo ""
    echo "Example:"
    echo "  $0 --genomeDir /path/to/genomeDir --rrnaRef /path/to/rrna.fa --all"
    exit 1
}

# ---------------------------
# Function to Check Dependencies
# ---------------------------
check_dependencies() {
    local dependencies=("STAR" "Minimap2" "StringTie" "gffcompare" "agat_convert_sp_gxf2gxf.pl" "agat_sp_add_introns.pl" "TransDecoder.LongOrfs" "TransDecoder.Predict" "blastp" "blastx" "pfam_scan.pl" "gffread" "conda")
    for cmd in "${dependencies[@]}"; do
        if ! command -v "$cmd" &> /dev/null; then
            echo "Error: $cmd is not installed or not in PATH."
            exit 1
        fi
    done
}

# ---------------------------
# Function for Step 2: Preprocessing - rRNA Removal
# ---------------------------
step2_rrna_removal() {
    echo "Starting Step 2: rRNA Removal"

    for SAMPLE in "${SAMPLES[@]}"; do
        echo "Processing Sample: ${SAMPLE}"

        # Paths to input FASTQ files
        READ1="${DATA_DIR}/${SAMPLE}_R1.fastq.gz"
        READ2="${DATA_DIR}/${SAMPLE}_R2.fastq.gz"

        # Output files
        NON_RRNA_READ1="${PREPROC_DIR}/${SAMPLE}_non_rrna_R1.fastq.gz"
        NON_RRNA_READ2="${PREPROC_DIR}/${SAMPLE}_non_rrna_R2.fastq.gz"

        # Short Reads - STAR Aligner
        echo "Removing rRNA from short reads using STAR"
        STAR --runThreadN ${THREADS} \
             --genomeDir "${RRNA_REF_INDEX}" \
             --readFilesIn "${READ1}" "${READ2}" \
             --readFilesCommand zcat \
             --outFilterType BySJout \
             --outFilterMultimapNmax 1 \
             --outFilterMismatchNmax 0 \
             --outReadsUnmapped Fastx \
             --outFileNamePrefix "${PREPROC_DIR}/${SAMPLE}_STAR_rrna_"

        # Extract non-rRNA reads
        zcat "${PREPROC_DIR}/${SAMPLE}_STAR_rrna_Unmapped.out.mate1.gz" > "${NON_RRNA_READ1}"
        zcat "${PREPROC_DIR}/${SAMPLE}_STAR_rrna_Unmapped.out.mate2.gz" > "${NON_RRNA_READ2}"

        # Long Reads - Minimap2 (if applicable)
        if [ -n "${LONG_READS_DIR}" ]; then
            LONG_READS="${LONG_READS_DIR}/${SAMPLE}_long.fastq.gz"
            NON_RRNA_LONG="${PREPROC_DIR}/${SAMPLE}_non_rrna_long.fastq.gz"
            echo "Removing rRNA from long reads using Minimap2"
            minimap2 -ax sr "${RRNA_REF}" "${LONG_READS}" | grep -v "^@" | awk '{if($3=="*") print $0}' | gzip > "${NON_RRNA_LONG}"
        fi

        echo "rRNA removal completed for Sample: ${SAMPLE}"
    done

    echo "Step 2 Completed: rRNA Removal"
}

# ---------------------------
# Function for Step 3: Read Alignment to Reference Genome
# ---------------------------
step3_read_alignment() {
    echo "Starting Step 3: Read Alignment to the Reference Genome"

    for SAMPLE in "${SAMPLES[@]}"; do
        echo "Aligning Sample: ${SAMPLE}"

        NON_RRNA_READ1="${PREPROC_DIR}/${SAMPLE}_non_rrna_R1.fastq.gz"
        NON_RRNA_READ2="${PREPROC_DIR}/${SAMPLE}_non_rrna_R2.fastq.gz"

        # Align reads using STAR
        echo "Running STAR aligner for Sample: ${SAMPLE}"
        STAR --runThreadN ${THREADS} \
             --genomeDir "${GENOME_DIR}" \
             --readFilesIn "${NON_RRNA_READ1}" "${NON_RRNA_READ2}" \
             --readFilesCommand zcat \
             --outSAMtype BAM SortedByCoordinate \
             --outFileNamePrefix "${ALIGN_DIR}/${SAMPLE}_"

        # The output BAM file
        ALIGNED_BAM="${ALIGN_DIR}/${SAMPLE}_Aligned.sortedByCoord.out.bam"

        echo "Alignment completed for Sample: ${SAMPLE}"
    done

    echo "Step 3 Completed: Read Alignment"
}

# ---------------------------
# Function for Step 4: Gene and Transcript Assembly
# ---------------------------
step4_gene_transcript_assembly() {
    echo "Starting Step 4: Gene and Transcript Assembly"

    for SAMPLE in "${SAMPLES[@]}"; do
        echo "Assembling transcripts for Sample: ${SAMPLE}"

        ALIGNED_BAM="${ALIGN_DIR}/${SAMPLE}_Aligned.sortedByCoord.out.bam"

        # Reference-Based Annotation (RB) with StringTie2 v2.1.1
        echo "Running StringTie2 v2.1.1 for Reference-Based Assembly (RB)"
        source activate stringtie211
        stringtie -p ${THREADS} \
                  -G "${GENOME_GTF}" \
                  -c 1.5 \
                  -f 0.02 \
                  -o "${ASSEMBLY_DIR}/${SAMPLE}_RB.gtf" \
                  "${ALIGNED_BAM}"
        source deactivate

        # Reference-Based Annotation (RB) with StringTie2 v2.2.1 (Mixed Reads)
        echo "Running StringTie2 v2.2.1 for Reference-Based Assembly (Mixed Reads)"
        source activate stringtie221
        stringtie --mix -p ${THREADS} \
                  -G "${GENOME_GTF}" \
                  -c 1.5 \
                  -f 0.02 \
                  -o "${ASSEMBLY_DIR}/${SAMPLE}_mix_RB.gtf" \
                  "${ALIGNED_BAM}"
        source deactivate

        # De Novo Annotation (DN) with StringTie2 v2.1.1
        echo "Running StringTie2 v2.1.1 for De Novo Assembly (DN)"
        source activate stringtie211
        stringtie -p ${THREADS} \
                  -c 1.5 \
                  -f 0.02 \
                  -o "${ASSEMBLY_DIR}/${SAMPLE}_DN.gtf" \
                  "${ALIGNED_BAM}"
        source deactivate

        # De Novo Annotation (DN) with StringTie2 v2.2.1 (Mixed Reads)
        echo "Running StringTie2 v2.2.1 for De Novo Assembly (Mixed Reads)"
        source activate stringtie221
        stringtie --mix -p ${THREADS} \
                  -c 1.5 \
                  -f 0.02 \
                  -o "${ASSEMBLY_DIR}/${SAMPLE}_mix_DN.gtf" \
                  "${ALIGNED_BAM}"
        source deactivate

        echo "Transcript assembly completed for Sample: ${SAMPLE}"
    done

    echo "Step 4 Completed: Gene and Transcript Assembly"
}

# ---------------------------
# Function for Step 5: Merge Annotation from Step 1
# ---------------------------
step5_merging_transcripts_step1() {
    echo "Starting Step 5: Merge Annotation from Step 1"

    # Activate StringTie 2.2.1 environment
    source activate stringtie221

    # Create a list of GTF files from Step 1 (Reference-Based, Short Reads)
    echo "Creating list of GTF files from Step 1 (RB, Short Reads)"
    > "${LIST_STEP1_GTF}"  # Empty the file
    for SAMPLE in "${SAMPLES[@]}"; do
        echo "${ASSEMBLY_DIR}/${SAMPLE}_RB.gtf" >> "${LIST_STEP1_GTF}"
    done

    # Merge annotations from Step 1
    echo "Merging Reference-Based Assemblies (Step 1)"
    stringtie --merge -p ${THREADS} \
             -G "${REFERENCE_GTF}" \
             -o "${MERGE_DIR}/merged_step5.gtf" \
             "${LIST_STEP1_GTF}"

    # Deactivate the environment
    source deactivate

    echo "Step 5 Completed: Merge Annotation from Step 1"
}

# ---------------------------
# Function for Step 6: Merge Annotation from Step 5 and Step 2
# ---------------------------
step6_merging_transcripts_step2() {
    echo "Starting Step 6: Merge Annotation from Step 5 and Step 2"

    # Activate StringTie 2.2.1 environment
    source activate stringtie221

    # Create a list of GTF files from Step 2 (Reference-Based, Mixed Reads)
    echo "Creating list of GTF files from Step 2 (RB, Mixed Reads)"
    > "${LIST_STEP2_GTF}"  # Empty the file
    for SAMPLE in "${SAMPLES[@]}"; do
        echo "${ASSEMBLY_DIR}/${SAMPLE}_mix_RB.gtf" >> "${LIST_STEP2_GTF}"
    done

    # Merge annotations from Step 5 and Step 2
    echo "Merging Step 5 and Step 2 Annotations"
    stringtie --merge -p ${THREADS} \
             -G "${MERGE_DIR}/merged_step5.gtf" \
             -o "${MERGE_DIR}/merged_step6.gtf" \
             "${LIST_STEP2_GTF}"

    # Deactivate the environment
    source deactivate

    echo "Step 6 Completed: Merge Annotation from Step 5 and Step 2"
}

# ---------------------------
# Function for Step 7: Isoform Comparison and Annotation
# ---------------------------
step7_isoform_comparison() {
    echo "Starting Step 7: Isoform Comparison and Annotation"

    COMP_OUTPUT_PREFIX="${ANNOTATION_DIR}/gffcomp"

    echo "Running gffcompare to compare assembled transcripts with reference"
    gffcompare -r "${GENOME_GTF}" \
               -G \
               --chr-stats \
               --strict-match \
               -o "${COMP_OUTPUT_PREFIX}" \
               "${MERGE_DIR}/merged_step6.gtf"

    echo "Step 7 Completed: Isoform Comparison and Annotation"
}

# ---------------------------
# Function for Step 8: GTF File Correction and Enhancement
# ---------------------------
step8_gtf_correction() {
    echo "Starting Step 8: GTF Correction and Enhancement"

    ANNOTATED_GTF="${ANNOTATION_DIR}/gffcomp.annotated.gtf"
    CORRECTED_GTF="${ANNOTATION_DIR}/corrected.gtf"
    CORRECTED_INTRONS_GTF="${ANNOTATION_DIR}/corrected_with_introns.gtf"

    echo "Converting and correcting GTF file using AGAT"
    agat_convert_sp_gxf2gxf.pl -i "${ANNOTATED_GTF}" -o "${CORRECTED_GTF}"

    echo "Adding intron features to GTF file using AGAT"
    agat_sp_add_introns.pl -i "${CORRECTED_GTF}" -o "${CORRECTED_INTRONS_GTF}"

    echo "Step 8 Completed: GTF Correction and Enhancement"
}

# ---------------------------
# Function for Step 9: Functional Annotation and Filtering
# ---------------------------
step9_functional_annotation() {
    echo "Starting Step 9: Functional Annotation and Filtering"

    # Paths to transcript and protein files
    TRANSCRIPTS_FA="${ANNOTATION_DIR}/transcripts.fa"
    LONGEST_ORFS_PEP="${FUNCTIONAL_DIR}/longest_orfs.pep"
    LONGEST_ORFS_FA="${FUNCTIONAL_DIR}/longest_orfs.fa"

    # Extract transcript sequences (requires gffread)
    echo "Extracting transcript sequences using gffread"
    gffread "${ANNOTATION_DIR}/corrected_with_introns.gtf" -g "${GENOME_REF}" -w "${TRANSCRIPTS_FA}"

    # ORF Prediction with TransDecoder
    echo "Running TransDecoder to predict ORFs"
    TransDecoder.LongOrfs -t "${TRANSCRIPTS_FA}" -m ${MIN_ORF_LENGTH}

    # Predict likely coding regions (optional)
    echo "Running TransDecoder.Predict to identify likely coding regions"
    TransDecoder.Predict -t "${TRANSCRIPTS_FA}" --single_best_orf -o "${FUNCTIONAL_DIR}/TransDecoder"

    # Extract longest ORFs
    echo "Extracting longest ORFs"
    cp "${FUNCTIONAL_DIR}/TransDecoder/longest_orfs.pep" "${LONGEST_ORFS_PEP}"
    cp "${FUNCTIONAL_DIR}/TransDecoder/longest_orfs.fa" "${LONGEST_ORFS_FA}"

    # Homology Search with BLASTp and BLASTx
    echo "Running BLASTp against NR database"
    blastp -query "${LONGEST_ORFS_PEP}" \
           -db "${BLAST_DB_NR}" \
           -outfmt 6 \
           -num_threads ${THREADS} \
           -out "${FUNCTIONAL_DIR}/blastp_results.out"

    echo "Running BLASTx against SwissProt database"
    blastx -query "${TRANSCRIPTS_FA}" \
           -db "${BLAST_DB_SwissProt}" \
           -outfmt 6 \
           -num_threads ${THREADS} \
           -out "${FUNCTIONAL_DIR}/blastx_results.out"

    # Domain Identification with PFAM Scan
    echo "Running PFAM Scan to identify protein domains"
    pfam_scan.pl -fasta "${LONGEST_ORFS_PEP}" \
                -dir "${PFAM_DB}" \
                -outfile "${FUNCTIONAL_DIR}/pfam_results.out"

    # Filtering Criteria
    echo "Filtering transcripts based on ORFs, PFAM domains, and BLAST hits"

    # Create a list of transcripts with ORFs
    grep ">" "${TRANSCRIPTS_FA}" | sed 's/>//' > "${FUNCTIONAL_DIR}/transcripts_with_orfs.txt"

    # Filter BLAST results for significant hits (e.g., e-value < 1e-5, identity > 30%)
    awk '$11 < 1e-5 && $3 > 30 {print $1}' "${FUNCTIONAL_DIR}/blastp_results.out" | sort | uniq > "${FUNCTIONAL_DIR}/blastp_filtered.txt"
    awk '$11 < 1e-5 && $3 > 30 {print $1}' "${FUNCTIONAL_DIR}/blastx_results.out" | sort | uniq > "${FUNCTIONAL_DIR}/blastx_filtered.txt"

    # Filter PFAM results for significant domains (assuming presence indicates functional)
    awk '{print $1}' "${FUNCTIONAL_DIR}/pfam_results.out" | sort | uniq > "${FUNCTIONAL_DIR}/pfam_filtered.txt"

    # Combine filters
    cat "${FUNCTIONAL_DIR}/blastp_filtered.txt" "${FUNCTIONAL_DIR}/blastx_filtered.txt" "${FUNCTIONAL_DIR}/pfam_filtered.txt" | sort | uniq > "${FUNCTIONAL_DIR}/functional_transcripts.txt"

    # Intersect with transcripts having ORFs
    comm -12 <(sort "${FUNCTIONAL_DIR}/transcripts_with_orfs.txt") <(sort "${FUNCTIONAL_DIR}/functional_transcripts.txt") > "${FUNCTIONAL_DIR}/final_filtered_transcripts.txt"

    # Generate final annotation GTF
    echo "Generating final filtered GTF file"
    grep -F -f "${FUNCTIONAL_DIR}/final_filtered_transcripts.txt" "${ANNOTATION_DIR}/corrected_with_introns.gtf" > "${FUNCTIONAL_DIR}/final_annotation.gtf"

    echo "Step 9 Completed: Functional Annotation and Filtering"
}

# ---------------------------
# Parse Command-Line Arguments
# ---------------------------
# Initialize variables with default values
GENOME_DIR=""
RRNA_REF=""
RRNA_REF_INDEX=""
GENOME_REF=""
GENOME_GTF=""
BLAST_DB_NR=""
BLAST_DB_SwissProt=""
PFAM_DB=""
DATA_DIR=""
OUTPUT_DIR=""
THREADS=8
STEPS=()
RUN_ALL=false
LONG_READS_DIR=""

# Use getopt for parsing long options
PARSED_OPTIONS=$(getopt -n "$0" -o "" --long genomeDir:,rrnaRef:,genomeRef:,genomeGTF:,blastDB_NR:,blastDB_SwissProt:,pfamDB:,dataDir:,outputDir:,threads:,steps:,all,help -- "$@")
if [ $? -ne 0 ]; then
    show_help
fi

eval set -- "$PARSED_OPTIONS"

while true; do
    case "$1" in
        --genomeDir)
            GENOME_DIR="$2"
            shift 2
            ;;
        --rrnaRef)
            RRNA_REF="$2"
            shift 2
            ;;
        --genomeRef)
            GENOME_REF="$2"
            shift 2
            ;;
        --genomeGTF)
            GENOME_GTF="$2"
            shift 2
            ;;
        --blastDB_NR)
            BLAST_DB_NR="$2"
            shift 2
            ;;
        --blastDB_SwissProt)
            BLAST_DB_SwissProt="$2"
            shift 2
            ;;
        --pfamDB)
            PFAM_DB="$2"
            shift 2
            ;;
        --dataDir)
            DATA_DIR="$2"
            shift 2
            ;;
        --outputDir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        --threads)
            THREADS="$2"
            shift 2
            ;;
        --steps)
            IFS=',' read -ra STEPS <<< "$2"
            shift 2
            ;;
        --all)
            RUN_ALL=true
            shift
            ;;
        --help)
            show_help
            ;;
        --)
            shift
            break
            ;;
        *)
            echo "Unknown option: $1"
            show_help
            ;;
    esac
done

# ---------------------------
# Validate Required Arguments
# ---------------------------
if [ "$RUN_ALL" = false ] && [ ${#STEPS[@]} -eq 0 ]; then
    echo "Error: You must specify --all or provide --steps to run."
    show_help
fi

# Check required arguments for all steps
REQUIRED_ARGS="genomeDir rrnaRef genomeRef genomeGTF blastDB_NR blastDB_SwissProt pfamDB dataDir outputDir"
for ARG in $REQUIRED_ARGS; do
    eval "VALUE=\$$ARG"
    if [ -z "$VALUE" ]; then
        echo "Error: --${ARG} is required."
        show_help
    fi
done

# Optional: Long reads directory (if you have long reads)
# Uncomment and set the LONG_READS_DIR variable if applicable
# LONG_READS_DIR="/path/to/long_reads"

# ---------------------------
# Set Up Output Directories
# ---------------------------
PREPROC_DIR="${OUTPUT_DIR}/preprocessing"
ALIGN_DIR="${OUTPUT_DIR}/alignment"
ASSEMBLY_DIR="${OUTPUT_DIR}/assembly"
MERGE_DIR="${OUTPUT_DIR}/merging"
ANNOTATION_DIR="${OUTPUT_DIR}/annotation"
FUNCTIONAL_DIR="${OUTPUT_DIR}/functional_annotation"

mkdir -p "${PREPROC_DIR}" "${ALIGN_DIR}" "${ASSEMBLY_DIR}" "${MERGE_DIR}" "${ANNOTATION_DIR}" "${FUNCTIONAL_DIR}"

# ---------------------------
# Set rRNA Reference Index (STAR requires genomeDir)
# ---------------------------
# Assuming the user provides the STAR index for rRNA as part of genomeDir
# Alternatively, you may need to build it here if not provided
# For simplicity, we'll assume it's provided
RRNA_REF_INDEX="${GENOME_DIR}/rRNA_STAR_index"

# ---------------------------
# Define Sample Names
# ---------------------------
# Extract sample names based on FASTQ files in dataDir
SAMPLES=($(ls "${DATA_DIR}"/*_R1.fastq.gz | sed 's/_R1.fastq.gz//g' | xargs -n 1 basename))

echo "Samples detected: ${SAMPLES[@]}"

# ---------------------------
# Check Dependencies
# ---------------------------
check_dependencies

# ---------------------------
# Run Selected Steps
# ---------------------------
if [ "$RUN_ALL" = true ]; then
    step2_rrna_removal
    step3_read_alignment
    step4_gene_transcript_assembly
    step5_merging_transcripts_step1
    step6_merging_transcripts_step2
    step7_isoform_comparison
    step8_gtf_correction
    step9_functional_annotation
else
    for STEP in "${STEPS[@]}"; do
        case "$STEP" in
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
                step5_merging_transcripts_step1
                ;;
            6)
                step6_merging_transcripts_step2
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
                echo "Warning: Step $STEP is not recognized and will be skipped."
                ;;
        esac
    done
fi

# =====================================================================
# Pipeline Completed
# =====================================================================

echo "RNA-Seq Bioinformatics Pipeline Completed Successfully!"
echo "Final annotated GTF file is located at ${FUNCTIONAL_DIR}/final_annotation.gtf"
