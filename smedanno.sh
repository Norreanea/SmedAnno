#!/bin/bash

# =====================================================================
# SmedAnno: Genome annotation pipeline (Hybrid Docker + Conda Version)
# =====================================================================
# Master script 
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then

	set -e  # Exit immediately if a command exits with a non-zero status
	set -o pipefail  # Pipeline returns the exit status of the last command to fail
	# set -x  # Enable debugging for troubleshooting

	# Define a trap for cleaning up temporary files on exit
	# The 'trap' command ensures that the cleanup_temp_files function is called when the script exits,
	# for any reason (success, error, or user interrupt).
	TEMP_FILES=()
	cleanup_temp_files() {
		echo_blue "Cleaning up temporary files..."
		for file in "${TEMP_FILES[@]}"; do
			rm -f -- "$file"
		done
	}
	trap cleanup_temp_files EXIT

	# ---------------------------
	# Create Linux-native temporary directory for STAR
	# ---------------------------
	mkdir -p /tmp/temp_star
	chmod 777 /tmp/temp_star

	# ---------------------------
	# Define color variables for echo messages
	# ---------------------------
	RED='\033[0;31m'
	GREEN='\033[0;32m'
	BLUE='\033[1;34m'
	NC='\033[0m' # No Color

	# ---------------------------
	# Helper functions for colored echo
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
	# Function to display help
	# ---------------------------
	show_help() {
		echo "Usage: SmedAnno [OPTIONS]"
		echo ""
		echo "General mandatory options (if applicable based on steps selected):"
		echo "  --genomeRef PATH              Path to genome reference FASTA file"
		echo "  --dataDir PATH                Path to input RNA-Seq data directory (must contain short_reads and/or mix_reads folders)"
		echo ""
		echo "Mandatory options for specific steps:"
		echo "  --finalGTF PATH               Path to final GTF file (optional, required for Functional Annotation only)"
		echo "  --alignDir PATH               Path to directory containing BAM files (optional, required for Gene and Transcript Assembly only)"
		echo "  --SR_RB_gtf_dir PATH          Directory containing SR_RB.gtf files (optional, required for Merging Assemblies)"
		echo "  --SR_DN_gtf_dir PATH          Directory containing SR_DN.gtf files (optional, required for Merging Assemblies)"
		echo "  --MR_RB_gtf_dir PATH          Directory containing MR_RB.gtf files (optional, required for Merging Assemblies)"
		echo "  --MR_DN_gtf_dir PATH          Directory containing MR_DN.gtf files (optional, required for Merging Assemblies)"
		echo ""
		echo "Long-Read Trimming Options (Filtlong for Step 0):"
		echo "  --longread_min_len N          Minimum length for long reads (default: 1000)"
		echo "  --longread_keep_percent N     Keep the best N percent of reads (default: 95)"
		echo "  --longread_target_bases N     Target a total number of bases, keeping the best reads."
		echo "  NOTE: These options provide basic filtering. For best results, use dedicated tools before using SmedAnno."
		echo ""
		echo "Optional options:"
		echo "  --genomeDir PATH              Path to STAR genome directory (will be created if not provided)"
		echo "  --genomeGTF PATH              Path to genome annotation GTF file (optional, required for Reference-Based assembly)"
		echo "  --rrnaRef PATH                Path to rRNA reference FASTA file"
		echo "  --outputDir PATH              Path to output directory (default: ./outputDir)"
		echo "  --threads N                   Number of CPU threads to use (default: 8)"
		echo "  --Stringtie2minReadCoverage      Minimum read coverage allowed for the predicted transcripts (default: 1)"
		echo "	--Stringtie2minIsoformAbundance  Minimum isoform abundance of the predicted transcripts (default: 0.01)"
		echo "  --minOrfLength N              Minimum ORF length for TransDecoder (default: 100)"
		echo "  --maxExonLength N             Maximum allowed exon length (default: 10000)"           
		echo "  --maxTranscriptLength N       Maximum allowed transcript length (default: 100000)"      
		echo "  --steps LIST                  Comma-separated list of steps to run (0-10)"
		echo "  --all                         Run all steps sequentially"
		echo "  --functionalMethods METHODS   Comma-separated list of functional annotation methods to apply (BLASTp,BLASTx,PFAM,INTERPRO; default: BLASTp,BLASTx,PFAM,INTERPRO)"
		echo "  --stringtie VERSION               Set stringtie version for both short and mix reads (if not using individual overrides)"
		echo "  --stringtie_short VERSION         Set stringtie version for short reads (default: 2.1.1)"
		echo "  --stringtie_mix VERSION           Set stringtie version for mixed reads (default: 2.2.1)"
		echo "  --trimQual N                 TrimGalore quality cutoff (Phred, default 20)"
		echo "  --trimAdapter SEQ            Adapter sequence to remove               (default: auto-detection)"
		echo "  --trimGzip                   Gzip-compress trimmed FASTQ              (TRUE/FALSE)"
		echo "  --trimLen N                  Minimum read length after trimming       (default: 20 bp.)"
		echo "  --genomeType <type>           Specify the type of genome being processed. Options: 'nuclear', 'mito', 'mixed'."
		echo "                                If not set, the script will auto-detect 'mixed' if headers match --mitoPattern."
		echo "                                default: 'nuclear'."
		echo "  --mitoPattern <regex>         Regular expression to identify mitochondrial headers for 'mixed' mode."
		echo "                                default: '[Mm]ito|[Mm]itochondria|mtDNA'."
		echo "  --geneticCodeNucl <code>      Genetic code for nuclear transcripts. Default: 'Universal'."
		echo "  --geneticCodeMito <code>      Genetic code for mitochondrial transcripts. Default: 'Mitochondrial-Vertebrates'."
		echo ""
		echo_blue "Available genetic codes:"
		echo "  Acetabularia, Candida, Ciliate, Dasycladacean, Euplotid, Hexamita, Mesodinium,"
		echo "  Mitochondrial-Ascidian, Mitochondrial-Chlorophycean, Mitochondrial-Echinoderm, "
		echo "  Mitochondrial-Flatworm, Mitochondrial-Invertebrates,  Mitochondrial-Protozoan,"
		echo "  Mitochondrial-Pterobranchia, Mitochondrial-Scenedesmus_obliquus, Mitochondrial-Thraustochytrium,"
		echo "  Mitochondrial-Trematode, Mitochondrial-Vertebrates, Mitochondrial-Yeast,"
		echo "  Pachysolen_tannophilus, Peritrich, SR1_Gracilibacteria, Tetrahymena, Universal"
		echo ""
		echo "Filtering and threshold options (for Step 10):"
		echo "  --maxDistance N               Max distance between genes to be considered fragmented (default: 1000)"
		echo "  --largeIntronThreshold N      Min size of an intron to be considered unusually large (default: 100000)"
		echo "  --blastEvalue FLOAT           E-value threshold for BLAST searches (default: 1e-5)"
		echo "  --blastIdentity N             Percent identity threshold for BLAST searches (default: 25)"
		echo "  --pfamCEvalue FLOAT           Conditional E-value threshold for Pfam (default: 1e-5)"
		echo "  --pfamBitScore N              Bit score threshold for Pfam (default: 10)"
		echo "  --pfamCoverage N              Coverage threshold for Pfam domains (default: 50)"
		echo "  --interproEvalue FLOAT        E-value threshold for InterProScan hits (default: 1e-5)"
		echo ""
		echo "Steps:"
		echo "  0 - Read Trimming / Quality-Filtering"
		echo "  1 - rRNA Removal"
		echo "  2 - Read Alignment to Reference Genome"
		echo "  3 - Gene and Transcript Assembly"
		echo "  4 - Merge Reference-Based Assemblies"
		echo "  5 - Merge De Novo Assemblies and Create Pre-Final Annotation"
		echo "  6 - Filter Transcripts with Excessively Long Exons or Genomic Spans"
		echo "  7 - Isoform Comparison and Annotation"
		echo "  8 - GTF File Correction and Enhancement"
		echo "  9 - Functional Annotation and Filtering"
		echo "  10 - Integrate Functional Annotation (including Overlapped Genes and Transcripts, Reversed Duplicates, Fragmmented and Chimeric Genes Identification )"
		echo ""
		echo "Examples:"
		echo "  Run all steps:"
		echo "    ./run_smedanno.sh --genomeRef /path/to/genome.fa --dataDir ./data --outputDir ./output --threads 4 --all"
		echo ""
		echo "  Run Steps 8, 9, 10 only starting from a GTF file:"
		echo "    ./run_smedanno.sh --finalGTF /path/to/your.gtf --genomeRef /path/to/genome.fa --outputDir ./output --steps 8,9,10"
		exit 1
	}

	# ---------------------------
	# Function to determine if a file is gzipped
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
	# Helper function to find read files (handles .fq, .fq.gz, .fastq, .fastq.gz)
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
	# Helper function
	# Checks if a value exists in an array
	# Usage: contains "value" "${array[@]}"
	# ---------------------------
	contains() {
		local seeking="$1"; shift
		local in=("$@")
		for element in "${in[@]}"; do
			if [[ "$element" == "$seeking" ]]; then
				return 0 # Found
			fi
		done
		return 1 # Not found
	}

	# ---------------------------
	# Function to check global dependencies
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
	# Function to create conda environments and install tools
	# ---------------------------
	setup_conda_environments() {
		echo_blue "Verifying and creating Conda environments as needed..."

		# Determine environment names based on specified versions
		if [ "$STRINGTIE_SHORT_VERSION" = "$STRINGTIE_MIX_VERSION" ]; then
			ENV_COMMON="stringtie_${STRINGTIE_SHORT_VERSION//./}"
			ENV_SHORT="${ENV_COMMON}"
			ENV_MIX="${ENV_COMMON}"
		else
			ENV_SHORT="stringtie_short_${STRINGTIE_SHORT_VERSION//./}"
			ENV_MIX="stringtie_mix_${STRINGTIE_MIX_VERSION//./}"
		fi

		# Function to create an environment if it doesn't exist
		create_env_if_needed() {
			local env_name="$1"
			local st_version="$2"
			if ! conda env list | grep -q "^${env_name}\s"; then
				echo_green "Creating Conda environment '${env_name}' for StringTie v${st_version}..."
				if ! mamba create -y -n "${env_name}" -c bioconda "stringtie=${st_version}"; then
					echo_red "Failed to create conda environment '${env_name}'. Please check the StringTie version and your network connection."
					exit 1
				fi
			else
				echo_green "Conda environment '${env_name}' already exists."
			fi
		}

		# Create the environments
		create_env_if_needed "${ENV_SHORT}" "${STRINGTIE_SHORT_VERSION}"
		if [ "${ENV_SHORT}" != "${ENV_MIX}" ]; then
			create_env_if_needed "${ENV_MIX}" "${STRINGTIE_MIX_VERSION}"
		fi

		echo_green "Conda environments are set up correctly."
	}

	# ---------------------------
	# Function to create only necessary directories
	# ---------------------------
	setup_directories() {
		echo_blue "Setting up required directories..."

		# Define all possible directory variables
		GENOME_DIR="${OUTPUT_DIR}/genome"
		PREPROC_DIR="${OUTPUT_DIR}/preprocessing"
		ALIGN_DIR_INTERNAL="${OUTPUT_DIR}/alignment" # Use internal variable to avoid conflict with user --alignDir
		ASSEMBLY_DIR="${OUTPUT_DIR}/assembly"
		SR_RB_GTF_DIR_INTERNAL="${ASSEMBLY_DIR}/SR_RB"
		SR_DN_GTF_DIR_INTERNAL="${ASSEMBLY_DIR}/SR_DN"
		MR_RB_GTF_DIR_INTERNAL="${ASSEMBLY_DIR}/MR_RB"
		MR_DN_GTF_DIR_INTERNAL="${ASSEMBLY_DIR}/MR_DN"
		MERGE_DIR="${OUTPUT_DIR}/merging"
		ANNOTATION_DIR="${OUTPUT_DIR}/annotation"
		FUNCTIONAL_DIR="${OUTPUT_DIR}/functional_annotation"
		LOGS_DIR="${OUTPUT_DIR}/logs"

		# Always create the base output and logs directory
		mkdir -p "${OUTPUT_DIR}" "${LOGS_DIR}"

		# Map steps to the directories they require
		declare -A DIRS_FOR_STEP
		DIRS_FOR_STEP[0]="${PREPROC_DIR}"
		DIRS_FOR_STEP[1]="${PREPROC_DIR} ${ALIGN_DIR_INTERNAL}"
		DIRS_FOR_STEP[2]="${GENOME_DIR} ${ALIGN_DIR_INTERNAL}"
		DIRS_FOR_STEP[3]="${ASSEMBLY_DIR} ${SR_RB_GTF_DIR_INTERNAL} ${SR_DN_GTF_DIR_INTERNAL} ${MR_RB_GTF_DIR_INTERNAL} ${MR_DN_GTF_DIR_INTERNAL}"
		DIRS_FOR_STEP[4]="${MERGE_DIR}"
		DIRS_FOR_STEP[5]="${MERGE_DIR}"
		DIRS_FOR_STEP[6]="${MERGE_DIR}"
		DIRS_FOR_STEP[7]="${ANNOTATION_DIR}"
		DIRS_FOR_STEP[8]="${ANNOTATION_DIR}"
		DIRS_FOR_STEP[9]="${FUNCTIONAL_DIR} ${ANNOTATION_DIR}"
		DIRS_FOR_STEP[10]="${FUNCTIONAL_DIR} ${ANNOTATION_DIR}"

		# Build a unique list of directories to create based on selected steps
		DIRS_TO_CREATE=()
		for STEP in "${STEPS_TO_RUN[@]}"; do
			DIRS_TO_CREATE+=(${DIRS_FOR_STEP[$STEP]})
		done
		
		# Create only the unique, required directories
		for DIR_PATH in $(echo "${DIRS_TO_CREATE[@]}" | tr ' ' '\n' | sort -u); do
			echo_green "Creating directory: \"${DIR_PATH}\""
			mkdir -p "${DIR_PATH}"
		done
	}

	# ---------------------------
	# Function to run pre-installed InterProScan
	# ---------------------------
	run_interproscan() {
		local protein_file="$1"
		local output_basename="$2"
		local threads="$3"
		local ipr_executable="/opt/interproscan/interproscan.sh" 

		if [ ! -f "${ipr_executable}" ]; then
			echo_red "Error: InterProScan executable not found at \"${ipr_executable}\"."
			exit 1
		fi

		echo_green "Running pre-installed InterProScan..."
		local temp_dir
		temp_dir=$(mktemp -d -p "${FUNCTIONAL_DIR}" interproscan_temp_XXXXXX)
		TEMP_FILES+=("$temp_dir") # Add to trap for cleanup

		"${ipr_executable}" -i "${protein_file}" -f tsv \
							--cpu "${threads}" -o "${output_basename}" --goterms --pathways --disable-precalc \
							-T "${temp_dir}" \
							-appl AntiFam,CDD,Coils,FunFam,Gene3D,Hamap,MobiDBLite,PANTHER,Pfam,PIRSF,PIRSR,PRINTS,ProSitePatterns,ProSiteProfiles,SFLD,SMART,TIGRFAM

		echo_green "InterProScan finished. Temp directory will be cleaned up on exit."
	}
	# ---------------------------
	# Function to validate user-provided parameters
	# ---------------------------
	validate_parameters() {
		echo_blue "Validating user-provided parameters..."
		
        # Helper functions
		is_positive_integer() {
			[[ "$1" =~ ^[1-9][0-9]*$ ]]
		}
		is_positive_float() {
			[[ "$1" =~ ^[0-9]+([.][0-9]+)?(e-?[0-9]+)?$ ]]
		}
		is_valid_version() {
			local version_to_check="$1"; shift
			contains "${version_to_check}" "$@"
		}
		is_version_ge() {
			[[ "$(printf '%s\n' "$1" "$2" | sort -V | head -n1)" == "$2" ]]
		}
		
		# StringTie version validation
		local valid_stringtie_versions=("2.0.1" "2.0.2" "2.0.3" "2.0.4" "2.0.5" "2.0.6" "2.1.0" "2.1.1" "2.1.2" "2.1.3" "2.1.4" "2.1.6" "2.1.7" "2.2.0" "2.2.1" "2.2.2" "2.2.3" "3.0.0" "3.0.1")
		local min_mix_version="2.1.6"
		
		# StringTie version validation
		if ! is_valid_version "${STRINGTIE_SHORT_VERSION}" "${valid_stringtie_versions[@]}"; then
			echo_red "Error: Invalid version for --stringtie_short. You provided '${STRINGTIE_SHORT_VERSION}'."
			echo_red "Please use one of the following available versions: ${valid_stringtie_versions[*]}"
			show_help
		fi

		if ! is_valid_version "${STRINGTIE_MIX_VERSION}" "${valid_stringtie_versions[@]}"; then
			echo_red "Error: Invalid version for --stringtie_mix. You provided '${STRINGTIE_MIX_VERSION}'."
			echo_red "Please use one of the following available versions: ${valid_stringtie_versions[*]}"
			show_help
		fi

		if ! is_version_ge "${STRINGTIE_MIX_VERSION}" "${min_mix_version}"; then
			echo_red "Error: Invalid version for --stringtie_mix. The '--mix' functionality requires StringTie v${min_mix_version} or newer, but you provided v${STRINGTIE_MIX_VERSION}."
			show_help
		fi

		# Conditionally validate --genomeType if the user provided it
		if [ -n "${GENOME_TYPE_USER}" ]; then
			local valid_genome_types=("nuclear" "mito" "mixed")
			if ! contains "${GENOME_TYPE_USER}" "${valid_genome_types[@]}"; then
				echo_red "Error: Invalid --genomeType '${GENOME_TYPE_USER}'. Must be one of: ${valid_genome_types[*]}"
				show_help
			fi
		fi
		# Genetic code validation
		local valid_genetic_codes_nucl=(
			"Acetabularia" "Candida" "Ciliate" "Dasycladacean" "Euplotid" 
			"Hexamita" "Mesodinium" "Pachysolen_tannophilus" "Peritrich" 
			"SR1_Gracilibacteria" "Tetrahymena" "Universal"
		)
		local valid_genetic_codes_mito=(
			"Mitochondrial-Ascidian" "Mitochondrial-Chlorophycean" "Mitochondrial-Echinoderm" 
			"Mitochondrial-Flatworm" "Mitochondrial-Invertebrates" "Mitochondrial-Protozoan" 
			"Mitochondrial-Pterobranchia" "Mitochondrial-Scenedesmus_obliquus" "Mitochondrial-Thraustochytrium" 
			"Mitochondrial-Trematode" "Mitochondrial-Vertebrates" "Mitochondrial-Yeast"
		)
		if ! contains "${GENETIC_CODE_NUCL}" "${valid_genetic_codes_nucl[@]}"; then
			echo_red "Error: Invalid --geneticCodeNucl '${GENETIC_CODE_NUCL}'."
			echo_red "Please use one of the valid nuclear codes. See --help."
			show_help
		fi
		if ! contains "${GENETIC_CODE_MITO}" "${valid_genetic_codes_mito[@]}"; then
			echo_red "Error: Invalid --geneticCodeMito '${GENETIC_CODE_MITO}'."
			echo_red "Please use one of the valid mitochondrial codes. See --help."
			show_help
		fi

		is_positive_integer "${THREADS}" || { echo_red "Error: --threads must be a positive integer."; show_help; }
		is_positive_integer "${MIN_ORF_LENGTH}" || { echo_red "Error: --minOrfLength must be a positive integer."; show_help; }
		is_positive_integer "${MAX_EXON_LENGTH}" || { echo_red "Error: --maxExonLength must be a positive integer."; show_help; }
		is_positive_integer "${MAX_TRANSCRIPT_LENGTH}" || { echo_red "Error: --maxTranscriptLength must be a positive integer."; show_help; }
		is_positive_integer "${MAX_DISTANCE}" || { echo_red "Error: --maxDistance must be a positive integer."; show_help; }
		is_positive_integer "${LARGE_INTRON_THRESHOLD}" || { echo_red "Error: --largeIntronThreshold must be a positive integer."; show_help; }
		[[ "${LARGE_INTRON_THRESHOLD}" -gt 10000000 ]] && { echo_red "Warning: --largeIntronThreshold seems unusually high (>10,000,000). Please verify."; }
		is_positive_integer "${PFAM_BITSCORE}" || { echo_red "Error: --pfamBitScore must be a positive integer."; show_help; }
		is_positive_integer "${LONGREAD_MIN_LEN}" || { echo_red "Error: --longread_min_len must be a positive integer."; show_help; }
		is_positive_integer "${LONGREAD_KEEP_PERCENT}" || { echo_red "Error: --longread_keep_percent must be a positive integer."; show_help; }
		is_positive_integer "${LONGREAD_TARGET_BASES}" || { echo_red "Error: --longread_target_bases must be a positive integer."; show_help; }
		
		is_positive_integer "${TRIM_QUAL}" || { echo_red "Error: --trimQual must be a positive integer."; show_help; }
		[[ "${TRIM_QUAL}" -gt 40 ]] && { echo_red "Error: --trimQual should be within the Phred33 range (0-40)."; show_help; }
		is_positive_integer "${TRIM_LEN}" || { echo_red "Error: --trimLen must be a positive integer."; show_help; }

		is_positive_float "${BLAST_EVALUE}" || { echo_red "Error: --blastEvalue must be a positive number (e.g., 0.05 or 1e-5)."; show_help; }
		is_positive_float "${PFAM_CEVALUE}" || { echo_red "Error: --pfamCEvalue must be a positive number."; show_help; }
		is_positive_float "${INTERPRO_EVALUE}" || { echo_red "Error: --interproEvalue must be a positive number."; show_help; }
		is_positive_float "${STRINGTIE2_COVERAGE}" || { echo_red "Error: --Stringtie2minReadCoverage must be a positive number."; show_help; }
		is_positive_float "${STRINGTIE2_ABUNDANCE}" || { echo_red "Error: --Stringtie2minIsoformAbundance must be a positive number."; show_help; }

		is_positive_integer "${BLAST_IDENTITY}" || { echo_red "Error: --blastIdentity must be a positive integer."; show_help; }
		[[ "${BLAST_IDENTITY}" -lt 0 || "${BLAST_IDENTITY}" -gt 100 ]] && { echo_red "Error: --blastIdentity must be between 0 and 100."; show_help; }
		is_positive_integer "${PFAM_COVERAGE}" || { echo_red "Error: --pfamCoverage must be a positive integer."; show_help; }
		[[ "${PFAM_COVERAGE}" -lt 0 || "${PFAM_COVERAGE}" -gt 100 ]] && { echo_red "Error: --pfamCoverage must be between 0 and 100."; show_help; }

		echo_green "All parameters validated successfully."
	}
	
	# ---------------------------
	# Step 0: Read trimming / quality-filtering
	# ---------------------------
	step0_preprocess_trimming() {
		echo_blue  "Starting Step 0: Read trimming / filtering"
		
		# --- Short Read Trimming ---
		local TG_OPTS="--paired --cores ${THREADS} --quality ${TRIM_QUAL}"
		[[ -n "${TRIM_ADAPTER}" ]] && TG_OPTS+=" --adapter \"${TRIM_ADAPTER}\""
		[[ -n "${TRIM_LEN}" ]]     && TG_OPTS+=" --length ${TRIM_LEN}"
		[[ -n "${TRIM_GZIP}" ]]    && TG_OPTS+=" --gzip"

		if [[ ${#SHORT_ONLY_SAMPLES[@]} -gt 0 || ${#MIX_SAMPLES[@]} -gt 0 ]]; then
			echo_green "Processing short reads with Trim Galore..."
			for SAMPLE in "${SHORT_ONLY_SAMPLES[@]}" "${MIX_SAMPLES[@]}"; do
				local READ_DIR
				[[ -d "${SHORT_ONLY_DIR}" ]] && contains "$SAMPLE" "${SHORT_ONLY_SAMPLES[@]}" && READ_DIR="$SHORT_ONLY_DIR"
				[[ -d "${MIX_SHORT_DIR}" ]] && contains "$SAMPLE" "${MIX_SAMPLES[@]}" && READ_DIR="$MIX_SHORT_DIR"

				R1=$(find_read_file "${READ_DIR}" "${SAMPLE}" "1")
				R2=$(find_read_file "${READ_DIR}" "${SAMPLE}" "2")
				[[ -z "$R1" || -z "$R2" ]] && { echo_red "Missing FASTQ for \"${SAMPLE}\""; continue; }
				
				trim_galore ${TG_OPTS} -o "${PREPROC_DIR}" "${R1}" "${R2}"
			done
		fi

		# --- Long Read Filtering (for mixed samples) ---
		if [[ ${#MIX_SAMPLES[@]} -gt 0 ]]; then
			echo_green "Processing long reads with Filtlong..."
			for SAMPLE in "${MIX_SAMPLES[@]}"; do
				local LONG_READ_IN; LONG_READ_IN=$(find_read_file "${MIX_LONG_DIR}" "${SAMPLE}" "_long_reads")
				[ -z "${LONG_READ_IN}" ] && { echo_green "No long reads found for mixed sample \"${SAMPLE}\""; continue; }

				local LONG_READ_OUT="${PREPROC_DIR}/${SAMPLE}_long_reads_filtered.fastq"
				local FL_OPTS="--min_length ${LONGREAD_MIN_LEN} --keep_percent ${LONGREAD_KEEP_PERCENT}"
				[[ -n "${LONGREAD_TARGET_BASES}" ]] && FL_OPTS+=" --target_bases ${LONGREAD_TARGET_BASES}"
				
				# Use presets based on simple filename check
				if [[ "${LONG_READ_IN}" == *"[Pp]ac[Bb]io"* || "${LONG_READ_IN}" == *"[Hh][Ii][Ff][Ii]"* ]]; then
					FL_OPTS+=" --pacbio_corr"
				else # Default to ONT settings
					FL_OPTS+=" --ont_guppy"
				fi
				
				echo_green "Filtering long reads for sample \"${SAMPLE}\"..."
				filtlong ${FL_OPTS} "${LONG_READ_IN}" > "${LONG_READ_OUT}"
			done
		fi
		echo_blue "Step 0 completed: trimmed reads are in \"${OUTPUT_PATH}/preprocessing\""
	}

	# ---------------------------
	# Step 1: Preprocessing - rRNA Removal
	# ---------------------------
	step1_rrna_removal() {
		echo_blue  "Starting Step 1: rRNA Removal"

		if [ -z "${RRNA_REF}" ]; then
			echo_green  "No rRNA reference provided. Skipping rRNA removal."
			return
		fi

		echo_green  "Calculating rRNA reference length for genomeSAindexNbases parameter"
		
		RRNA_LENGTH=$(gawk '/^>/ {next} {len+=length($0)} END {print len}' "${RRNA_REF}")
		RRNA_SA_INDEX_NBASES=$(python3 -c "import math; l=math.log(${RRNA_LENGTH},2); n=int(l/2 -1); print(min(14, n))")

		echo_green  "rRNA reference length: ${RRNA_LENGTH}"
		echo_green  "Calculated genomeSAindexNbases for rRNA reference: ${RRNA_SA_INDEX_NBASES}"

		local RRNA_REF_INDEX="${OUTPUT_DIR}/rRNA_STAR_index"
		if [ ! -d "${RRNA_REF_INDEX}" ] || [ ! -f "${RRNA_REF_INDEX}/SA" ]; then
			echo_green  "Building STAR index for rRNA reference at \"${RRNA_REF_INDEX}\""
			mkdir -p "${RRNA_REF_INDEX}"
			STAR --runThreadN "${THREADS}" \
				 --runMode genomeGenerate \
				 --genomeDir "${RRNA_REF_INDEX}" \
				 --genomeFastaFiles "${RRNA_REF}" \
				 --genomeSAindexNbases "${RRNA_SA_INDEX_NBASES}" \
				 --outTmpDir /tmp/temp_star/temp
		else
			echo_green  "rRNA STAR index already exists at \"${RRNA_REF_INDEX}\", skipping index generation."
		fi

		if compgen -G "./Log.out" > /dev/null; then
			mv ./Log.out "${LOGS_DIR}/rRNA_genomeGenerate_Log.out"
		fi

		# Process short-only samples
		if [ "${#SHORT_ONLY_SAMPLES[@]}" -gt 0 ]; then
			for SAMPLE in "${SHORT_ONLY_SAMPLES[@]}"; do
				echo_green  "Processing Short-Only Sample: \"${SAMPLE}\""
				local READ1; READ1=$(find_read_file "${SHORT_ONLY_DIR}" "${SAMPLE}" "1")
				local READ2; READ2=$(find_read_file "${SHORT_ONLY_DIR}" "${SAMPLE}" "2")
				[ -z "${READ1}" ] || [ -z "${READ2}" ] && { echo_red "Error: One or both short read files for sample \"${SAMPLE}\" not found."; exit 1; }
				local READ_FILES_CMD; READ_FILES_CMD=$([ "$(is_gzipped "${READ1}")" == "yes" ] && echo "zcat" || echo "cat")

				local NON_RRNA_READ1="${PREPROC_DIR}/${SAMPLE}_non_rrna_R1.fq.gz"
				local NON_RRNA_READ2="${PREPROC_DIR}/${SAMPLE}_non_rrna_R2.fq.gz"

				echo_green  "Removing rRNA from short reads using STAR for Sample: \"${SAMPLE}\""
				STAR --runThreadN "${THREADS}" \
					 --genomeDir "${RRNA_REF_INDEX}" \
					 --readFilesIn "${READ1}" "${READ2}" \
					 --readFilesCommand "${READ_FILES_CMD}" \
					 --outFilterType BySJout --outFilterMultimapNmax 1 --outFilterMismatchNmax 0 \
					 --outReadsUnmapped Fastx --outSAMtype BAM SortedByCoordinate \
					 --outFileNamePrefix "${PREPROC_DIR}/${SAMPLE}_STAR_rrna_" \
					 --outTmpDir /tmp/temp_star/temp

				gzip "${PREPROC_DIR}/${SAMPLE}_STAR_rrna_Unmapped.out.mate1"
				gzip "${PREPROC_DIR}/${SAMPLE}_STAR_rrna_Unmapped.out.mate2"
				mv "${PREPROC_DIR}/${SAMPLE}_STAR_rrna_Unmapped.out.mate1.gz" "${NON_RRNA_READ1}"
				mv "${PREPROC_DIR}/${SAMPLE}_STAR_rrna_Unmapped.out.mate2.gz" "${NON_RRNA_READ2}"
			done
		fi

		# Process mixed samples
		if [ "${#MIX_SAMPLES[@]}" -gt 0 ]; then
			for SAMPLE in "${MIX_SAMPLES[@]}"; do
				echo_green "Processing Mixed Sample: \"${SAMPLE}\""
				local MIX_READ1; MIX_READ1=$(find_read_file "${MIX_SHORT_DIR}" "${SAMPLE}" "1")
				local MIX_READ2; MIX_READ2=$(find_read_file "${MIX_SHORT_DIR}" "${SAMPLE}" "2")
				[ -z "${MIX_READ1}" ] || [ -z "${MIX_READ2}" ] && { echo_red "Error: Short read files for mixed sample \"${SAMPLE}\" not found."; exit 1; }
				local MIX_READ_FILES_CMD; MIX_READ_FILES_CMD=$([ "$(is_gzipped "${MIX_READ1}")" == "yes" ] && echo "zcat" || echo "cat")

				echo_green  "Running STAR aligner for Mixed Sample: \"${SAMPLE}\" (Short Reads)"
				STAR --runThreadN "${THREADS}" \
					 --genomeDir "${RRNA_REF_INDEX}" \
					 --readFilesIn "${MIX_READ1}" "${MIX_READ2}" \
					 --readFilesCommand "${MIX_READ_FILES_CMD}" \
					 --outFilterType BySJout --outFilterMultimapNmax 1 --outFilterMismatchNmax 0 \
					 --outReadsUnmapped Fastx --outSAMtype BAM SortedByCoordinate \
					 --outFileNamePrefix "${PREPROC_DIR}/${SAMPLE}_STAR_rrna_" \
					 --outTmpDir /tmp/temp_star/temp

				gzip "${PREPROC_DIR}/${SAMPLE}_STAR_rrna_Unmapped.out.mate1"
				gzip "${PREPROC_DIR}/${SAMPLE}_STAR_rrna_Unmapped.out.mate2"
				mv "${PREPROC_DIR}/${SAMPLE}_STAR_rrna_Unmapped.out.mate1.gz" "${PREPROC_DIR}/${SAMPLE}_non_rrna_R1.fq.gz"
				mv "${PREPROC_DIR}/${SAMPLE}_STAR_rrna_Unmapped.out.mate2.gz" "${PREPROC_DIR}/${SAMPLE}_non_rrna_R2.fq.gz"

				local MIX_LONG_READ; MIX_LONG_READ=$(find_read_file "${MIX_LONG_DIR}" "${SAMPLE}" "_long_reads")
				[ -z "${MIX_LONG_READ}" ] && { echo_green "No long reads found for Mixed Sample \"${SAMPLE}\""; continue; }

				echo_green  "Removing rRNA from mixed sample's long reads using minimap2 for Sample: \"${SAMPLE}\""
				local ALIGNED_BAM_TEMP; ALIGNED_BAM_TEMP=$(mktemp -p "${ALIGN_DIR_INTERNAL}" -u --suffix .bam)
				TEMP_FILES+=("$ALIGNED_BAM_TEMP")
				
				local MINIMAP_CMD="minimap2 -t ${THREADS} -ax map-pb \"${RRNA_REF}\""
				local LONG_READ_INPUT=$([ "$(is_gzipped "${MIX_LONG_READ}")" == "yes" ] && echo "<(zcat \"${MIX_LONG_READ}\")" || echo "\"${MIX_LONG_READ}\"")
				eval "${MINIMAP_CMD} ${LONG_READ_INPUT}" | samtools view -Sb - > "${ALIGNED_BAM_TEMP}"

				if [ -f "${ALIGNED_BAM_TEMP}" ]; then
					echo_green  "Extracting unmapped (non-rRNA) long reads for Sample: \"${SAMPLE}\""
					samtools view -b -f 4 "${ALIGNED_BAM_TEMP}" | samtools fastq - > "${PREPROC_DIR}/${SAMPLE}_non_rrna_long_reads.fq.gz"
				fi
			done
		fi
	}

	# ---------------------------
	# Step 2: Read Alignment to Reference Genome
	# ---------------------------
	step2_read_alignment() {
		echo_blue  "Starting Step 2: Read Alignment to the Reference Genome"

		if [ ! -d "${GENOME_DIR}" ] || [ ! -f "${GENOME_DIR}/SA" ]; then
			echo_green  "Genome index not found. Generating genome index as part of Step 2."
			local GENOME_LENGTH; GENOME_LENGTH=$(gawk '/^>/ {next} {len+=length($0)} END {print len}' "${GENOME_REF}")
			local SA_INDEX_NBASES; SA_INDEX_NBASES=$(python3 -c "import math; l=math.log(${GENOME_LENGTH},2); n=int(l/2 -1); print(min(14, n))")

			echo_green  "Genome length: ${GENOME_LENGTH}, Calculated genomeSAindexNbases: ${SA_INDEX_NBASES}"
			mkdir -p "${GENOME_DIR}"
			STAR --runThreadN "${THREADS}" \
				 --runMode genomeGenerate \
				 --genomeDir "${GENOME_DIR}" \
				 --genomeFastaFiles "${GENOME_REF}" \
				 --genomeSAindexNbases "${SA_INDEX_NBASES}" \
				 --outTmpDir /tmp/temp_star/temp
		else
			echo_green  "Genome index already exists at \"${GENOME_DIR}\", skipping index generation."
		fi
		
		if compgen -G "./Log.out" > /dev/null; then
			mv ./Log.out "${LOGS_DIR}/ref_genomeGenerate_Log.out"
		fi

		# Process short-only samples
		if [ "${#SHORT_ONLY_SAMPLES[@]}" -gt 0 ]; then
			for SAMPLE in "${SHORT_ONLY_SAMPLES[@]}"; do
				echo_green  "Aligning Short-Only Sample: \"${SAMPLE}\""
				local READ1; READ1=$([ -n "${RRNA_REF}" ] && echo "${PREPROC_DIR}/${SAMPLE}_non_rrna_R1.fq.gz" || find_read_file "${SHORT_ONLY_DIR}" "${SAMPLE}" "1")
				local READ2; READ2=$([ -n "${RRNA_REF}" ] && echo "${PREPROC_DIR}/${SAMPLE}_non_rrna_R2.fq.gz" || find_read_file "${SHORT_ONLY_DIR}" "${SAMPLE}" "2")
				[ -z "${READ1}" ] || [ -z "${READ2}" ] && { echo_red  "Error: Read files for sample \"${SAMPLE}\" not found."; exit 1; }
				local READ_FILES_CMD; READ_FILES_CMD=$([ "$(is_gzipped "${READ1}")" == "yes" ] && echo "zcat" || echo "cat")

				echo_green  "Running STAR aligner for Short-Only Sample: \"${SAMPLE}\""
				STAR --runThreadN "${THREADS}" \
					 --genomeDir "${GENOME_DIR}" \
					 --readFilesIn "${READ1}" "${READ2}" \
					 --readFilesCommand "${READ_FILES_CMD}" \
					 --outSAMtype BAM SortedByCoordinate \
					 --outFileNamePrefix "${ALIGN_DIR_INTERNAL}/${SAMPLE}_STAR_" \
					 --outTmpDir /tmp/temp_star/temp

				local ALIGNED_SORTED_BAM="${ALIGN_DIR_INTERNAL}/${SAMPLE}_STAR_Aligned.sortedByCoord.out.bam"
				[ ! -f "${ALIGNED_SORTED_BAM}" ] && { echo_red "Error: STAR did not produce the expected sorted BAM file for sample \"${SAMPLE}\"."; exit 1; }

				echo_green  "Adding XS attribute to BAM file for Sample: \"${SAMPLE}\""
				samtools view -h "${ALIGNED_SORTED_BAM}" | \
					gawk -v strType=2 -f ./tagXSstrandedData.awk | \
					samtools view -bSq 10 -F 4 - > "${ALIGN_DIR_INTERNAL}/${SAMPLE}_STAR_Aligned.sortedByCoord.outXS.bam"
				rm "${ALIGNED_SORTED_BAM}"
			done
		fi

		# Process mixed samples
		if [ "${#MIX_SAMPLES[@]}" -gt 0 ]; then
			for SAMPLE in "${MIX_SAMPLES[@]}"; do
				echo_green  "Aligning Mixed Sample: \"${SAMPLE}\""
				local MIX_READ1; MIX_READ1=$([ -n "${RRNA_REF}" ] && echo "${PREPROC_DIR}/${SAMPLE}_non_rrna_R1.fq.gz" || find_read_file "${MIX_SHORT_DIR}" "${SAMPLE}" "1")
				local MIX_READ2; MIX_READ2=$([ -n "${RRNA_REF}" ] && echo "${PREPROC_DIR}/${SAMPLE}_non_rrna_R2.fq.gz" || find_read_file "${MIX_SHORT_DIR}" "${SAMPLE}" "2")
				[ -z "${MIX_READ1}" ] || [ -z "${MIX_READ2}" ] && { echo_red "Error: Short read files for mixed sample \"${SAMPLE}\" not found."; exit 1; }
				local MIX_READ_FILES_CMD; MIX_READ_FILES_CMD=$([ "$(is_gzipped "${MIX_READ1}")" == "yes" ] && echo "zcat" || echo "cat")

				echo_green  "Running STAR aligner for Mixed Sample: \"${SAMPLE}\" (Short Reads)"
				STAR --runThreadN "${THREADS}" \
					 --genomeDir "${GENOME_DIR}" \
					 --readFilesIn "${MIX_READ1}" "${MIX_READ2}" \
					 --readFilesCommand "${MIX_READ_FILES_CMD}" \
					 --outFilterType BySJout --outFilterMultimapNmax 1 --outFilterMismatchNmax 0 \
					 --outSAMtype BAM SortedByCoordinate \
					 --outFileNamePrefix "${ALIGN_DIR_INTERNAL}/${SAMPLE}_STAR_" \
					 --outTmpDir /tmp/temp_star/temp
				
				local ALIGNED_SHORT_SORTED_BAM="${ALIGN_DIR_INTERNAL}/${SAMPLE}_STAR_Aligned.sortedByCoord.out.bam"
				[ ! -f "${ALIGNED_SHORT_SORTED_BAM}" ] && { echo_red "Error: STAR did not produce the expected sorted BAM file for mixed sample \"${SAMPLE}\"."; exit 1; }

				echo_green  "Adding XS attribute to short BAM file for Sample: \"${SAMPLE}\""
				samtools view -h "${ALIGNED_SHORT_SORTED_BAM}" | gawk -v strType=2 -f ./tagXSstrandedData.awk | samtools view -bSq 10 -F 4 - > "${ALIGN_DIR_INTERNAL}/${SAMPLE}_STAR_Aligned.sortedByCoord.outXS.bam"
				rm "${ALIGNED_SHORT_SORTED_BAM}"
				
				local MIX_LONG_READ; MIX_LONG_READ=$([ -n "${RRNA_REF}" ] && echo "${PREPROC_DIR}/${SAMPLE}_non_rrna_long_reads.fq.gz" || find_read_file "${MIX_LONG_DIR}" "${SAMPLE}" "_long_reads")
				if [ -n "${MIX_LONG_READ}" ] && [ -f "${MIX_LONG_READ}" ]; then
					echo_green  "Running minimap2 aligner for Mixed Sample: \"${SAMPLE}\" (Long Reads)"
					local ALIGNED_LONG_BAM="${ALIGN_DIR_INTERNAL}/${SAMPLE}_long_aligned.bam"
					minimap2 -t "${THREADS}" -ax splice "${GENOME_REF}" "${MIX_LONG_READ}" | samtools sort -@ "${THREADS}" -o "${ALIGNED_LONG_BAM}" -
				else
					echo_green  "No long reads found or file does not exist for Mixed Sample: \"${SAMPLE}\""
				fi
			done
		fi
		echo_blue  "Step 2 Completed: Read Alignment"
	}

	# ---------------------------
	# Step 3: Gene and Transcript Assembly
	# ---------------------------
	step3_gene_transcript_assembly() {
		echo_blue  "Starting Step 3: Gene and Transcript Assembly"

		# Use the user-provided alignDir if available, otherwise use the internal one
		local current_align_dir=${ALIGN_DIR:-$ALIGN_DIR_INTERNAL}
		
		# Process short-only samples
		if [ "${#SHORT_ONLY_SAMPLES[@]}" -gt 0 ]; then
			for SAMPLE in "${SHORT_ONLY_SAMPLES[@]}"; do
				echo_green  "Assembling transcripts for Short-Only Sample: \"${SAMPLE}\""
				local ALIGNED_SORTED_BAM="${current_align_dir}/${SAMPLE}_STAR_Aligned.sortedByCoord.outXS.bam"
				
				if [ -n "${GENOME_GTF}" ]; then
					echo_green  "Running StringTie2 for Reference-Based Assembly (RB, SR)"
					conda run -n ${ENV_SHORT} stringtie -p "${THREADS}" -G "${GENOME_GTF}" -c "${STRINGTIE2_COVERAGE}" -f "${STRINGTIE2_ABUNDANCE}" -o "${SR_RB_GTF_DIR_INTERNAL}/${SAMPLE}_SR_RB.gtf" "${ALIGNED_SORTED_BAM}"
				fi

				echo_green  "Running StringTie2 for De Novo Assembly (DN, SR)"
				conda run -n ${ENV_SHORT} stringtie -p "${THREADS}" -c "${STRINGTIE2_COVERAGE}" -f "${STRINGTIE2_ABUNDANCE}" -o "${SR_DN_GTF_DIR_INTERNAL}/${SAMPLE}_SR_DN.gtf" "${ALIGNED_SORTED_BAM}"
			done
		fi

		# Process mixed samples
		if [ "${#MIX_SAMPLES[@]}" -gt 0 ]; then
			for SAMPLE in "${MIX_SAMPLES[@]}"; do
				echo_green  "Assembling transcripts for Mixed Sample: \"${SAMPLE}\""
				local ALIGNED_SHORT_SORTED_BAM="${current_align_dir}/${SAMPLE}_STAR_Aligned.sortedByCoord.outXS.bam"
				local ALIGNED_LONG_BAM="${current_align_dir}/${SAMPLE}_long_aligned.bam"
				[ ! -f "${ALIGNED_LONG_BAM}" ] && { echo_red "No aligned long reads BAM found for Sample: \"${SAMPLE}\", skipping mixed assembly."; continue; }

				if [ -n "${GENOME_GTF}" ]; then
					echo_green  "Running StringTie2 for Reference-Based Assembly (RB, MR)"
					conda run -n ${ENV_MIX} stringtie --mix -p "${THREADS}" -G "${GENOME_GTF}" -c "${STRINGTIE2_COVERAGE}" -f "${STRINGTIE2_ABUNDANCE}" -o "${MR_RB_GTF_DIR_INTERNAL}/${SAMPLE}_MR_RB.gtf" "${ALIGNED_SHORT_SORTED_BAM}" "${ALIGNED_LONG_BAM}"
				fi

				echo_green  "Running StringTie2 for De Novo Assembly (DN, MR)"
				conda run -n ${ENV_MIX} stringtie --mix -p "${THREADS}" -c "${STRINGTIE2_COVERAGE}" -f "${STRINGTIE2_ABUNDANCE}" -o "${MR_DN_GTF_DIR_INTERNAL}/${SAMPLE}_MR_DN.gtf" "${ALIGNED_SHORT_SORTED_BAM}" "${ALIGNED_LONG_BAM}"
			done
		fi
		echo_blue  "Step 3 Completed: Gene and Transcript Assembly"
	}

	# ---------------------------
	# Step 4: Merge Reference-Based Assemblies
	# ---------------------------
	step4_merging_transcripts_I() {
		echo_blue "Starting Step 4: Merge Reference-Based Assemblies"
		if [ -z "${GENOME_GTF}" ]; then echo_green "No reference GTF provided. Skipping Reference-Based merging."; return; fi
	
		local target_dir_sr=${SR_RB_GTF_DIR:-$SR_RB_GTF_DIR_INTERNAL}
		local target_dir_mr=${MR_RB_GTF_DIR:-$MR_RB_GTF_DIR_INTERNAL}
		
		# Merge SR
		if [ -d "${target_dir_sr}" ]; then
			local LIST_SR_RB_GTF="${MERGE_DIR}/SR_RB_gtf_list.txt"
			find "${target_dir_sr}" -name "*.gtf" -type f -size +1c > "${LIST_SR_RB_GTF}"
			if [ -s "${LIST_SR_RB_GTF}" ]; then
				echo_green "Merging Reference-Based Assemblies for Short Reads (SR)"
				conda run -n ${ENV_MIX} stringtie --merge -p "${THREADS}" -G "${GENOME_GTF}" -o "${MERGE_DIR}/merged_SR_RB.gtf" "${LIST_SR_RB_GTF}"
			fi
		fi

		# Merge MR
		if [ -d "${target_dir_mr}" ]; then
			local LIST_MR_RB_GTF="${MERGE_DIR}/MR_RB_gtf_list.txt"
			find "${target_dir_mr}" -name "*.gtf" -type f -size +1c > "${LIST_MR_RB_GTF}"
			if [ -s "${LIST_MR_RB_GTF}" ]; then
				echo_green "Merging Reference-Based Assemblies for Mixed Reads (MR)"
				conda run -n ${ENV_MIX} stringtie --merge -p "${THREADS}" -G "${GENOME_GTF}" -o "${MERGE_DIR}/merged_MR_RB.gtf" "${LIST_MR_RB_GTF}"
			fi
		fi
		echo_blue "Step 4 Completed: Merge Reference-Based Assemblies"
	}

	# ---------------------------
	# Step 5: Merge De Novo Assemblies and Create Pre-Final Annotation
	# ---------------------------
	step5_merging_transcripts_II() {
		echo_blue "Starting Step 5: Merge De Novo Assemblies and Create Pre-Final Annotation"

		local target_dir_sr_dn=${SR_DN_GTF_DIR:-$SR_DN_GTF_DIR_INTERNAL}
		local target_dir_mr_dn=${MR_DN_GTF_DIR:-$MR_DN_GTF_DIR_INTERNAL}
		
		# Merge SR De Novo
		if [ -d "${target_dir_sr_dn}" ]; then
			local LIST_SR_DN_GTF="${MERGE_DIR}/SR_DN_gtf_list.txt"
			find "${target_dir_sr_dn}" -name "*.gtf" -type f -size +1c > "${LIST_SR_DN_GTF}"
			if [ -s "${LIST_SR_DN_GTF}" ]; then
				echo_green "Merging De Novo Assemblies for Short Reads (SR)"
				conda run -n ${ENV_MIX} stringtie --merge -p "${THREADS}" -o "${MERGE_DIR}/merged_SR_DN.gtf" "${LIST_SR_DN_GTF}"
			fi
		fi
		
		# Merge MR De Novo
		if [ -d "${target_dir_mr_dn}" ]; then
			local LIST_MR_DN_GTF="${MERGE_DIR}/MR_DN_gtf_list.txt"
			find "${target_dir_mr_dn}" -name "*.gtf" -type f -size +1c > "${LIST_MR_DN_GTF}"
			if [ -s "${LIST_MR_DN_GTF}" ]; then
				echo_green "Merging De Novo Assemblies for Mixed Reads (MR)"
				conda run -n ${ENV_MIX} stringtie --merge -p "${THREADS}" -o "${MERGE_DIR}/merged_MR_DN.gtf" "${LIST_MR_DN_GTF}"
			fi
		fi

		echo_green "Creating pre-final annotation GTF"
		local merge_list=()
		[ -s "${MERGE_DIR}/merged_SR_RB.gtf" ] && merge_list+=("${MERGE_DIR}/merged_SR_RB.gtf")
		[ -s "${MERGE_DIR}/merged_MR_RB.gtf" ] && merge_list+=("${MERGE_DIR}/merged_MR_RB.gtf")
		[ -s "${MERGE_DIR}/merged_SR_DN.gtf" ] && merge_list+=("${MERGE_DIR}/merged_SR_DN.gtf")
		[ -s "${MERGE_DIR}/merged_MR_DN.gtf" ] && merge_list+=("${MERGE_DIR}/merged_MR_DN.gtf")
		
		if [ ${#merge_list[@]} -gt 0 ]; then
			echo_green "Merging all available assemblies into pre-final annotation"
			conda run -n ${ENV_MIX} stringtie --merge -p "${THREADS}" -o "${MERGE_DIR}/prefinal_annotation.gtf" "${merge_list[@]}"
		else
			echo_red "No merged GTF files available to create a pre-final annotation. Skipping."
		fi
		echo_blue "Step 5 Completed: Merge De Novo Assemblies and Create Pre-Final Annotation"
	}

	# ---------------------------
	# Step 6: Filter Out Transcripts with Excessively Long Exons or Genomic Spans
	# ---------------------------
	step6_filter_transcripts() {
		echo_blue "Starting Step 6: Filtering Transcripts"
		local PREFINAL_GTF="${MERGE_DIR}/prefinal_annotation.gtf"
		local FILTERED_GTF="${MERGE_DIR}/filtered_annotation.gtf"
		[ ! -f "${PREFINAL_GTF}" ] && { echo_red "Error: \"${PREFINAL_GTF}\" not found. Skipping filtering."; return; }

		echo_green "Identifying transcripts with exons > ${MAX_EXON_LENGTH} nt or spans > ${MAX_TRANSCRIPT_LENGTH} nt"
		
		local transcripts_to_exclude_file; transcripts_to_exclude_file=$(mktemp -p "${MERGE_DIR}" transcripts_to_exclude_XXXXXX.txt)
		TEMP_FILES+=("$transcripts_to_exclude_file")

		# Use gawk to find transcript IDs to exclude in a single pass
		gawk -v max_exon="${MAX_EXON_LENGTH}" -v max_span="${MAX_TRANSCRIPT_LENGTH}" '
			function get_id(line) {
				if (match(line, /transcript_id "([^"]+)"/, arr)) { return arr[1]; }
				return "";
			}
			$3 == "exon" {
				if (($5 - $4 + 1) > max_exon) { bad[get_id($0)] = 1; }
			}
			$3 == "transcript" {
				if (($5 - $4 + 1) > max_span) { bad[get_id($0)] = 1; }
			}
			END { for (t_id in bad) { print t_id; } }
		' "${PREFINAL_GTF}" | sort -u > "${transcripts_to_exclude_file}"

		if [ -s "${transcripts_to_exclude_file}" ]; then
			echo_red "Found $(wc -l < "${transcripts_to_exclude_file}") transcripts to exclude. Filtering GTF."
			gawk 'NR==FNR {exclude[$1]=1; next} !/^#/ { if (match($0, /transcript_id "([^"]+)"/, arr)) { if (!(arr[1] in exclude)) print; } else { print; } } /^#/ {print}' \
				"${transcripts_to_exclude_file}" "${PREFINAL_GTF}" > "${FILTERED_GTF}"
		else
			echo_green "No transcripts exceeded length thresholds. Copying to filtered location."
			cp "${PREFINAL_GTF}" "${FILTERED_GTF}"
		fi
		echo_blue "Step 6 Completed: Filtering Transcripts"
	}

	# ---------------------------
	# Step 7: Isoform Comparison and Annotation
	# ---------------------------
	step7_isoform_comparison() {
		echo_blue  "Starting Step 7: Isoform Comparison and Annotation"

		local GTF_TO_COMPARE; GTF_TO_COMPARE=$([ -n "${finalGTF}" ] && echo "${finalGTF}" || echo "${MERGE_DIR}/filtered_annotation.gtf")
		[ ! -f "${GTF_TO_COMPARE}" ] && { echo_red "GTF for comparison not found. Skipping Step 7."; return; }

		local COMP_OUTPUT_PREFIX="${ANNOTATION_DIR}/gffcomp"
		echo_green  "Running gffcompare on \"${GTF_TO_COMPARE}\""
		
		declare -a GFFCOMPARE_ARGS
		GFFCOMPARE_ARGS=("-o" "${COMP_OUTPUT_PREFIX}" "--chr-stats" "--strict-match" "${GTF_TO_COMPARE}")
		[ -n "${GENOME_GTF}" ] && GFFCOMPARE_ARGS=("-r" "${GENOME_GTF}" "-G" "${GFFCOMPARE_ARGS[@]}")
		
		gffcompare "${GFFCOMPARE_ARGS[@]}"
		echo_blue  "Step 7 Completed: Isoform Comparison and Annotation"
	}

	# ---------------------------
	# Step 8: GTF File Correction and Enhancement
	# ---------------------------
	step8_gtf_correction() {
		echo_blue  "Starting Step 8: GTF Correction and Enhancement"

		local GTF_TO_CORRECT; GTF_TO_CORRECT=$([ -f "${ANNOTATION_DIR}/gffcomp.combined.gtf" ] && echo "${ANNOTATION_DIR}/gffcomp.combined.gtf" || echo "${MERGE_DIR}/filtered_annotation.gtf")
		[ -n "${finalGTF}" ] && [ -f "${finalGTF}" ] && GTF_TO_CORRECT="${finalGTF}"
		[ ! -f "${GTF_TO_CORRECT}" ] && { echo_red "No suitable GTF found for correction. Skipping Step 8."; return; }

		local CORRECTED_GTF="${ANNOTATION_DIR}/corrected.gtf"
		local CORRECTED_INTRONS_GTF="${ANNOTATION_DIR}/corrected_with_introns.gtf"

		echo_green  "Running AGAT to correct and enhance GTF file: \"${GTF_TO_CORRECT}\""
		agat_convert_sp_gxf2gxf.pl -g "${GTF_TO_CORRECT}" -o "${CORRECTED_GTF}"
		agat_sp_add_introns.pl -g "${CORRECTED_GTF}" -o "${CORRECTED_INTRONS_GTF}"
		
		if compgen -G "./*.agat.log" > /dev/null; then
			mv ./*.agat.log "${LOGS_DIR}/"
		fi
		echo_blue  "Step 8 Completed: GTF Correction and Enhancement"
	}

	# ---------------------------
	# Step 9: Functional Annotation and Filtering
	# ---------------------------
	step9_functional_annotation() {
		echo_blue  "Starting Step 9: Functional Annotation and Filtering"
		
		local ANNOTATED_GTF=$([ -n "${finalGTF}" ] && echo "${finalGTF}" || echo "${ANNOTATION_DIR}/corrected_with_introns.gtf")
		[ ! -f "${ANNOTATED_GTF}" ] && { echo_red "No valid GTF file found for functional annotation. Skipping Step 9."; return; }
		echo_green "Using GTF file for annotation: \"${ANNOTATED_GTF}\""

		local SANITIZED_GTF="${ANNOTATION_DIR}/sanitized.final.gtf"
		echo_green "Sanitizing GTF file to correct non-standard whitespace..."
		sed 's/\xC2\xA0/\t/g' "${ANNOTATED_GTF}" > "${SANITIZED_GTF}"

		echo_green "Detecting genome composition from sanitized GTF..."
		local has_nucl=false; local has_mito=false
		[ -n "$(gawk -v pattern="${MITO_PATTERN}" '!/^#/ && $1 !~ pattern' "${SANITIZED_GTF}")" ] && has_nucl=true
		[ -n "$(gawk -v pattern="${MITO_PATTERN}" '!/^#/ && $1 ~ pattern' "${SANITIZED_GTF}")" ] && has_mito=true
		
		local GENOME_TYPE_DETECTED
		if [[ "$has_nucl" == true && "$has_mito" == true ]]; then GENOME_TYPE_DETECTED="mixed"
		elif [[ "$has_nucl" == true ]]; then GENOME_TYPE_DETECTED="nuclear"
		elif [[ "$has_mito" == true ]]; then GENOME_TYPE_DETECTED="mito"
		else GENOME_TYPE_DETECTED="nuclear"; echo_red "Warning: Could not classify contigs. Defaulting to 'nuclear'."
		fi
		# Override auto-detection if user specified the type
		[ -n "${GENOME_TYPE_USER}" ] && GENOME_TYPE="${GENOME_TYPE_USER}" || GENOME_TYPE="${GENOME_TYPE_DETECTED}"
		echo_blue "Genome type set to: \"${GENOME_TYPE}\""

		local NUCL_TRANSCRIPTS_FA="${FUNCTIONAL_DIR}/nuclear_transcripts.fa"
		local MITO_TRANSCRIPTS_FA="${FUNCTIONAL_DIR}/mito_transcripts.fa"
		local LONGEST_ORFS_PEP="${FUNCTIONAL_DIR}/longest_orfs.pep"

		if [[ "$GENOME_TYPE" == "mixed" ]]; then
			echo_green "Creating separate nuclear and mitochondrial transcript FASTA files..."
			local nucl_gtf_tmp; nucl_gtf_tmp=$(mktemp -p "${FUNCTIONAL_DIR}" nucl_XXXXXX.gtf)
			local mito_gtf_tmp; mito_gtf_tmp=$(mktemp -p "${FUNCTIONAL_DIR}" mito_XXXXXX.gtf)
			TEMP_FILES+=("$nucl_gtf_tmp" "$mito_gtf_tmp")

			gawk -v pattern="${MITO_PATTERN}" '!/^#/{ if ($1 !~ pattern) print > "'"${nucl_gtf_tmp}"'"; else print > "'"${mito_gtf_tmp}"'"}' "${SANITIZED_GTF}"
			gffread -w "${NUCL_TRANSCRIPTS_FA}" -g "${GENOME_REF}" "${nucl_gtf_tmp}"
			gffread -w "${MITO_TRANSCRIPTS_FA}" -g "${GENOME_REF}" "${mito_gtf_tmp}"

		elif [[ "$GENOME_TYPE" == "nuclear" ]]; then
			echo_green "Creating nuclear-only transcript FASTA file..."
			gffread -w "${NUCL_TRANSCRIPTS_FA}" -g "${GENOME_REF}" "${SANITIZED_GTF}"
		elif [[ "$GENOME_TYPE" == "mito" ]]; then
			echo_green "Creating mitochondrial-only transcript FASTA file..."
			gffread -w "${MITO_TRANSCRIPTS_FA}" -g "${GENOME_REF}" "${SANITIZED_GTF}"
		fi

		run_transdecoder() {
			local input_fasta="$1"; local genetic_code="$2"; local output_prefix="$3"
			local work_dir="${FUNCTIONAL_DIR}/${output_prefix}_transdecoder_work"
			local genetic_code_cap=${genetic_code^} # Capitalize first letter
			echo_green "Running TransDecoder for '${output_prefix}' with genetic code '${genetic_code_cap}'..."
			mkdir -p "${work_dir}"
			TransDecoder.LongOrfs -t "${input_fasta}" --genetic_code "${genetic_code_cap}" -O "${work_dir}"
			TransDecoder.Predict -t "${input_fasta}" --genetic_code "${genetic_code_cap}" --single_best_only -O "${work_dir}"
		}

		if [[ "$GENOME_TYPE" == "nuclear" || "$GENOME_TYPE" == "mixed" ]]; then
			[ -s "${NUCL_TRANSCRIPTS_FA}" ] && run_transdecoder "${NUCL_TRANSCRIPTS_FA}" "${GENETIC_CODE_NUCL}" "nuclear"
		fi
		if [[ "$GENOME_TYPE" == "mito" || "$GENOME_TYPE" == "mixed" ]]; then
			[ -s "${MITO_TRANSCRIPTS_FA}" ] && run_transdecoder "${MITO_TRANSCRIPTS_FA}" "${GENETIC_CODE_MITO}" "mito"
		fi
		
		# Consolidate peptide predictions
		cat "${FUNCTIONAL_DIR}"/*/transcripts.fa.transdecoder.pep 2>/dev/null > "${LONGEST_ORFS_PEP}"
		local LONGEST_ORFS_PEP_CLEAN="${FUNCTIONAL_DIR}/longest_orfs.pep.clean"
		sed 's/\*//g' "${LONGEST_ORFS_PEP}" > "${LONGEST_ORFS_PEP_CLEAN}"

		# Download databases if needed
		if [ -z "${BLAST_DB_SwissProt}" ]; then
			BLAST_DB_SwissProt="${FUNCTIONAL_DIR}/blast_dbs/swissprot"
			if [ ! -f "${BLAST_DB_SwissProt}.pin" ]; then
				echo_green  "SwissProt database not found. Downloading..."
				mkdir -p "${FUNCTIONAL_DIR}/blast_dbs"
				wget -qO - "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz" | gunzip -c > "${FUNCTIONAL_DIR}/blast_dbs/swissprot.fasta"
				makeblastdb -in "${FUNCTIONAL_DIR}/blast_dbs/swissprot.fasta" -dbtype prot -out "${BLAST_DB_SwissProt}"
			fi
		fi

		if [ -z "${PFAM_DB}" ]; then
			PFAM_DB="${FUNCTIONAL_DIR}/pfam_db"
			if [ ! -f "${PFAM_DB}/Pfam-A.hmm.h3m" ]; then
				echo_green  "Pfam database not found. Downloading..."
				mkdir -p "${PFAM_DB}"
				wget -qO - "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz" | gunzip -c > "${PFAM_DB}/Pfam-A.hmm"
				hmmpress "${PFAM_DB}/Pfam-A.hmm"
			fi
		fi

		# Run homology searches if requested and peptides were predicted
		if [ -s "${LONGEST_ORFS_PEP_CLEAN}" ]; then
			[[ "${FUNCTIONAL_METHODS}" == *"BLASTp"* ]] && echo_green "Running BLASTp..." && blastp -query "${LONGEST_ORFS_PEP_CLEAN}" -db "${BLAST_DB_SwissProt}" -outfmt 6 -num_threads "${THREADS}" -out "${FUNCTIONAL_DIR}/blastp_results.out"
			[[ "${FUNCTIONAL_METHODS}" == *"PFAM"* ]] && echo_green "Running PFAM Scan..." && hmmscan --cpu "${THREADS}" --domtblout "${FUNCTIONAL_DIR}/pfam_results.out" "${PFAM_DB}/Pfam-A.hmm" "${LONGEST_ORFS_PEP_CLEAN}"
			[[ "${FUNCTIONAL_METHODS}" == *"INTERPRO"* ]] && echo_green "Running InterProScan..." && run_interproscan "${LONGEST_ORFS_PEP_CLEAN}" "${FUNCTIONAL_DIR}/interpro.tsv" "${THREADS}"
		fi
		
		local ALL_TRANSCRIPTS_FA="${FUNCTIONAL_DIR}/all_transcripts.fa"
		cat "${NUCL_TRANSCRIPTS_FA}" "${MITO_TRANSCRIPTS_FA}" 2>/dev/null > "${ALL_TRANSCRIPTS_FA}"
		if [ -s "${ALL_TRANSCRIPTS_FA}" ]; then
			[[ "${FUNCTIONAL_METHODS}" == *"BLASTx"* ]] && echo_green "Running BLASTx..." && blastx -query "${ALL_TRANSCRIPTS_FA}" -db "${BLAST_DB_SwissProt}" -outfmt 6 -num_threads "${THREADS}" -out "${FUNCTIONAL_DIR}/blastx_results.out"
		fi
		echo_blue  "Step 9 Completed: Functional Annotation"
	}

	# ---------------------------
	# Step 10: Integrate Functional Annotation (R Script)
	# ---------------------------
	step10_integrate_functional_annotation() {
		echo_blue "Starting Step 10: Integrate Functional Annotation (R Script)"
		
		local R_CMD="Rscript functional_annotation.R"
		R_CMD+=" --annotation_dir \"${ANNOTATION_DIR}\""
		R_CMD+=" --functional_dir \"${FUNCTIONAL_DIR}\""
		R_CMD+=" --output_dir \"${OUTPUT_DIR}\""
		[ -n "$GENOME_GTF" ] && R_CMD+=" --genome_ref \"$GENOME_GTF\""
		[ -n "$MAX_DISTANCE" ] && R_CMD+=" --maxDistance $MAX_DISTANCE"
		[ -n "$LARGE_INTRON_THRESHOLD" ] && R_CMD+=" --largeIntronThreshold $LARGE_INTRON_THRESHOLD"
		[ -n "$BLAST_EVALUE" ] && R_CMD+=" --blastEvalue $BLAST_EVALUE"
		[ -n "$BLAST_IDENTITY" ] && R_CMD+=" --blastIdentity $BLAST_IDENTITY"
		[ -n "$PFAM_CEVALUE" ] && R_CMD+=" --pfamCEvalue $PFAM_CEVALUE"
		[ -n "$PFAM_BITSCORE" ] && R_CMD+=" --pfamBitScore $PFAM_BITSCORE"
		[ -n "$PFAM_COVERAGE" ] && R_CMD+=" --pfamCoverage $PFAM_COVERAGE"
		[ -n "$INTERPRO_EVALUE" ] && R_CMD+=" --interproEvalue $INTERPRO_EVALUE"

		echo_green "Executing R script..."
		
		eval "${R_CMD}" 2>&1 | tee "${LOGS_DIR}/functional_annotation.log"
		
		if [ ${PIPESTATUS[0]} -ne 0 ]; then
			echo_red "Error: Functional annotation R script failed. Check \"${OUTPUT_PATH}/logs/functional_annotation.log\" for details."
			exit 1
		else
			echo_green "Functional annotation R script completed successfully."
		fi

		echo_blue "Step 10 Completed."
	}
	
	### MAIN SCRIPT LOGIC ###
	
	# ---------------------------
	# Argument Parsing and Default Value Setup
	# ---------------------------
	OUTPUT_DIR="./outputDir"
	THREADS=8
	STRINGTIE2_COVERAGE=1
	STRINGTIE2_ABUNDANCE="0.01"
	MIN_ORF_LENGTH=100
	MAX_EXON_LENGTH=10000
	MAX_TRANSCRIPT_LENGTH=100000
	TRIM_QUAL=20
	TRIM_LEN=20
	MAX_DISTANCE=1000
	LARGE_INTRON_THRESHOLD=100000
	BLAST_EVALUE="1e-5"
	BLAST_IDENTITY=25
	PFAM_CEVALUE="1e-5"
	PFAM_BITSCORE=10
	PFAM_COVERAGE=50
	INTERPRO_EVALUE="1e-5"
	FUNCTIONAL_METHODS="BLASTp,BLASTx,PFAM,INTERPRO"
	DEFAULT_CONDA_SHORT="2.1.1"
	DEFAULT_CONDA_MIX="2.2.1"
	STRINGTIE_SHORT_VERSION="${DEFAULT_CONDA_SHORT}"
	STRINGTIE_MIX_VERSION="${DEFAULT_CONDA_MIX}"
	STRINGTIE_VERSION=""
	GENOME_TYPE_USER="" # Store user-specified genome type
	MITO_PATTERN="[Mm]ito|[Mm]itochondria|mtDNA"
	GENETIC_CODE_NUCL="Universal"
	GENETIC_CODE_MITO="Mitochondrial-Vertebrates"
	LONGREAD_MIN_LEN=1000
	LONGREAD_KEEP_PERCENT=95
	RUN_ALL=false
	STEPS=()
	
	# Use a standard while/case loop for argument parsing
	while [[ "$#" -gt 0 ]]; do
		case $1 in
			--genomeDir) GENOME_DIR="$2"; shift 2 ;;
			--rrnaRef) RRNA_REF="$2"; shift 2 ;;
			--genomeRef) GENOME_REF="$2"; shift 2 ;;
			--genomeGTF) GENOME_GTF="$2"; shift 2 ;;
			--alignDir) ALIGN_DIR="$2"; shift 2 ;;
			--finalGTF) finalGTF="$2"; shift 2 ;;
			--SR_RB_gtf_dir) SR_RB_GTF_DIR="$2"; shift 2 ;;
			--SR_DN_gtf_dir) SR_DN_GTF_DIR="$2"; shift 2 ;;
			--MR_RB_gtf_dir) MR_RB_GTF_DIR="$2"; shift 2 ;;
			--MR_DN_gtf_dir) MR_DN_GTF_DIR="$2"; shift 2 ;;
			--blastDB_SwissProt) BLAST_DB_SwissProt="$2"; shift 2 ;;
			--pfamDB) PFAM_DB="$2"; shift 2 ;;
			--genomeType) GENOME_TYPE_USER="$2"; shift 2;;
			--mitoPattern) MITO_PATTERN="$2"; shift 2;;
			--geneticCodeNucl) GENETIC_CODE_NUCL="$2"; shift 2;;
			--geneticCodeMito) GENETIC_CODE_MITO="$2"; shift 2;;
			--dataDir) DATA_DIR="$2"; shift 2 ;;
			--outputDir) OUTPUT_DIR="$2"; shift 2 ;;
			--threads) THREADS="$2"; shift 2 ;;
			--Stringtie2minReadCoverage) STRINGTIE2_COVERAGE="$2"; shift 2 ;;
			--Stringtie2minIsoformAbundance)  STRINGTIE2_ABUNDANCE="$2"; shift 2 ;;
			--minOrfLength) MIN_ORF_LENGTH="$2"; shift 2 ;;
			--maxExonLength) MAX_EXON_LENGTH="$2"; shift 2 ;;
			--maxTranscriptLength) MAX_TRANSCRIPT_LENGTH="$2"; shift 2 ;;
			--steps) IFS=',' read -ra STEPS <<< "$2"; shift 2 ;;
			--all) RUN_ALL=true; shift ;;
			--functionalMethods) FUNCTIONAL_METHODS="$2"; shift 2 ;;
			--stringtie) STRINGTIE_VERSION="$2"; shift 2 ;;
			--stringtie_short) STRINGTIE_SHORT_VERSION="$2"; shift 2 ;;
			--stringtie_mix) STRINGTIE_MIX_VERSION="$2"; shift 2 ;;
			--trimQual)   TRIM_QUAL="$2";  shift 2 ;;
			--trimAdapter) TRIM_ADAPTER="$2"; shift 2 ;;
			--trimGzip)   TRIM_GZIP="--gzip"; shift ;;
			--trimLen)    TRIM_LEN="$2";   shift 2 ;;
			--maxDistance) MAX_DISTANCE="$2"; shift 2 ;;
			--largeIntronThreshold) LARGE_INTRON_THRESHOLD="$2"; shift 2 ;;
			--blastEvalue) BLAST_EVALUE="$2"; shift 2 ;;
			--blastIdentity) BLAST_IDENTITY="$2"; shift 2 ;;
			--pfamCEvalue) PFAM_CEVALUE="$2"; shift 2 ;;
			--pfamBitScore) PFAM_BITSCORE="$2"; shift 2 ;;
			--pfamCoverage) PFAM_COVERAGE="$2"; shift 2 ;;
			--interproEvalue) INTERPRO_EVALUE="$2"; shift 2 ;;
			--longread_min_len) LONGREAD_MIN_LEN="$2"; shift 2 ;;
			--longread_keep_percent) LONGREAD_KEEP_PERCENT="$2"; shift 2 ;;
			--longread_target_bases) LONGREAD_TARGET_BASES="$2"; shift 2 ;;
			-h|--help) show_help ;;
			*) echo "Unknown parameter passed: $1"; show_help ;;
		esac
	done

	# Post-processing for --stringtie override
	if [ -n "$STRINGTIE_VERSION" ]; then
		[ "$STRINGTIE_SHORT_VERSION" = "$DEFAULT_CONDA_SHORT" ] && STRINGTIE_SHORT_VERSION="$STRINGTIE_VERSION"
		[ "$STRINGTIE_MIX_VERSION" = "$DEFAULT_CONDA_MIX" ] && STRINGTIE_MIX_VERSION="$STRINGTIE_VERSION"
	fi
	
	# Determine steps to run
	STEP_ORDER=(0 1 2 3 4 5 6 7 8 9 10)
	if [ "$RUN_ALL" = true ]; then
		STEPS_TO_RUN=("${STEP_ORDER[@]}")
	elif [ ${#STEPS[@]} -eq 0 ]; then
		echo_red "Error: You must specify --all or provide a list with --steps."
		show_help
	else
		STEPS_TO_RUN=("${STEPS[@]}")
	fi

	# Let user know the host path for the output
	#HOST_OUTPUT_PATH="${OUTPUT_PATH:-${OUTPUT_DIR}}"
	echo_green "Pipeline output will be written to host directory: \"${OUTPUT_PATH}\""

	# ---------------------------
	# Setup and Pre-flight Checks
	# ---------------------------
	validate_parameters
	check_dependencies
	setup_directories
	setup_conda_environments
	
	# Define Sample Names
	SHORT_ONLY_DIR="${DATA_DIR}/short_reads"
	SHORT_ONLY_SAMPLES=()
	if [ -d "${SHORT_ONLY_DIR}" ]; then
		mapfile -t SHORT_ONLY_SAMPLES < <(find "${SHORT_ONLY_DIR}" -maxdepth 1 -type f \( -name "*1.f*q*" -o -name "*2.f*q*" \) -printf "%f\n" | sed -E 's/(_R)?[12]\..*//;s/(_)[12]\..*//' | sort | uniq)
	fi
	MIX_DIR="${DATA_DIR}/mix_reads"; MIX_SHORT_DIR="${MIX_DIR}/short_reads"; MIX_LONG_DIR="${MIX_DIR}/long_reads"
	MIX_SAMPLES=()
	if [ -d "${MIX_SHORT_DIR}" ]; then
		mapfile -t MIX_SAMPLES < <(find "${MIX_SHORT_DIR}" -maxdepth 1 -type f \( -name "*1.f*q*" -o -name "*2.f*q*" \) -printf "%f\n" | sed -E 's/(_R)?[12]\..*//;s/(_)[12]\..*//' | sort | uniq)
	fi
	[ "${#SHORT_ONLY_SAMPLES[@]}" -gt 0 ] && echo_green "Short-read-only Samples detected: ${SHORT_ONLY_SAMPLES[*]}"
	[ "${#MIX_SAMPLES[@]}" -gt 0 ] && echo_green "Mixed-read Samples detected: ${MIX_SAMPLES[*]}"

	# Define Step Functions Mapping
	declare -A STEP_FUNCTIONS=(
		[0]=step0_preprocess_trimming    [1]=step1_rrna_removal
		[2]=step2_read_alignment         [3]=step3_gene_transcript_assembly
		[4]=step4_merging_transcripts_I  [5]=step5_merging_transcripts_II
		[6]=step6_filter_transcripts     [7]=step7_isoform_comparison
		[8]=step8_gtf_correction         [9]=step9_functional_annotation
		[10]=step10_integrate_functional_annotation
	)

	echo_blue "Beginning pipeline execution..."
	for STEP in "${STEP_ORDER[@]}"; do
		if contains "$STEP" "${STEPS_TO_RUN[@]}"; then
			
			# --- Find required input file for the current step ---
			INPUT_GTF="" # Reset for each step
			if [ "$STEP" -ge 6 ]; then
				if [ -n "${finalGTF}" ] && [ -f "${finalGTF}" ]; then
					INPUT_GTF="${finalGTF}"
				elif [ -f "${ANNOTATION_DIR}/final_annotation.gtf" ]; then
					INPUT_GTF="${ANNOTATION_DIR}/final_annotation.gtf"
				elif [ -f "${ANNOTATION_DIR}/corrected_with_introns.gtf" ]; then
					INPUT_GTF="${ANNOTATION_DIR}/corrected_with_introns.gtf"
				elif [ -f "${MERGE_DIR}/filtered_annotation.gtf" ]; then
					INPUT_GTF="${MERGE_DIR}/filtered_annotation.gtf"
				elif [ -f "${MERGE_DIR}/prefinal_annotation.gtf" ]; then
					INPUT_GTF="${MERGE_DIR}/prefinal_annotation.gtf"
				fi
			fi

			# --- Conditional Validation ---
			case "$STEP" in
				0)
					[ -z "${DATA_DIR}" ] && { echo_red "Error: --dataDir is required for Step 0."; exit 1; }
					;;
				1)
					[ -z "${DATA_DIR}" ] && { echo_red "Error: --dataDir is required for Step 1."; exit 1; }
					[ -z "${RRNA_REF}" ] && { echo_green "Warning: --rrnaRef not provided. Skipping Step 1."; continue 2; }
					;;
				2)
					[ -z "${DATA_DIR}" ] && { echo_red "Error: --dataDir is required for Step 2."; exit 1; }
					[ -z "${GENOME_REF}" ] && { echo_red "Error: --genomeRef is required for Step 2."; exit 1; }
					;;
				3)
					[ -z "${ALIGN_DIR}" ] && [ ! -d "${ALIGN_DIR_INTERNAL}" ] && { echo_red "Error: --alignDir must be provided for Step 3 if not running from Step 2."; exit 1; }
					;;
				4|5)
					# These steps can be skipped if a suitable downstream file exists
					[ -n "$INPUT_GTF" ] && { echo_green "Found downstream GTF \"$INPUT_GTF\". Skipping Step $STEP."; continue 2; }
					[ -z "${SR_RB_GTF_DIR}" ] && [ -z "${SR_DN_GTF_DIR}" ] && [ -z "${MR_RB_GTF_DIR}" ] && [ -z "${MR_DN_GTF_DIR}" ] && \
					[ ! -d "${SR_RB_GTF_DIR_INTERNAL}" ] && [ ! -d "${SR_DN_GTF_DIR_INTERNAL}" ] && [ ! -d "${MR_RB_GTF_DIR_INTERNAL}" ] && \
					[ ! -d "${MR_DN_GTF_DIR_INTERNAL}" ] && [ ! -n "$INPUT_GTF" ] && { echo_red "Error: GTF directories not found for merging steps."; exit 1; }
					;;
				6|7|8)
					[ -z "$INPUT_GTF" ] && { echo_red "Error: A valid input GTF is required for Step $STEP. Provide with --finalGTF or run previous steps."; exit 1; }
					;;
				9|10)
					[ -z "$INPUT_GTF" ] && { echo_red "Error: A valid input GTF is required for Step $STEP. Provide with --finalGTF or run previous steps."; exit 1; }
					[ -z "${GENOME_REF}" ] && { echo_red "Error: --genomeRef must be provided for Step $STEP."; exit 1; }
					;;
			esac

			# --- Execute Step ---
			STEP_FUNCTION="${STEP_FUNCTIONS[$STEP]}"
			echo_green "Executing Step $STEP: ${STEP_FUNCTION}"
			# Pass the determined input GTF to the relevant steps
			if [ "$STEP" -ge 6 ]; then
				"${STEP_FUNCTION}" "${INPUT_GTF}"
			else
				"${STEP_FUNCTION}"
			fi
		fi
	done
	# for STEP in "${STEP_ORDER[@]}"; do
		# if contains "$STEP" "${STEPS_TO_RUN[@]}"; then
			# # --- Conditional Validation ---
			# case "$STEP" in
				# 0|1|2)
					# [ -z "${DATA_DIR}" ] && { echo_red "Error: --dataDir is required for Step ${STEP}."; exit 1; }
					# [ ! -d "${DATA_DIR}/short_reads" ] && [ ! -d "${DATA_DIR}/mix_reads" ] && { echo_red "Error: --dataDir must contain 'short_reads' and/or 'mix_reads'."; exit 1; }
					# ;;
				# 3)
					# [ -z "${ALIGN_DIR}" ] && [ ! -d "${ALIGN_DIR_INTERNAL}" ] && { echo_red "Error: --alignDir must be provided for Step 3 if not running from Step 2."; exit 1; }
					# ;;
				# 4|5|6)
				    # [ -z "${SR_RB_GTF_DIR}" ] && [ -z "${SR_DN_GTF_DIR}" ] && [ -z "${MR_RB_GTF_DIR}" ] && [ -z "${MR_DN_GTF_DIR}" ] && [ ! -d "${SR_RB_GTF_DIR_INTERNAL}" ] && [ ! -d "${SR_DN_GTF_DIR_INTERNAL}" ] && [ ! -d "${MR_RB_GTF_DIR_INTERNAL}" ] && [ ! -d "${MR_DN_GTF_DIR_INTERNAL}" ] && { echo_red "Error: GTF directories not found for merging steps."; exit 1; }
					# ;;
				# 7|8|9|10)
					# [ -z "${GENOME_REF}" ] && { echo_red "Error: --genomeRef must be provided for annotation steps (7-10)."; exit 1; }
					# ;;
			# esac

			# # --- Execute Step ---
			# STEP_FUNCTION="${STEP_FUNCTIONS[$STEP]}"
			# echo_green "Executing Step $STEP: ${STEP_FUNCTION}"
			# "${STEP_FUNCTION}"
		# fi
	# done
	
	# =====================================================================
	# Pipeline Completed
	# =====================================================================
	echo_green  "SmedAnno pipeline completed successfully!"

fi
