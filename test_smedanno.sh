#!/bin/bash

# =================================================================
# SmedAnno Automated Test Suite
# =================================================================
# This script runs a series of tests to validate the SmedAnno pipeline,
# covering different scenarios.
# =================================================================

# Configuration 
# Define colors for output
GREEN='\033[0;32m'
RED='\033[0;31m'
BLUE='\033[1;34m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Base paths relative to the SmedAnno-main directory 
BASE_DIR="../"
GENOME_REF="${BASE_DIR}/SM_genome_reduced/smed_chr4_reduced.fa"
GENOME_GTF="${BASE_DIR}/SM_genome_reduced/SM_HC_filtered_chr4.gtf"
DATA_DIR="${BASE_DIR}/SM_genome_reduced/" # Contains the 'short_reads' folder
ALIGN_DIR="${BASE_DIR}/output_new/alignment/"
PREEXISTING_GTF="${BASE_DIR}/output_new/annotation/corrected_with_introns.gtf"
RUN_SCRIPT="./run_smedanno.sh"

# Helper function
# This function runs a test case, captures its output, and reports success or failure.
# Arguments:
#   $1: Test number (e.g., "1")
#   $2: Test description string
#   $3: Expected outcome ("PASS" or "FAIL")
#   $@: The full run_smedanno.sh command to execute (from the 4th argument onwards)
run_test() {
    local test_num="$1"
    local description="$2"
    local expected_outcome="$3"
    shift 3
    local command_to_run=("$@")

    local output_dir="${BASE_DIR}/test_output_${test_num}"
    local log_file="${BASE_DIR}/test_log_${test_num}.txt"

    echo -e "\n${BLUE}=================================================================${NC}"
    echo -e "${BLUE}▶️  STARTING TEST ${test_num}: ${description}${NC}"
    echo -e "${BLUE}=================================================================${NC}"
    echo "   Output will be in: ${output_dir}"
    echo "   Log will be in:    ${log_file}"
    echo "   Command: ${command_to_run[*]} --outputDir ${output_dir}"
    
    # Add the unique output directory to the command
    command_to_run+=("--outputDir" "${output_dir}")

    # Execute the command and capture exit code
    if "${command_to_run[@]}" &> "${log_file}"; then
        # Command succeeded (exit code 0)
        if [ "$expected_outcome" == "PASS" ]; then
            echo -e "${GREEN}✅ TEST ${test_num} PASSED (Completed Successfully as Expected)${NC}"
        else
            echo -e "${RED}❌ TEST ${test_num} FAILED (Completed Successfully but Was Expected to Fail)${NC}"
        fi
    else
        # Command failed (non-zero exit code)
        if [ "$expected_outcome" == "FAIL" ]; then
            echo -e "${GREEN}✅ TEST ${test_num} PASSED (Failed as Expected)${NC}"
            echo -e "${YELLOW}   (This is a successful test of the validation logic. See log for error details.)${NC}"
        else
            echo -e "${RED}❌ TEST ${test_num} FAILED (Crashed Unexpectedly)${NC}"
            echo -e "${RED}   Please check the log file for details: ${log_file}${NC}"
        fi
    fi
}


# =================================================================
# ▶️  RUNNING TEST SUITE
# =================================================================

# Valid pipeline runs 

# Test 1: Full DE NOVO pipeline from FASTQ
run_test "1" "Full DE NOVO pipeline from FASTQ (steps 0-10)" "PASS" \
    "$RUN_SCRIPT" --all --dataDirShort "$DATA_DIR" --genomeRef "$GENOME_REF" --threads 4

# Test 2: Full REFERENCE-BASED pipeline from FASTQ
run_test "2" "Full REFERENCE-BASED pipeline from FASTQ (steps 0-10)" "PASS" \
    "$RUN_SCRIPT" --all --dataDirShort "$DATA_DIR" --genomeRef "$GENOME_REF" --genomeGTF "$GENOME_GTF" --threads 4

# Test 3: DE NOVO assembly from existing BAM files + change StringTie2 version
run_test "3" "DE NOVO assembly from existing BAMs (steps 3,5-10) and change StringTie2 version" "PASS" \
    "$RUN_SCRIPT" --steps "3,5,6,7,8,9,10" --alignDirShort "$ALIGN_DIR" --genomeRef "$GENOME_REF" --threads 4 --stringtie 3.0.1

# Test 4: REFERENCE-BASED assembly from existing BAM files
run_test "4" "REFERENCE-BASED assembly from existing BAMs (steps 3,4,5-10)" "PASS" \
    "$RUN_SCRIPT" --steps "3,4,5,6,7,8,9,10" --alignDirShort "$ALIGN_DIR" --genomeRef "$GENOME_REF" --genomeGTF "$GENOME_GTF" --threads 4

# Test 5: Annotation-only pipeline from a pre-existing GTF
run_test "5" "Annotation-only from a pre-existing GTF (steps 9,10)" "PASS" \
    "$RUN_SCRIPT" --steps "9,10" --finalGTF "$PREEXISTING_GTF" --genomeRef "$GENOME_REF" --threads 4


# Invalid runs 

# Test 6: Missing --dataDir for alignment
run_test "6" "INVALID: Missing --dataDir for alignment (step 2)" "FAIL" \
    "$RUN_SCRIPT" --steps "2" --genomeRef "$GENOME_REF"

# Test 7: Missing --alignDir for assembly
run_test "7" "INVALID: Missing --alignDir for assembly (step 3)" "FAIL" \
    "$RUN_SCRIPT" --steps "3" --genomeRef "$GENOME_REF"

# Test 8: Missing --genomeRef for functional annotation
run_test "8" "INVALID: Missing --genomeRef for annotation (step 9)" "FAIL" \
    "$RUN_SCRIPT" --steps "9" --finalGTF "$PREEXISTING_GTF"

# Test 9: Invalid value for --threads
run_test "9" "INVALID: Using a non-numeric value for --threads" "FAIL" \
    "$RUN_SCRIPT" --steps "2" --dataDirShort "$DATA_DIR" --genomeRef "$GENOME_REF" --threads "four"

echo -e "\n${BLUE}=================================================================${NC}"
echo -e "${BLUE}🎉 Test Suite Finished! 🎉${NC}"
echo -e "${BLUE}=================================================================${NC}"

