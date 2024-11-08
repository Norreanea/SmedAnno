#!/bin/bash

OUTPUT_DIR="test_data/mix_reads/long_reads/Sample1"

# Create output directory
mkdir -p ${OUTPUT_DIR}
pbsim --strategy wgs\
      --method qshmm\
      --qshmm QSHMM-RSII.model\
      --depth 20\
      --genome test_data/genome.fa\
	  --prefix ${OUTPUT_DIR}
echo "Long-read simulation completed."