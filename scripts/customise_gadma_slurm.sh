#!/bin/bash

# Usage check
if [ "$#" -ne 4 ]; then
  echo "Usage: $0 POP1 POP2 process_id param_file"
  exit 1
fi

# Arguments
POP1=$1
POP2=$2
PROCESS_ID=$3
PARAM_FILE=$4

# Date in DD-M-YY format
DATE=$(date +%d-%-m-%y)

# Template file
TEMPLATE_FILE="run_multithread_gadma-tmp.sh"

# Output filename
OUTPUT_FILE="submit_gadma_${PARAM_FILE}-${POP1}-${POP2}_p${PROCESS_ID}_${DATE}.sh"

# Replace values
sed -e "s/--job-name=\"\"/--job-name=\"gadma_${PARAM_FILE}-${POP1}-${POP2}-p${PROCESS_ID}\"/" \
    -e "s/current_gadma-_%A_%a.o/${PARAM_FILE}_gadma-${POP1}-${POP2}_p${PROCESS_ID}_${DATE}_%A_%a.o/" \
    -e "s/current_gadma-_%A_%a.e/${PARAM_FILE}_gadma-${POP1}-${POP2}_p${PROCESS_ID}_${DATE}_%A_%a.e/" \
    -e "s/param_file/${PARAM_FILE}-${POP1}-${POP2}_process${PROCESS_ID}/" \
    "$TEMPLATE_FILE" > "./hpc_files/${OUTPUT_FILE}"

echo "Customised SLURM script created: $OUTPUT_FILE"
