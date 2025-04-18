#!/bin/bash

# Usage check
if [ "$#" -ne 3 ]; then
  echo "Usage: $0 POP1 POP2 process_id"
  exit 1
fi

# Arguments
POP1=$1
POP2=$2
PROCESS_ID=$3

# Date in DD-M-YY format
DATE=$(date +%d-%-m-%y)

# Template file
TEMPLATE_FILE="run_multithread_gadma-tmp.sh"

# Output filename
OUTPUT_FILE="submit_gadma_${POP1}-${POP2}_p${PROCESS_ID}_${DATE}.sh"

# Replace values
sed -e "s/--job-name=\"\"/--job-name=\"${POP1}-${POP2}-gadma-p${PROCESS_ID}\"/" \
    -e "s/current_gadma-_%A_%a.o/current_gadma-${POP1}-${POP2}_p${PROCESS_ID}_${DATE}_%A_%a.o/" \
    -e "s/current_gadma-_%A_%a.e/current_gadma-${POP1}-${POP2}_p${PROCESS_ID}_${DATE}_%A_%a.e/" \
    -e "s/param_file/params-nomask-${POP1}-${POP2}_process${PROCESS_ID}/" \
    "$TEMPLATE_FILE" > "./hpc_files/${OUTPUT_FILE}"

echo "Customised SLURM script created: $OUTPUT_FILE"
