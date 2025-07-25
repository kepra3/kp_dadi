#!/bin/bash

# Usage: ./customize_gadma_params.sh <pop1> <pop2> <process_number> <num_repeats> <num_processes> <template_file> <output_file>

# Check for the correct number of arguments
if [ "$#" -ne 7 ]; then
    echo "Usage: $0 <pop1> <pop2> <process_number> <num_repeats> <num_processes> <template_file> <output_file>"
    exit 1
fi

# Assign input arguments to variables
pop1="$1"
pop2="$2"
PROCESS_NUM="$3"
NUM_REPEATS="$4"
NUM_PROCESSES="$5"
TEMPLATE_FILE="$6"
OUTPUT_FILE="$7"

# Get the current date in DD-M-YY format
CURRENT_DATE=$(date +%d-%-m-%y)

# Perform replacements in the output file
sed -e "s/pop1/$pop1/g" \
    -e "s/pop2/$pop2/g" \
    -e "s/processX/process${PROCESS_NUM}_${CURRENT_DATE}/g" \
    -e "s/^Number of repeats: .*/Number of repeats: ${NUM_REPEATS}/" \
    -e "s/^Number of processes: .*/Number of processes: ${NUM_PROCESSES}/" \
    "$TEMPLATE_FILE" > "${OUTPUT_FILE}"

echo "Updated parameter file created: $OUTPUT_FILE"