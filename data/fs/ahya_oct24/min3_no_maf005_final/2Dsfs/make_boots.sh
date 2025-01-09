#!/bin/bash

# Check if the correct number of arguments is provided
if [ "$#" -ne 6 ]; then
    echo "Usage: $0 input_file.sfs integer1 integer2 folded/unfolded \"pop1\" \"pop2\""
    exit 1
fi

input_file="$1"
integer1="$2"
integer2="$3"
fold_status="$4"
pop1="$5"
pop2="$6"

# Validate that integer1 and integer2 are integers
if ! [[ "$integer1" =~ ^[0-9]+$ ]] || ! [[ "$integer2" =~ ^[0-9]+$ ]]; then
    echo "Error: integer1 and integer2 must be integers."
    exit 1
fi

# Validate fold_status
if [[ "$fold_status" != "folded" && "$fold_status" != "unfolded" ]]; then
    echo "Error: fold_status must be either 'folded' or 'unfolded'."
    exit 1
fi

base_name="${input_file%.sfs}"
output_dir="${base_name}_boots"

# Ensure the input file exists
if [ ! -f "$input_file" ]; then
    echo "Error: File '$input_file' not found."
    exit 1
fi

# Create the output directory if it doesn't exist
if [ ! -d "$output_dir" ]; then
    mkdir "$output_dir"
    echo "Created directory: $output_dir"
fi

# Initialize line counter
line_number=1

# Read the input file line by line
while IFS= read -r line; do
    # Define the output file name within the output directory
    output_file="${output_dir}/${base_name}_${line_number}_boot.fs"

    # Write the parameters and the line to the output file
    {
        echo "$integer1 $integer2 $fold_status \"$pop1\" \"$pop2\""
        echo "$line"
    } > "$output_file"

    # Increment the line counter
    ((line_number++))
done < "$input_file"

echo "Processing complete. Created $((line_number - 1)) files in '$output_dir'."