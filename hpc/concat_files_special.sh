#!/bin/bash

# Check if the directory is provided as an argument
if [ $# -ne 1 ]; then
    echo "Usage: $0 <directory>"
    exit 1
fi

# Assign the directory to a variable
directory="$1"

# Check if the argument is a valid directory
if [ ! -d "$directory" ]; then
    echo "Error: $directory is not a valid directory"
    exit 1
fi

# Create the output file name by appending .txt to the directory name
output_file="$(basename "$directory").txt"

# Empty the output file if it exists
> "$output_file"

# Get sorted list of files (optional for predictable order)
files=("$directory"/*)
file_count=${#files[@]}

if [ "$file_count" -eq 0 ]; then
    echo "No files found in $directory"
    exit 1
fi

# Process the first file: append entire contents
cat "${files[0]}" >> "$output_file"

# Process remaining files: extract and append second line only
for ((i=1; i<file_count; i++)); do
    sed -n '2p' "${files[$i]}" >> "$output_file"
done

# Print success message
echo "Output written to $output_file"

