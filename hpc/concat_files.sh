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

# Concatenate all files in the directory into the output file
cat "$directory"/* > "$output_file"

# Print success message
echo "All files in $directory have been concatenated into $output_file"
