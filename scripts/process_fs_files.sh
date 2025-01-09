#!/bin/bash

# Check if the directory argument is provided
if [ -z "$1" ]; then
    echo "Usage: $0 /path/to/directory"
    exit 1
fi

# Assign the first argument to the directory variable
directory="$1"
pops=$2

# Check if the provided argument is a valid directory
if [ ! -d "$directory" ]; then
    echo "Error: '$directory' is not a valid directory."
    exit 1
fi

# Loop through each .fs file in the specified directory
for file in "$directory"/*.fs; do
    # Check if the file exists to handle the case where no .fs files are found
    if [ -e "$file" ]; then
        # Extract the base name of the file without the directory path
        base_name=$(basename "$file" .fs)
        # Run the Python script with the specified parameters
        python plot_fs.py $file folded projection no -o "projected_$2"
        echo "Processed: ${base_name}"
    else
        echo "No .fs files found in the directory."
        break
    fi
done

