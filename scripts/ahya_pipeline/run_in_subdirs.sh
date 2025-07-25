#!/bin/bash

# Set the root directory where your subdirectories are
BASE_DIR="$1"  # allow user to pass a directory argument
SCRIPT_TO_RUN="$2"  # allow user to pass the script to run in each subdir

if [[ -z "$BASE_DIR" || -z "$SCRIPT_TO_RUN" ]]; then
    echo "Usage: $0 /path/to/base_dir /path/to/script_to_run.sh"
    exit 1
fi

# Make sure script is executable
chmod +x "$SCRIPT_TO_RUN"

# Loop through all subdirectories of BASE_DIR
find "$BASE_DIR" -mindepth 1 -maxdepth 1 -type d | while read -r SUBDIR; do
    echo "Running script in: $SUBDIR"
    (cd "$SUBDIR" && "$SCRIPT_TO_RUN")
done