#!/bin/bash

# Usage: ./extract_best_models.sh
# Searches all subdirectories for output files and extracts final best models.

OUTPUT_FILE="all_extracted_models.txt"
> "$OUTPUT_FILE"  # clear output file at the start

find . -type f -name "GADMA.log" | while read -r FILE; do
    LAST_LINE=$(awk '/^All best by log-likelihood models/ { line=NR } END { print line }' "$FILE")
    
    if [ -n "$LAST_LINE" ]; then
        echo "----- $FILE -----" >> "$OUTPUT_FILE"
        awk -v start="$LAST_LINE" 'NR >= start && NR <= start + 6' "$FILE" >> "$OUTPUT_FILE"
        echo "" >> "$OUTPUT_FILE"
    fi
done

echo "Done! Extracted models saved to $OUTPUT_FILE"
