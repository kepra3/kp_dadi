#!/bin/bash
#
# Script to plot residuals for 1het model using parameters from two different CSV files
# Author: Generated for kp_dadi analysis
# Usage: ./plot_1het_residuals.sh

set -e  # Exit on any error

# Load conda and activate environment
module load anaconda3/2023.09-0
source $EBROOTANACONDA3/etc/profile.d/conda.sh
conda activate gadma_env
export PATH="/home/uqkprat2/.conda/envs/gadma_env/bin:$PATH"


# Define file paths
csv_new="../../results/proj_1het_best_models_dadi_params_2025-08-20_11-02-17.csv"
csv_og="../../results/proj_1het_best_models_dadi_params_2025-08-20_11-02-47.csv"
plot_script="../compare_model.py"

# Check if files exist
if [[ ! -f "$csv_new" ]]; then
    echo "Error: $csv_new not found!"
    exit 1
fi

if [[ ! -f "$csv_og" ]]; then
    echo "Error: $csv_og not found!"
    exit 1
fi

if [[ ! -f "$plot_script" ]]; then
    echo "Error: $plot_script not found!"
    exit 1
fi

echo "Starting residual plotting for 1het model..."

# Function to extract parameters from CSV and run plotting
plot_residuals() {
    local csv_file="$1"
    local suffix="$2"
    
    echo "Processing $csv_file with suffix $suffix"
    
    # Skip header and process each line
    tail -n +2 "$csv_file" | while IFS=',' read -r line; do
        # Remove quotes and parse the line
        line=$(echo "$line" | sed 's/"//g')
        
        # Parse fields
        group=$(echo "$line" | cut -d',' -f1)
        run=$(echo "$line" | cut -d',' -f2)
        loglik=$(echo "$line" | cut -d',' -f3)
        
        # Extract parameters from column 5 onwards (skip Groups, Run, LogLikelihood, theta)
        params=$(echo "$line" | cut -d',' -f5- | tr ',' ' ')
        
        # Define the fs file path
        fs_file="${group}_projected0.8"
        
        echo "Plotting residuals for $group ($run) with suffix $suffix"
        echo "Parameters: $params"
        
        # Run the plotting script
        # Using 'none' for mask, 'folded' for fold type, and 3.0 for resid_range
        cd ../../scripts
        python compare_model.py "$fs_file" "1het" "none" "folded" 3.0 -o $params
        
        # Move/rename the output files to include the suffix
        cd ../plots/models/
        if [[ -f "${fs_file}_1het_none_residual.pdf" ]]; then
            mv "${fs_file}_1het_none_residual.pdf" "${fs_file}_1het_${suffix}_residual.pdf"
            echo "Created: ${fs_file}_1het_${suffix}_residual.pdf"
        fi
        
        if [[ -f "${fs_file}_1het_none_all.pdf" ]]; then
            mv "${fs_file}_1het_none_all.pdf" "${fs_file}_1het_${suffix}_all.pdf"
            echo "Created: ${fs_file}_1het_${suffix}_all.pdf"
        fi
        
        if [[ -f "${fs_file}_data.pdf" ]]; then
            # Only rename data file once (for the first run)
            if [[ ! -f "${fs_file}_data_${suffix}.pdf" ]]; then
                cp "${fs_file}_data.pdf" "${fs_file}_data_${suffix}.pdf"
                echo "Created: ${fs_file}_data_${suffix}.pdf"
            fi
        fi
        
        cd ../../scripts/ahya_pipeline/
        
        echo "Completed plotting for $group with suffix $suffix"
        echo "---"
    done
}

# Create plots directory if it doesn't exist
mkdir -p ../../plots/models/

# Plot residuals for both parameter sets
echo "Plotting residuals for 'new' parameters (first CSV file)..."
plot_residuals "$csv_new" "new"

echo "Plotting residuals for 'og' parameters (second CSV file)..."
plot_residuals "$csv_og" "og"

echo "All residual plots completed!"
echo "Check ../../plots/models/ for the generated PDF files"
