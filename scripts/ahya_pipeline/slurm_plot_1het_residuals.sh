#!/bin/bash --login
#SBATCH --job-name="plot_1het_residuals"
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=2:00:00
#SBATCH --account=a_riginos
#SBATCH --partition=general
#SBATCH --array=1-20
#SBATCH -o plot_1het_residuals_%A_%a.o
#SBATCH -e plot_1het_residuals_%A_%a.e

module load anaconda3/2023.09-0

source $EBROOTANACONDA3/etc/profile.d/conda.sh
conda activate gadma_env

export PATH="/home/uqkprat2/.conda/envs/gadma_env/bin:$PATH"

echo "Conda envs:"
conda info --envs
echo "Python path:"
which python
echo "Python version:"
python --version
echo "Dadi location:"
python -c "import dadi; print(dadi.__file__)"

cd /scratch/user/uqkprat2/analysis/kp_dadi/scripts

# Define file paths
csv_new="../results/proj_1het_best_models_dadi_params_2025-08-20_11-02-17.csv"
csv_og="../results/proj_1het_best_models_dadi_params_2025-08-20_11-02-47.csv"

# Create plots directory if it doesn't exist
mkdir -p ../plots/models/

# Define arrays for population pairs and suffixes
# We have 10 population pairs, each will be run with both "new" and "og" parameters
# So array indices 1-10 will be "new" and 11-20 will be "og"

populations=("group1-group2" "group1-group3" "group1-group4" "group2-group3" "group2-group4" "group3-group4" "group1-Amil" "group2-Amil" "group3-Amil" "group4-Amil")

# Determine which population and which parameter set based on SLURM_ARRAY_TASK_ID
if [ $SLURM_ARRAY_TASK_ID -le 10 ]; then
    # Use "new" parameters (first CSV file)
    csv_file="$csv_new"
    suffix="new"
    pop_index=$((SLURM_ARRAY_TASK_ID - 1))
else
    # Use "og" parameters (second CSV file)
    csv_file="$csv_og"
    suffix="og"
    pop_index=$((SLURM_ARRAY_TASK_ID - 11))
fi

population=${populations[$pop_index]}

echo "Array task ID: $SLURM_ARRAY_TASK_ID"
echo "Population: $population"
echo "Parameter set: $suffix"
echo "CSV file: $csv_file"

# Extract parameters for this specific population from the CSV file
params=$(awk -F',' -v pop="$population" '
    NR==1 {next}  # Skip header
    {
        # Remove quotes from first column for comparison
        gsub(/"/, "", $1)
    }
    $1 == pop {
        # Extract from 5th column onwards and join with spaces
        result = ""
        for(i=5; i<=NF; i++) {
            if(i > 5) result = result " "
            # Remove quotes from parameter values
            gsub(/"/, "", $i)
            result = result $i
        }
        print result
        exit
    }
' "$csv_file")

if [ -z "$params" ]; then
    echo "Error: No parameters found for population: $population in $csv_file"
    exit 1
fi

echo "Parameters: $params"

# Define the fs file path
fs_file="${population}_projected0.8"

echo "Plotting residuals for $population with suffix $suffix"

# Run the plotting script
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
    # Rename data file with suffix
    mv "${fs_file}_data.pdf" "${fs_file}_data_${suffix}.pdf"
    echo "Created: ${fs_file}_data_${suffix}.pdf"
fi

echo "Completed plotting for $population with suffix $suffix"
