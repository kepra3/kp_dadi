#!/bin/bash --login
#SBATCH --job-name="plot_hetsc_hetsym_residuals"
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=2:00:00
#SBATCH --account=a_riginos
#SBATCH --partition=general
#SBATCH --array=1-13
#SBATCH -o plot_hetsc_hetsym_residuals_%A_%a.o
#SBATCH -e plot_hetsc_hetsym_residuals_%A_%a.e

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
csv_hetsc="../results/proj_hetsc_best_models_dadi_params.csv"
csv_hetsym="../results/proj_hetsym_best_models_dadi_params_2025-09-16_10-17-14.csv"

# Create plots directory if it doesn't exist
mkdir -p ../plots/models/

# Define arrays for population pairs and model types
# hetsc populations (7 groups)
hetsc_populations=("group1-group2" "group1-group3" "group1-group4" "group1-Amil" "group2-Amil" "group3-Amil" "group4-Amil")

# hetsym populations (6 groups)  
hetsym_populations=("group1-group2" "group1-group4" "group2-group3" "group2-group4" "group3-group4" "group1-Amil")

# Determine which population and model type based on SLURM_ARRAY_TASK_ID
if [ $SLURM_ARRAY_TASK_ID -le 7 ]; then
    # Use hetsc parameters (tasks 1-7)
    csv_file="$csv_hetsc"
    model_type="hetsc"
    pop_index=$((SLURM_ARRAY_TASK_ID - 1))
    population=${hetsc_populations[$pop_index]}
else
    # Use hetsym parameters (tasks 8-13)
    csv_file="$csv_hetsym"
    model_type="hetsym"
    pop_index=$((SLURM_ARRAY_TASK_ID - 8))
    population=${hetsym_populations[$pop_index]}
fi

echo "Array task ID: $SLURM_ARRAY_TASK_ID"
echo "Population: $population"
echo "Model type: $model_type"
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

echo "Plotting residuals for $population with model $model_type"

# Run the plotting script with appropriate model type
python compare_model.py "$fs_file" "$model_type" "none" "folded" 3.0 -o $params

# Move/rename the output files to include the model type
cd ../plots/models/

if [[ -f "${fs_file}_${model_type}_none_residual.pdf" ]]; then
    mv "${fs_file}_${model_type}_none_residual.pdf" "${fs_file}_${model_type}_residual.pdf"
    echo "Created: ${fs_file}_${model_type}_residual.pdf"
fi

if [[ -f "${fs_file}_${model_type}_none_all.pdf" ]]; then
    mv "${fs_file}_${model_type}_none_all.pdf" "${fs_file}_${model_type}_all.pdf"
    echo "Created: ${fs_file}_${model_type}_all.pdf"
fi

if [[ -f "${fs_file}_data.pdf" ]]; then
    # Rename data file with model type
    mv "${fs_file}_data.pdf" "${fs_file}_data_${model_type}.pdf"
    echo "Created: ${fs_file}_data_${model_type}.pdf"
fi

echo "Completed plotting for $population with model $model_type"