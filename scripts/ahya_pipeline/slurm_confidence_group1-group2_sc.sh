#!/bin/bash --login
#SBATCH --job-name="confidence_group1-group2_sc"
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=10:00:00
#SBATCH --account=a_riginos
#SBATCH --partition=general
#SBATCH -o confidence-group1-group2_sc_%A.o
#SBATCH -e confidence-group1-group2_sc_%A.e

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

pop="group1-group2"

echo "Population pair is: $pop"

# Extract optimised parameters from the results file (starting from 5th column)
results_file="../results/proj_hetc_models_dadi_params.csv"

if [ ! -f "$results_file" ]; then
    echo "Error: Results file $results_file not found!"
    exit 1
fi

echo "Reading optimised parameters from: $results_file"

# Extract the header to find column positions
header=$(head -n 1 "$results_file")
echo "Header: $header"

# Extract parameters for the specified population pair
# Fixed to handle quoted CSV values
opt=$(awk -F',' -v pop="$pop" '
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
' "$results_file")

if [ -z "$opt" ]; then
    echo "Error: No parameters found for population pair: $pop"
    echo "Available population pairs in file:"
    awk -F',' 'NR>1 {gsub(/"/, "", $1); print $1}' "$results_file" | sort | uniq
    exit 1
fi

echo "Optimised params are: $opt"

for i in 1 0.01 0.001 0.0001 0.00001 0.000001; do
    echo "Running confidence intervals for hetsc with epsilon ${i} for ${pop}"
    python confidence_intervals.py "../data/fs/${pop}_projected0.8.fs" "../data/fs/projected_${pop}/" \
        GIM  hetsc ${i} -o ${opt}
done
