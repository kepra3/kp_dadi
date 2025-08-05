#!/bin/bash --login
#SBATCH --job-name="hetsc_reformat-results"
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=2G
#SBATCH --time=24:00:00
#SBATCH --account=a_riginos
#SBATCH --partition=general
#SBATCH -o reformat_results-proj_hetsc_%A_%a.o
#SBATCH -e reformat_results-proj_hetsc_%A_%a.e


# module load python/3.9.5-gcccore-10.3.0
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

cd /scratch/user/uqkprat2/analysis/kp_dadi/scripts/ahya_pipeline

python reformat-gadma-results_custom.py '../../results/gadma/proj_hetsc' '../../results/gadma_proj_hetsc_results_combined.txt' '../custom_model_het_sc.py' 0.8