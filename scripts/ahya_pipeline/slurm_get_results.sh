#!/bin/bash --login
#SBATCH --job-name="get_results_g1-g4_hetsym"
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=2G
#SBATCH --time=2:00:00
#SBATCH --account=a_riginos
#SBATCH --partition=general
#SBATCH -o get_results-g1-g4-hetsym_%A_%a.o
#SBATCH -e get_results-g1-g4-hetsym_%A_%a.e


# module load python/3.9.5-gcccore-10.3.0 ##### SBATCH --qos=debug
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

./run_in_subdirs.sh '../../results/gadma/proj_hetsym_22-8-25' "/scratch/user/uqkprat2/analysis/kp_dadi/scripts/ahya_pipeline/extract_best_gadma_results.sh"
