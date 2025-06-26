#!/bin/bash --login
#SBATCH --job-name="reformat-results"   	# job name
#SBATCH --nodes=1               # use 1 node
#SBATCH --ntasks-per-node=1     # use 1 for single and multi core jobs
#SBATCH --cpus-per-task=3       # number of cores per job
#SBATCH --mem=1G               # RAM per job given in megabytes (M), gigabytes (G), or terabytes (T)
#SBATCH --time=24:00:00         # walltime
#SBATCH --account=a_riginos     # group account name
#SBATCH --partition=general     # queue name
#SBATCH -o reformat_results-full_%A_%a.o # standard output
#SBATCH -e reformat_results-full_%A_%a.e # standard error

module load anaconda3/2023.09-0
module load python/3.9.5-gcccore-10.3.0
source $EBROOTANACONDA3/etc/profile.d/conda.sh
conda activate gadma_env

cd /scratch/user/uqkprat2/analysis/kp_dadi/scripts

python reformat-gadma-results.py '../results/gadma/full' '../results/gadma_full_results_combined.txt'