#!/bin/bash --login
#SBATCH --job-name="gadma_params-full_2het-group1-Amil-p6"   	# job name
#SBATCH --nodes=1               # use 1 node
#SBATCH --ntasks-per-node=1     # use 1 for single and multi core jobs
#SBATCH --cpus-per-task=8       # number of cores per job
#SBATCH --mem=25G               # RAM per job given in megabytes (M), gigabytes (G), or terabytes (T)
#SBATCH --time=100:00:00         # walltime
#SBATCH --account=a_riginos     # group account name
#SBATCH --partition=general     # queue name
#SBATCH -o params-full_2het_gadma-group1-Amil_p6_05-6-25_%A_%a.o # standard output
#SBATCH -e params-full_2het_gadma-group1-Amil_p6_05-6-25_%A_%a.e # standard error

module load anaconda3/2023.09-0
module load python/3.9.5-gcccore-10.3.0
source $EBROOTANACONDA3/etc/profile.d/conda.sh
conda activate gadma_env

gadma -p params-full_2het-group1-Amil_process6
