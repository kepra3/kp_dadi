#!/bin/bash --login
#SBATCH --job-name="LRT-g1-amil"
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=10:00:00
#SBATCH --account=a_riginos
#SBATCH --partition=general
#SBATCH -o LRT-g1-amil_%A_%a.o
#SBATCH -e LRT-g1-amil_%A_%a.e


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

cd /scratch/user/uqkprat2/analysis/kp_dadi/scripts

#pop=$1
#opt=$2
#nested=$3
#indices=$4

#echo "Population pair is: $pop"
#echo "Optimised params are: $opt"


python lrt_godambe.py group1-Amil_projected0.8 1het hetsc ../data/fs/projected_group1-Amil no \
  --opt_full 0.29053 300 0.60554 4.98481 3.92607 0.02339 0.25356 0.14269 12.10501 4.80591 0.70558 6.31901 0.35056 0.56677 \
  --opt_nested 0.39307 300 0.67909 7.38505 4.3923 0.05602 11.65843 3.60885 1.73737 14.99999 0.54413 \
  --nested_indices 5 6 12

