#!/bin/bash --login
#SBATCH --job-name="LRT-g1-amil_hetsc-og"
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=10:00:00
#SBATCH --account=a_riginos
#SBATCH --partition=general
#SBATCH -o LRT-g1-amil_hetsc-og_%A_%a.o
#SBATCH -e LRT-g1-amil_hetsc-og_%A_%a.e


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
  --opt_full 0.32369 200 0.6451 6.08667 4.19974 0.01749 0.11777 0.09341 12.2419 4.22317 1.10576 10 0.14524 0.57673 \
  --opt_nested 0.39307 300 0.67909 7.38505 4.3923 0.05602 11.65843 3.60885 1.73737 14.99999 0.54413 \
  --nested_indices 5 6 12
