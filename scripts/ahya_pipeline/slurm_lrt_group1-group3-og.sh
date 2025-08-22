#!/bin/bash --login
#SBATCH --job-name="LRT-g1-g3_hetsc-og"
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=10:00:00
#SBATCH --account=a_riginos
#SBATCH --partition=general
#SBATCH -o LRT-g1-g3_hetsc-og_%A_%a.o
#SBATCH -e LRT-g1-g3_hetsc-og_%A_%a.e


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


python lrt_godambe.py group1-group3_projected0.8 1het hetsc ../data/fs/projected_group1-group3 no \
  --opt_full 1.50434 24.1785 1.79904 27.34257 22.19025 0.05686 0.08521 0.09195 19.59638 34.76409 2.35863 9.99995 0.17654 0.46391 \
  --opt_nested 1.96862 299.91741 1.29427 28.6625 15.04335 0.16929 14.13291 35.62112 1.79627 2.92929 0.31897 \
  --nested_indices 5 6 12
