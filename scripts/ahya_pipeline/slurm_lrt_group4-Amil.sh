#!/bin/bash --login
#SBATCH --job-name="LRT-g4-amil"
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=10:00:00
#SBATCH --account=a_riginos
#SBATCH --partition=general
#SBATCH -o LRT-g4-amil_%A_%a.o
#SBATCH -e LRT-g4-amil_%A_%a.e


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


python lrt_godambe.py group4-Amil_projected0.8 1het hetsc ../data/fs/projected_group4-Amil no \
  --opt_full 0.80351 300 1.37358 14.67336 4.66591 0.02643 0.11572 0.15849 18.45095 9.70126 0.81308 6.04452 0.28779 0.5612 \
  --opt_nested 0.80155 300 1.14557 15.43704 4.42904 0.36292 15.25317 6.68107 0.33435 2.26313 0.42123 \
  --nested_indices 5 6 12

