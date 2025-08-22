#!/bin/bash --login
#SBATCH --job-name="LRT-g3-amil_hetsc-og"
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=10:00:00
#SBATCH --account=a_riginos
#SBATCH --partition=general
#SBATCH -o LRT-g3-amil_hetsc-og_%A_%a.o
#SBATCH -e LRT-g3-amil_hetsc-og_%A_%a.e


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

python lrt_godambe.py group3-Amil_projected0.8 1het hetsc ../data/fs/projected_group3-Amil no \
  --opt_full 1 200 1.50213 14.8254 5.98913 0.02281 0.06197 0.09379 17.48163 12.3383 1.59577 10 0.23576 0.59581 \
  --opt_nested 0.73348 300 1.03507 11.7859 4.73052 0.31988 13.18054 6.96225 0.41145 2.50226 0.45673 \
  --nested_indices 5 6 12

