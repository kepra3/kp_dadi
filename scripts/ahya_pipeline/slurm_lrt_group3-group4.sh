#!/bin/bash --login
#SBATCH --job-name="LRT-g3-g4_hetsym"
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=10:00:00
#SBATCH --account=a_riginos
#SBATCH --partition=general
#SBATCH -o LRT-g3-g4_hetsym_%A_%a.o
#SBATCH -e LRT-g3-g4_hetsym_%A_%a.e


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


python lrt_godambe.py group3-group4_projected0.8 1het 1het_sym ../data/fs/projected_group3-group4 no \
  --opt_full 6.01079 15.33701 1.66028 19.81554 37.17038 11.54377 5.6171 0.02302 4.2244 7.41015 11.63298 4.28013 0.54319 0.16216 \
  --opt_nested 5.42866 34.8866 1.37454 14.3881 24.86756 1.52713 0.50544 0.1912 \
  --nested_indices 7 8 9 10 11 13

