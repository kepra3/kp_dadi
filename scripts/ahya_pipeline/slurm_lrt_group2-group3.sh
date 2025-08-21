#!/bin/bash --login
#SBATCH --job-name="LRT-g2-g3_hetsym"
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=10:00:00
#SBATCH --account=a_riginos
#SBATCH --partition=general
#SBATCH -o LRT-g2-g3_hetsym_%A_%a.o
#SBATCH -e LRT-g2-g3_hetsym_%A_%a.e


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


python lrt_godambe.py group2-group3_projected0.8 1het 1het_sym ../data/fs/projected_group2-group3 no \
  --opt_full 9.8856 9.65009 1.9536 30.59704 36.38299 8.16505 5.01348 0.01365 6.17711 3.42416 13.09436 14.99969 0.64942 0.29096 \
  --opt_nested 9.5228 11.78008 1.69886 25.50125 25.83647 0.53095 0.53591 0.22488 \
  --nested_indices 7 8 9 10 11 13
