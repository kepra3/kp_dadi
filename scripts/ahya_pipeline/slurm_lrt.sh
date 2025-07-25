#!/bin/bash --login
#SBATCH --job-name="LRT-g2-g4"
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=10:00:00
#SBATCH --account=a_riginos
#SBATCH --partition=general
#SBATCH -o LRT-g2-g4_%A_%a.o
#SBATCH -e LRT-g2-g4_%A_%a.e


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

pop=$1
opt=$2
nested=$3
indices=$4

echo "Population pair is: $pop"
echo "Optimised params are: $opt"


 python lrt_godambe.py ../data/fs/group2-group4_0.8.fs 1het 1het_sym ../data/fs/projected_group2-group4 no \
  --opt_full nu_1, nu_2, t1, nu11, nu12, me1_12, me1_21, t2, nu21, nu22, me2_12, me2_21, P1, P2 \
  --opt_nested nu_1, nu_2, t1, nu11, nu12, me1_12, me1_21, P1 \
  --nested_indices 7 8 9 10 11 13

