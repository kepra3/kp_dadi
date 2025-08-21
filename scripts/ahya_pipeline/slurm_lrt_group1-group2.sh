#!/bin/bash --login
#SBATCH --job-name="LRT-g1-g2_hetsc"
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=10:00:00
#SBATCH --account=a_riginos
#SBATCH --partition=general
#SBATCH -o LRT-g1-g2_hetsc_%A_%a.o
#SBATCH -e LRT-g1-g2_hetsc_%A_%a.e


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


python lrt_godambe.py group1-group2_projected0.8 1het hetsc ../data/fs/projected_group1-group2 no  \
  --opt_full 1.68837 16.02573 2.70899 33.49254 32.81391 0.01289 0.44575 0.01283 34.22929 29.41194 15 14.21391 0.35562 0.43557 \
  --opt_nested 1.32826 300 1.22731 22.51267 14.62929 0.41298 22.37374 24.04221 0.57007 1.61637 0.28934 \
  --nested_indices  5 6 12

