#!/bin/bash --login
#SBATCH --job-name="LRT-g1-g4_hetsc"
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=10:00:00
#SBATCH --account=a_riginos
#SBATCH --partition=general
#SBATCH -o LRT-g1-g4_hetsc_%A_%a.o
#SBATCH -e LRT-g1-g4_hetsc_%A_%a.e


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


python lrt_godambe.py group1-group4_projected0.8 1het hetsc ../data/fs/projected_group1-group4 no \
  --opt_full 1.03823 34.30998 1.94183 19.01666 30.18014 0 0.63472 0.01701 14.77849 6.56549 15 3.47854 0.25364 0.3462 \
  --opt_nested 0.91597 299.90144 1.10548 14.37829 18.73503 0.42471 19.05739 19.71916 0.5353 1.95219 0.24297 \
  --nested_indices  5 6 12

