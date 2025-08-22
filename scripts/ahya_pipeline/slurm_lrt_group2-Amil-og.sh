#!/bin/bash --login
#SBATCH --job-name="LRT-g2-amil_hetsc-og"
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=10:00:00
#SBATCH --account=a_riginos
#SBATCH --partition=general
#SBATCH -o LRT-g2-amil_hetsc-og_%A_%a.o
#SBATCH -e LRT-g2-amil_hetsc-og_%A_%a.e


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


python lrt_godambe.py group2-Amil_projected0.8 1het hetsc ../data/fs/projected_group2-Amil no \
  --opt_full 0.81888 200 1.30798 13.30744 5.52148 0.02178 0.04701 0.1554 16.69237 9.66337 0.91669 6.21186 0.15831 0.57265 \
  --opt_nested 0.93373 94.19885 1.24583 14.22307 5.49178 0.14058 14.30639 8.0265 0.94924 5.65778 0.51016 \
  --nested_indices 5 6 12

