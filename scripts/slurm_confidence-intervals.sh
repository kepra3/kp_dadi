#!/bin/bash --login
#SBATCH --job-name="confidence_g2-g4"
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=3:00:00
#SBATCH --account=a_riginos
#SBATCH --partition=general
#SBATCH -o confidence-g2-g4_%A_%a.o
#SBATCH -e confidence-g2-g4_%A_%a.e


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

echo "Population pair is: $pop"
echo "Optimised params are: $opt"

for i in 0.0001; do
    echo "Running confidence intervals for 1het with epsilon ${i} for ${pop}"
    python confidence_intervals.py "../data/fs/${pop}_projected0.8.fs" "../data/fs/projected_${pop}/" \
        GIM  1het ${i} -o ${opt}
done
