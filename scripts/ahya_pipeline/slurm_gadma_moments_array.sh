#!/bin/bash --login
#SBATCH --job-name=gadma_moments
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=25G
#SBATCH --time=100:00:00
#SBATCH --account=a_senv_mbos
#SBATCH --partition=general
#SBATCH --array=0-79
#SBATCH -o /scratch/user/uqkprat2/analysis/kp_dadi/logs/gadma_moments/%x_%A_%a.o
#SBATCH -e /scratch/user/uqkprat2/analysis/kp_dadi/logs/gadma_moments/%x_%A_%a.e
#
# GADMA with Engine: moments, same settings as the dadi proj_* runs.
# 10 pairs x 8 processes = 80 tasks; each task = 8 GADMA repeats on 4 processes (as before).
#
# Submit (one array per model), from anywhere:
#   mkdir -p /scratch/user/uqkprat2/analysis/kp_dadi/logs/gadma_moments
#   for m in 1het hetsc hetsym; do
#     sbatch --job-name=gadma_moments_$m --export=ALL,MODEL=$m,RUNTAG=$(date +%d-%-m-%y) slurm_gadma_moments_array.sh
#   done

set -euo pipefail
: "${MODEL:?set MODEL=1het|hetsc|hetsym}"
: "${RUNTAG:?set RUNTAG, e.g. 02-11-26}"
R=/scratch/user/uqkprat2/analysis/kp_dadi
PAIRS=(group1-group2 group1-group3 group1-group4 group1-Amil group2-group3
       group2-group4 group2-Amil group3-group4 group3-Amil group4-Amil)
PAIR=${PAIRS[$((SLURM_ARRAY_TASK_ID / 8))]}
PROC=$((SLURM_ARRAY_TASK_ID % 8 + 1))
POP1=${PAIR%-*}; POP2=${PAIR#*-}

RUNDIR=$R/results/gadma/proj_${MODEL}_moments_${RUNTAG}
mkdir -p "$RUNDIR/params"
PF=$RUNDIR/params/params-proj_${MODEL}_moments-${PAIR}_process${PROC}
sed -e "s/RUNTAG/${RUNTAG}/g" \
    -e "s/pop1/${POP1}/g" -e "s/pop2/${POP2}/g" \
    -e "s/processX/process${PROC}/g" \
    -e "s/^Number of repeats: .*/Number of repeats: 8/" \
    -e "s/^Number of processes: .*/Number of processes: 4/" \
    "$R/scripts/ahya_pipeline/params-proj_${MODEL}_moments" > "$PF"

module load anaconda3/2023.09-0
source $EBROOTANACONDA3/etc/profile.d/conda.sh
conda activate gadma_env
export OMP_NUM_THREADS=1

echo "MODEL=$MODEL PAIR=$PAIR PROC=$PROC params=$PF"
grep -E '^(Output directory|Input file|Engine|Custom filename|Number of)' "$PF"
gadma -p "$PF"
