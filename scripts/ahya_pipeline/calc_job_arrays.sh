#!/bin/bash

# Set number of parameters (n) and maximum job array size
n=$1  # Example number of parameters
max_array_size=1000  # Slurm's maximum job array size, can be adjusted

# Calculate total combinations (2^n)
total_combinations=$((2**n))

echo "Total combinations: $total_combinations"

# Calculate how many Slurm job arrays are needed
job_arrays_needed=$(( (total_combinations + max_array_size - 1) / max_array_size ))

echo "Job arrays needed: $job_arrays_needed"

# Loop over each job array
for job_id in $(seq 0 $((job_arrays_needed - 1))); do
    # Calculate the start and end indices for this job array
    start_index=$((job_id * max_array_size))
    end_index=$(((job_id + 1) * max_array_size - 1))

    # Ensure the end index does not exceed total combinations
    if [ $end_index -ge $total_combinations ]; then
        end_index=$((total_combinations - 1))
    fi

    echo "Submitting job array for indices: $start_index-$end_index"

    # Example Slurm submission (modify based on your actual job script)
    # sbatch --array=${start_index}-${end_index} my_slurm_script.sh
done