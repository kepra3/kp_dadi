#!/bin/bash


SFS=$1
MODEL=$2

# Function to generate all combinations of parameters
generate_combinations() {
  local depth=$1
  local prefix=$2
  if [ "$depth" -eq 0 ]; then
    echo "$prefix"
  else
    generate_combinations $((depth - 1)) "$prefix 1"
    generate_combinations $((depth - 1)) "$prefix 10"
  fi
}

# Number of parameters
if [ "${MODEL}" == "no_mig" ]; then
	n=3
elif [ "${MODEL}" == "asym_mig" ]; then
	n=5
elif [ "${MODEL}" == "anc_asym_mig" ]; then
	n=6
elif [ "${MODEL}" == "sec_cont_asym_mig" ]; then
        n=6
elif [ "${MODEL}" == "split_bottle" ]; then
	n=6
elif [ "${MODEL}" == "split_bottle_asym_mig" ]; then
	n=10
elif [ "${MODEL}" == "split_bottle_anc_asym_mig" ]; then
        n=8
elif [ "${MODEL}" == "split_bottle_sec_asym_mig" ]; then
        n=8
elif [ "${MODEL}" == "split_sizechange" ]; then
        n=6
elif [ "${MODEL}" == "split_sizechange_asym_mig" ]; then
        n=10
elif [ "${MODEL}" == "split_sizechange_anc_asym_mig" ]; then
        n=8
elif [ "${MODEL}" == "split_sizechange_sec_asym_mig" ]; then
        n=8
else
	echo "Unknown model: $MODEL"
	exit 1
fi
echo "Number of parameters: $n"

# Generate all parameter combinations
combinations=$(generate_combinations "$n" "")

# Debug: Print all combinations
echo "Generated combinations:"
echo "$combinations"

# Run the Python script for each combination
while read -r combination; do
  # Trim any leading/trailing whitespace
  trimmed_combination=$(echo "$combination" | xargs)

  echo "Running with parameters: $trimmed_combination"
  python optimise_manual.py $SFS $MODEL no no 3 ../results/opt1 -p $trimmed_combination
done <<< "$combinations"