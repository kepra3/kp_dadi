import os
import re

# Set your base directory
base_dir = "../results/hpc_gadma"
output_file = "../results/gadma_results_combined.txt"

# Match param format like: 1.043(nu_1)
param_pattern = re.compile(r"([\d.eE+-]+)\(([^)]+)\)")
theta_pattern = re.compile(r"theta\s*=\s*([\d.eE+-]+)")

with open(output_file, "w") as out:
    # Write header
    out.write("Directory\tRun\tLogLikelihood\tModelParameters\tModelValues\n")

    for root, dirs, files in os.walk(base_dir):
        if "all_extracted_models.txt" in files:
            dir_name = os.path.basename(root)
            file_path = os.path.join(root, "all_extracted_models.txt")

            with open(file_path) as f:
                lines = f.readlines()[2:]  # Skip the first two header lines

            for line in lines:
                run_match = re.search(r"Run\s+(\d+)", line)
                logL_match = re.search(r"Run\s+\d+\s+(-?\d+\.?\d*)", line)

                if not run_match or not logL_match:
                    continue

                run = f"Run {run_match.group(1)}"
                logL = logL_match.group(1)

                # Extract parameters and their names
                params = param_pattern.findall(line)
                theta_match = theta_pattern.search(line)

                if theta_match:
                    params.append((theta_match.group(1), "theta"))

                if not params:
                    continue

                values = [v for v, _ in params]
                names = [n for _, n in params]

                names_str = f"({', '.join(names)})"
                values_str = f"({', '.join(values)})"

                out.write(f"{dir_name}\t{run}\t{logL}\t{names_str}\t{values_str}\n")
