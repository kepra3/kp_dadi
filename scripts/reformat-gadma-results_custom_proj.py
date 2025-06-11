import os
import re
import argparse
import dadi
import numpy as np
import importlib.util

parser = argparse.ArgumentParser(
    prog="Reformat GADMA Results",
    description="Condense GADMA results and calculate theta for each run.",
    usage="%(prog)s <base_dir> <output_file> <model_py>"
)
parser.add_argument("base_dir", help="Base directory containing model results.")
parser.add_argument("output_file", help="Output file to write the reformatted results.")
parser.add_argument("model_py", help="Path to the custom model Python file.")

args = parser.parse_args()

base_dir = args.base_dir
output_file = args.output_file
model_py = args.model_py

# Load model function dynamically
spec = importlib.util.spec_from_file_location('module', model_py)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
model_func = module.model_func

# Regex for extracting parameters inside parentheses
param_pattern = re.compile(r"([a-zA-Z0-9_]+)=([-\d.eE+]+)")
run_pattern = re.compile(r"Run\s+(\d+)")
logL_pattern = re.compile(r"Run\s+\d+\s+(-?\d+\.?\d*)")

with open(output_file, "w") as out:
    out.write("Directory\tRun\tLogLikelihood\tModelParameters\tModelValues\tTheta\n")
    for root, dirs, files in os.walk(base_dir):
        if "all_extracted_models.txt" in files:
            dir_name = os.path.basename(root)
            file_path = os.path.join(root, "all_extracted_models.txt")
            print(f"Processing file: {file_path}")

            # Extract group name for fs file (matches e.g. group1-group2 anywhere after an underscore)
            group_match = re.search(r'_([^_]+-[^_]+)', dir_name)
            if not group_match:
                print(f"Could not extract group name from directory: {dir_name}")
                continue
            group_name = group_match.group(1)
            fs_file = f"/scratch/user/uqkprat2/analysis/kp_dadi/data/fs/{group_name}_projected0.8.fs"

            # Load fs file for dadi
            try:
                data = dadi.Spectrum.from_file(fs_file)
                ns = data.sample_sizes
                pts = [60, 70, 100]
                func_ex = dadi.Numerics.make_extrap_log_func(model_func)
            except Exception as e:
                print(f"Could not load fs file {fs_file}: {e}")
                continue

            with open(file_path) as f:
                lines = f.readlines()[2:]  # Skip the first two header lines
                print(f"Number of lines after header: {len(lines)}")
            for line in lines:
                run_match = run_pattern.search(line)
                logL_match = logL_pattern.search(line)

                if not run_match or not logL_match:
                    print("No run or logL match found, skipping line.")
                    continue

                run = f"Run {run_match.group(1)}"
                logL = logL_match.group(1)

                # Extract parameter string inside parentheses
                param_str_match = re.search(r"\(([^)]+)\)", line)
                if not param_str_match:
                    print("No parameter string found, skipping line.")
                    continue
                param_str = param_str_match.group(1)
                params = param_pattern.findall(param_str)

                if not params:
                    print("No parameters found, skipping line.")
                    continue

                names = [n for n, v in params]
                values = [v for n, v in params]
                values_float = [float(v) for v in values]

                # Calculate theta using dadi
                try:
                    model = func_ex(values_float, ns, pts)
                    theta = dadi.Inference.optimal_sfs_scaling(model, data)
                except Exception as e:
                    print(f"Error calculating theta for {run}: {e}")
                    theta = "NA"

                names_str = f"({', '.join(names)})"
                values_str = f"({', '.join(values)})"

                print(f"Writing: {dir_name}\t{run}\t{logL}\t{names_str}\t{values_str}\t{theta}")
                out.write(f"{dir_name}\t{run}\t{logL}\t{names_str}\t{values_str}\t{theta}\n")
