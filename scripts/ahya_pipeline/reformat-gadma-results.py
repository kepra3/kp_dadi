import os
import re
import argparse

parser = argparse.ArgumentParser(
    prog="Reformat GADMA Results",
    description="A script to reformat GADMA results from multiple directories into a single output file.",
    usage="%(prog)s <base_dir> <output_file>"
)

parser.add_argument("base_dir", help="Name of the base directory containing model results.")
parser.add_argument("output_file", help=" Output file to write the reformatted results.")

args = parser.parse_args()

# Set your base directory
base_dir = args.base_dir
output_file = args.output_file

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
            print(f"Processing file: {file_path}")

            with open(file_path) as f:
                lines = f.readlines()[2:]  # Skip the first two header lines
                print(f"Number of lines after header: {len(lines)}")
            for line in lines:
                run_match = re.search(r"Run\s+(\d+)", line)
                logL_match = re.search(r"Run\s+\d+\s+(-?\d+\.?\d*)", line)

                if not run_match or not logL_match:
                    print("No run or logL match found, skipping line.")
                    continue

                run = f"Run {run_match.group(1)}"
                logL = logL_match.group(1)

                # Extract parameters and their names
                params = param_pattern.findall(line)
                theta_match = theta_pattern.search(line)

                if theta_match:
                    print("No parameters found, skipping line.")
                    params.append((theta_match.group(1), "theta"))

                if not params:
                    continue

                values = [v for v, _ in params]
                names = [n for _, n in params]

                names_str = f"({', '.join(names)})"
                values_str = f"({', '.join(values)})"

                print(f"Writing: {dir_name}\t{run}\t{logL}\t{names_str}\t{values_str}")
                out.write(f"{dir_name}\t{run}\t{logL}\t{names_str}\t{values_str}\n")
