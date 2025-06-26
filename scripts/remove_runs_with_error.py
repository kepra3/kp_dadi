import argparse
import os

# Argument parser for command line arguments
parser = argparse.ArgumentParser(description="Filter output based on error identifiers.")
parser.add_argument("error_file", help="Path to the .e file containing error identifiers.")
parser.add_argument("output_file", help="Path to the .txt file to filter based on error identifiers.")
parser.add_argument("filtered_output_file", help="Path to the output file for filtered results.")
args = parser.parse_args() 

# Check if the provided files exist
if not os.path.exists(args.error_file):
    raise FileNotFoundError(f"Error file '{args.error_file}' does not exist.")
if not os.path.exists(args.output_file):
    raise FileNotFoundError(f"Output file '{args.output_file}' does not exist.")

# Read identifiers from .e file
bad_ids = set()
with open(args.error_file) as ef:
    for line in ef:
        if "WARNING" in line:
            next_line = next(ef, "")
            #print(f"DEBUG: Next line after warning: {next_line.strip()}")
            if next_line.startswith("Writing:"):
                # Extract unique identifier, e.g., directory and run
                parts = next_line.strip().split("\t")
                #print(f"DEBUG: Split parts: {parts}")
                if len(parts) > 2:
                    dir_part = parts[0].replace("Writing: ", "").strip()
                    run_part = parts[1].strip()
                    #print(f"DEBUG: Adding bad_id: ({dir_part}, {run_part})")
                    bad_ids.add((dir_part, run_part))

print(f"Identifiers to remove (total {len(bad_ids)}): {list(bad_ids)[:]} ...")

# Filter .txt file
with open(args.output_file) as inf, open(args.filtered_output_file, "w") as outf:
    header = inf.readline()
    outf.write(header)
    for line in inf:
        parts = line.strip().split("\t")
        if len(parts) > 2:
            key = (parts[0].strip(), parts[1].strip())
            if key in bad_ids:
                #print(f"DEBUG: Skipping line due to bad_id: {key}")
            else:
                outf.write(line)