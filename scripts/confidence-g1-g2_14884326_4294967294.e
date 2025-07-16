Traceback (most recent call last):
  File "/scratch/user/uqkprat2/analysis/kp_dadi/scripts/confidence_intervals.py", line 156, in <module>
    main(args.filepath, args.bootpath, args.function, args.model, args.eps, args.opt_params, PTS)
  File "/scratch/user/uqkprat2/analysis/kp_dadi/scripts/confidence_intervals.py", line 113, in main
    for i, label in enumerate(p_labels + ["theta"]):
TypeError: can only concatenate str (not "list") to str
Traceback (most recent call last):
  File "/scratch/user/uqkprat2/analysis/kp_dadi/scripts/confidence_intervals.py", line 156, in <module>
    main(args.filepath, args.bootpath, args.function, args.model, args.eps, args.opt_params, PTS)
  File "/scratch/user/uqkprat2/analysis/kp_dadi/scripts/confidence_intervals.py", line 113, in main
    for i, label in enumerate(p_labels + "theta"):
TypeError: can only concatenate str (not "list") to str
Traceback (most recent call last):
  File "/scratch/user/uqkprat2/analysis/kp_dadi/scripts/confidence_intervals.py", line 156, in <module>
    # Define optimisation bounds
  File "/scratch/user/uqkprat2/analysis/kp_dadi/scripts/confidence_intervals.py", line 114, in main
    print(f"Results written to {out_name}")
IndexError: list index out of range
Traceback (most recent call last):
  File "/scratch/user/uqkprat2/analysis/kp_dadi/scripts/confidence_intervals.py", line 160, in <module>
    main(args.filepath, args.bootpath, args.function, args.model, args.eps, args.opt_params, PTS)
  File "/scratch/user/uqkprat2/analysis/kp_dadi/scripts/confidence_intervals.py", line 113, in main
    out.write(f"{label}\t{opt[i]}\t{low[i]}\t{upp[i]}\n")
IndexError: list index out of range
