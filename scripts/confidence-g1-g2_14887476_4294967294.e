Traceback (most recent call last):
  File "/scratch/user/uqkprat2/analysis/kp_dadi/scripts/confidence_intervals.py", line 159, in <module>
    main(args.filepath, args.bootpath, args.function, args.model, args.eps, args.opt_params, PTS)
  File "/scratch/user/uqkprat2/analysis/kp_dadi/scripts/confidence_intervals.py", line 117, in main
    out.write(f"{label}\t{opt[i]}\t{low[i]}\t{upp[i]}\n")
IndexError: list index out of range
Traceback (most recent call last):
  File "/scratch/user/uqkprat2/analysis/kp_dadi/scripts/confidence_intervals.py", line 159, in <module>
    main(args.filepath, args.bootpath, args.function, args.model, args.eps, args.opt_params, PTS)
  File "/scratch/user/uqkprat2/analysis/kp_dadi/scripts/confidence_intervals.py", line 112, in main
    out.write(f"{label}\t{opt[i]}\t{low[i]}\t{upp[i]}\n")
IndexError: list index out of range
slurmstepd: error: *** JOB 14887476 ON bun129 CANCELLED AT 2025-07-16T17:15:12 ***
