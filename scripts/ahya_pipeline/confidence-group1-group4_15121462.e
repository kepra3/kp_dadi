WARNING:Inference:Model is masked in some entries where data is not.
WARNING:Inference:Number of affected entries is 280. Sum of data in those entries is 318486:
WARNING:Inference:Model is masked in some entries where data is not.
WARNING:Inference:Number of affected entries is 257. Sum of data in those entries is 151347:
Traceback (most recent call last):
  File "/scratch/user/uqkprat2/analysis/kp_dadi/scripts/confidence_intervals.py", line 164, in <module>
    main(args.filepath, args.bootpath, args.function, args.model, args.eps, args.opt_params, PTS)
  File "/scratch/user/uqkprat2/analysis/kp_dadi/scripts/confidence_intervals.py", line 73, in main
    param_confidence_intervals = Godambe.GIM_uncert(func_ex, PTS, all_boot, opt, fs,
  File "/home/uqkprat2/.conda/envs/gadma_env/lib/python3.9/site-packages/dadi/Godambe.py", line 289, in GIM_uncert
    GIM, H, J, cU = get_godambe(func_ex, grid_pts, all_boot, p0, data, eps, log,
  File "/home/uqkprat2/.conda/envs/gadma_env/lib/python3.9/site-packages/dadi/Godambe.py", line 246, in get_godambe
    J_inv = numpy.linalg.inv(J)
  File "/home/uqkprat2/.conda/envs/gadma_env/lib/python3.9/site-packages/numpy/linalg/linalg.py", line 561, in inv
    ainv = _umath_linalg.inv(a, signature=signature, extobj=extobj)
  File "/home/uqkprat2/.conda/envs/gadma_env/lib/python3.9/site-packages/numpy/linalg/linalg.py", line 112, in _raise_linalgerror_singular
    raise LinAlgError("Singular matrix")
numpy.linalg.LinAlgError: Singular matrix
Traceback (most recent call last):
  File "/scratch/user/uqkprat2/analysis/kp_dadi/scripts/confidence_intervals.py", line 164, in <module>
    main(args.filepath, args.bootpath, args.function, args.model, args.eps, args.opt_params, PTS)
  File "/scratch/user/uqkprat2/analysis/kp_dadi/scripts/confidence_intervals.py", line 73, in main
    param_confidence_intervals = Godambe.GIM_uncert(func_ex, PTS, all_boot, opt, fs,
  File "/home/uqkprat2/.conda/envs/gadma_env/lib/python3.9/site-packages/dadi/Godambe.py", line 289, in GIM_uncert
    GIM, H, J, cU = get_godambe(func_ex, grid_pts, all_boot, p0, data, eps, log,
  File "/home/uqkprat2/.conda/envs/gadma_env/lib/python3.9/site-packages/dadi/Godambe.py", line 246, in get_godambe
    J_inv = numpy.linalg.inv(J)
  File "/home/uqkprat2/.conda/envs/gadma_env/lib/python3.9/site-packages/numpy/linalg/linalg.py", line 561, in inv
    ainv = _umath_linalg.inv(a, signature=signature, extobj=extobj)
  File "/home/uqkprat2/.conda/envs/gadma_env/lib/python3.9/site-packages/numpy/linalg/linalg.py", line 112, in _raise_linalgerror_singular
    raise LinAlgError("Singular matrix")
numpy.linalg.LinAlgError: Singular matrix
Traceback (most recent call last):
  File "/scratch/user/uqkprat2/analysis/kp_dadi/scripts/confidence_intervals.py", line 164, in <module>
    main(args.filepath, args.bootpath, args.function, args.model, args.eps, args.opt_params, PTS)
  File "/scratch/user/uqkprat2/analysis/kp_dadi/scripts/confidence_intervals.py", line 73, in main
    param_confidence_intervals = Godambe.GIM_uncert(func_ex, PTS, all_boot, opt, fs,
  File "/home/uqkprat2/.conda/envs/gadma_env/lib/python3.9/site-packages/dadi/Godambe.py", line 289, in GIM_uncert
    GIM, H, J, cU = get_godambe(func_ex, grid_pts, all_boot, p0, data, eps, log,
  File "/home/uqkprat2/.conda/envs/gadma_env/lib/python3.9/site-packages/dadi/Godambe.py", line 246, in get_godambe
    J_inv = numpy.linalg.inv(J)
  File "/home/uqkprat2/.conda/envs/gadma_env/lib/python3.9/site-packages/numpy/linalg/linalg.py", line 561, in inv
    ainv = _umath_linalg.inv(a, signature=signature, extobj=extobj)
  File "/home/uqkprat2/.conda/envs/gadma_env/lib/python3.9/site-packages/numpy/linalg/linalg.py", line 112, in _raise_linalgerror_singular
    raise LinAlgError("Singular matrix")
numpy.linalg.LinAlgError: Singular matrix
Traceback (most recent call last):
  File "/scratch/user/uqkprat2/analysis/kp_dadi/scripts/confidence_intervals.py", line 164, in <module>
    main(args.filepath, args.bootpath, args.function, args.model, args.eps, args.opt_params, PTS)
  File "/scratch/user/uqkprat2/analysis/kp_dadi/scripts/confidence_intervals.py", line 73, in main
    param_confidence_intervals = Godambe.GIM_uncert(func_ex, PTS, all_boot, opt, fs,
  File "/home/uqkprat2/.conda/envs/gadma_env/lib/python3.9/site-packages/dadi/Godambe.py", line 289, in GIM_uncert
    GIM, H, J, cU = get_godambe(func_ex, grid_pts, all_boot, p0, data, eps, log,
  File "/home/uqkprat2/.conda/envs/gadma_env/lib/python3.9/site-packages/dadi/Godambe.py", line 246, in get_godambe
    J_inv = numpy.linalg.inv(J)
  File "/home/uqkprat2/.conda/envs/gadma_env/lib/python3.9/site-packages/numpy/linalg/linalg.py", line 561, in inv
    ainv = _umath_linalg.inv(a, signature=signature, extobj=extobj)
  File "/home/uqkprat2/.conda/envs/gadma_env/lib/python3.9/site-packages/numpy/linalg/linalg.py", line 112, in _raise_linalgerror_singular
    raise LinAlgError("Singular matrix")
numpy.linalg.LinAlgError: Singular matrix
Traceback (most recent call last):
  File "/scratch/user/uqkprat2/analysis/kp_dadi/scripts/confidence_intervals.py", line 164, in <module>
    main(args.filepath, args.bootpath, args.function, args.model, args.eps, args.opt_params, PTS)
  File "/scratch/user/uqkprat2/analysis/kp_dadi/scripts/confidence_intervals.py", line 73, in main
    param_confidence_intervals = Godambe.GIM_uncert(func_ex, PTS, all_boot, opt, fs,
  File "/home/uqkprat2/.conda/envs/gadma_env/lib/python3.9/site-packages/dadi/Godambe.py", line 289, in GIM_uncert
    GIM, H, J, cU = get_godambe(func_ex, grid_pts, all_boot, p0, data, eps, log,
  File "/home/uqkprat2/.conda/envs/gadma_env/lib/python3.9/site-packages/dadi/Godambe.py", line 246, in get_godambe
    J_inv = numpy.linalg.inv(J)
  File "/home/uqkprat2/.conda/envs/gadma_env/lib/python3.9/site-packages/numpy/linalg/linalg.py", line 561, in inv
    ainv = _umath_linalg.inv(a, signature=signature, extobj=extobj)
  File "/home/uqkprat2/.conda/envs/gadma_env/lib/python3.9/site-packages/numpy/linalg/linalg.py", line 112, in _raise_linalgerror_singular
    raise LinAlgError("Singular matrix")
numpy.linalg.LinAlgError: Singular matrix
Traceback (most recent call last):
  File "/scratch/user/uqkprat2/analysis/kp_dadi/scripts/confidence_intervals.py", line 164, in <module>
    main(args.filepath, args.bootpath, args.function, args.model, args.eps, args.opt_params, PTS)
  File "/scratch/user/uqkprat2/analysis/kp_dadi/scripts/confidence_intervals.py", line 73, in main
    param_confidence_intervals = Godambe.GIM_uncert(func_ex, PTS, all_boot, opt, fs,
  File "/home/uqkprat2/.conda/envs/gadma_env/lib/python3.9/site-packages/dadi/Godambe.py", line 289, in GIM_uncert
    GIM, H, J, cU = get_godambe(func_ex, grid_pts, all_boot, p0, data, eps, log,
  File "/home/uqkprat2/.conda/envs/gadma_env/lib/python3.9/site-packages/dadi/Godambe.py", line 246, in get_godambe
    J_inv = numpy.linalg.inv(J)
  File "/home/uqkprat2/.conda/envs/gadma_env/lib/python3.9/site-packages/numpy/linalg/linalg.py", line 561, in inv
    ainv = _umath_linalg.inv(a, signature=signature, extobj=extobj)
  File "/home/uqkprat2/.conda/envs/gadma_env/lib/python3.9/site-packages/numpy/linalg/linalg.py", line 112, in _raise_linalgerror_singular
    raise LinAlgError("Singular matrix")
numpy.linalg.LinAlgError: Singular matrix
