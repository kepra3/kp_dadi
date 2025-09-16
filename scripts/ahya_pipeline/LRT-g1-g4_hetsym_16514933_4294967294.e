Traceback (most recent call last):
  File "/scratch/user/uqkprat2/analysis/kp_dadi/scripts/lrt_godambe.py", line 130, in <module>
    main(
  File "/scratch/user/uqkprat2/analysis/kp_dadi/scripts/lrt_godambe.py", line 79, in main
    adj_nested = Godambe.LRT_adjust(func_ex_full, PTS, all_boot, opt_nested_buffered, fs, nested_indices, multinom=True, eps=eps)
  File "/home/uqkprat2/.conda/envs/gadma_env/lib/python3.9/site-packages/dadi/Godambe.py", line 381, in LRT_adjust
    GIM, H, J, cU = get_godambe(diff_func, grid_pts, all_boot, p_nested, data, 
  File "/home/uqkprat2/.conda/envs/gadma_env/lib/python3.9/site-packages/dadi/Godambe.py", line 221, in get_godambe
    hess = -get_hess(func, p0, eps, args=[data])
  File "/home/uqkprat2/.conda/envs/gadma_env/lib/python3.9/site-packages/dadi/Godambe.py", line 123, in get_hess
    hess[ii][jj] = hessian_elem(func, f0, p0, ii, jj, eps, args=args, one_sided=one_sided)
  File "/home/uqkprat2/.conda/envs/gadma_env/lib/python3.9/site-packages/dadi/Godambe.py", line 42, in hessian_elem
    fpp = func(pwork, *args)
  File "/home/uqkprat2/.conda/envs/gadma_env/lib/python3.9/site-packages/dadi/Godambe.py", line 211, in func
    cache[key] = func_ex(params, ns, grid_pts)
  File "/home/uqkprat2/.conda/envs/gadma_env/lib/python3.9/site-packages/dadi/Godambe.py", line 378, in diff_func
    return func_ex(full_params, ns, grid_pts)
  File "/home/uqkprat2/.conda/envs/gadma_env/lib/python3.9/site-packages/dadi/Godambe.py", line 368, in <lambda>
    func_ex = lambda p, ns, pts: p[-1]*func_multi(p[:-1], ns, pts)
  File "/home/uqkprat2/.conda/envs/gadma_env/lib/python3.9/site-packages/dadi/Numerics.py", line 396, in extrap_func
    result_l = list(map(partial_func, pts_l))
  File "/scratch/user/uqkprat2/analysis/kp_dadi/scripts/demo_models_kp.py", line 1513, in custom_model_2p_m0
    phiN2 = Integration.two_pops(phi1, xx, T=t2, nu1=nu21, nu2=nu22, m12=0, m21=0)
  File "/home/uqkprat2/.conda/envs/gadma_env/lib/python3.9/site-packages/dadi/Integration.py", line 319, in two_pops
    return _two_pops_const_params(phi, xx, T, nu1, nu2, m12, m21,
  File "/home/uqkprat2/.conda/envs/gadma_env/lib/python3.9/site-packages/dadi/Integration.py", line 898, in _two_pops_const_params
    raise ValueError('A population size is 0. Has the model been '
ValueError: A population size is 0. Has the model been mis-specified?
