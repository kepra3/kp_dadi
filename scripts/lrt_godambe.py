#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import numpy as np
from dadi import Spectrum, Numerics, Inference, Godambe
import plot_fs
import SETTINGS
import sys
import glob, os
from scipy.stats import chi2

def main(data, full_model, nested_model, bootpath, mask, opt_full, opt_nested, nested_indices, PTS):
    # Load data
    fs = Spectrum.from_file(f"../data/fs/{data}.fs")
    
    # Apply masking
    plot_fs.apply_mask(fs, mask)


    # Prepare models
    model_fun_full = SETTINGS.get_settings(full_model)
    model_fun_nested = SETTINGS.get_settings(nested_model)
    func_ex_full = Numerics.make_extrap_log_func(model_fun_full)
    func_ex_nested = Numerics.make_extrap_log_func(model_fun_nested)


    # Simulate models
    sim_full = func_ex_full(opt_full, fs.sample_sizes, PTS)
    sim_nested = func_ex_nested(opt_nested, fs.sample_sizes, PTS)

    
    # Log-likelihoods
    ll_full = Inference.ll_multinom(sim_full, fs)
    ll_nested = Inference.ll_multinom(sim_nested, fs)
    print(f"Full model log-likelihood: {ll_full}")
    print(f"Nested model log-likelihood: {ll_nested}")

    # Load bootstraps
    bootstrap_files = sorted(glob.glob(os.path.join(bootpath, "*.fs")))
    all_boot = [Spectrum.from_file(f) for f in bootstrap_files]
    print(f"Loaded {len(all_boot)} bootstrap spectra from {bootpath}")

    # Calculate eps from the full model optimised parameters
    #eps_array = np.array([abs(p) / 1000 if p != 0 else 1e-6 for p in opt_full])
    eps = 0.001

    # Godambe adjustment
    adj = Godambe.LRT_adjust(func_ex_full, PTS, all_boot, opt_full, fs, nested_indices, multinom=True, eps=eps)
    print(f"Godambe adjustment: {adj}")

    # LRT statistic
    D = adj * 2 * (ll_full - ll_nested)
    print(f"Adjusted LRT statistic: {D}")

    # p-value from chi-squared distribution with 1 degree of freedom
    p = 1 - chi2.cdf(D, df=1)
    print(f"LRT p-value (chi^2, 1 df): {p}")

    # Weighted chi-squared p-value
    # If a single parameter is on the. boundary of parameter space the null distribution should be adjusted
    #p = Godambe.sum_chi2_ppf(D, weights)
    #print(f"LRT p-value: {p}")

    return D, p

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Log-likelihood ratio test with Godambe adjustment using dadi.")
    parser.add_argument("data", help="Name of the fs file (without .fs).")
    parser.add_argument("full_model", help="Full model name.")
    parser.add_argument("nested_model", help="Nested model name.")
    parser.add_argument("bootpath", help="Path to bootstrap spectra files.")
    parser.add_argument("mask", help="Mask type ('low', 'mid', 'both', or 'none').")
    #parser.add_argument("fold", choices=["folded", "unfolded"], help="Fold type.")
    parser.add_argument("--opt_full", nargs='+', type=float, required=True, help="Optimised parameters for full model.")
    parser.add_argument("--opt_nested", nargs='+', type=float, required=True, help="Optimised parameters for nested model.")
    parser.add_argument("--nested_indices", nargs='+', type=int, required=True, help="Indices of nested parameters in full model.")
    #parser.add_argument("--weights", nargs='+', type=float, required=True, help="Weights for chi-squared mixture (e.g., 0.5 0.5).")
    args = parser.parse_args()

    PTS = SETTINGS.SET_PTS

    main(
        args.data,
        args.full_model,
        args.nested_model,
        args.bootpath,
        args.mask,
        args.opt_full,
        args.opt_nested,
        args.nested_indices,
        PTS,
    )