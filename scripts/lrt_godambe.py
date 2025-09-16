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

    # Verify all indices are covered
    print(f"opt_full length: {len(opt_full)}")
    print(f"opt_nested length: {len(opt_nested)}")
    print(f"nested_indices: {nested_indices}")
    print(f"Expected nested length: {len(opt_full) - len(nested_indices)}")
    print(f"opt_full: {opt_full}")
    print(f"opt_nested: {opt_nested}")

    # Always calculate Godambe adjustment with full model
    adj_full = Godambe.LRT_adjust(func_ex_full, PTS, all_boot, opt_full, fs, nested_indices, multinom=True, eps=eps)
    print(f"Godambe adjustment with full model: {adj_full}")

    # Always calculate LRT statistic with full model
    D1 = adj_full * 2 * (ll_full - ll_nested)
    print(f"Adjusted LRT statistic with full model: {D1}")

    # Weighted chi-squared p-value with full model
    weights = (0.5, 0.5)
    p1 = Godambe.sum_chi2_ppf(D1, weights)
    print(f"LRT p-value (chi^2, 1 df) with full model: {round(p1, 5)}")

    # Handle nested model calculations based on model type
    if full_model == "1het_sym":
        print("Skipping nested-based model calculations for 1het_sym model because pop sizes cannot be set to zero in Hessian calculations")
        D = D1
        p = p1
        model_choice = "full"
    else:
        # Create opt_nested_buffered from opt_full structure for other models
        opt_nested_buffered = [1.0] * len(opt_full)  # Start with full model structure

        # Set nested indices to 0 (these are the parameters being tested)
        for idx in nested_indices:
            opt_nested_buffered[idx] = 0

        # Fill remaining positions with opt_nested values
        nested_param_idx = 0
        for i in range(len(opt_full)):
            if i not in nested_indices:
                if nested_param_idx < len(opt_nested):
                    opt_nested_buffered[i] = opt_nested[nested_param_idx]
                    nested_param_idx += 1

        # Check that we used all nested parameters
        if nested_param_idx != len(opt_nested):
            raise ValueError(f"Mismatch: used {nested_param_idx} nested params, but got {len(opt_nested)}")

        print(f"opt_nested_buffered: {opt_nested_buffered}")

        # Calculate nested model Godambe adjustment
        adj_nested = Godambe.LRT_adjust(func_ex_full, PTS, all_boot, opt_nested_buffered, fs, nested_indices, multinom=True, eps=eps)
        print(f"Godambe adjustment with nested model: {adj_nested}")

        # Calculate LRT statistic with nested model
        D2 = adj_nested * 2 * (ll_full - ll_nested)
        print(f"Adjusted LRT statistic with nested model: {D2}")

        # Weighted chi-squared p-value with nested model
        p2 = Godambe.sum_chi2_ppf(D2, weights)
        print(f"LRT p-value (chi^2, 1 df) with nested model: {round(p2, 5)}")

        # Choose most conservative method
        if D1 < D2:
            D = D1
            p = p1
            model_choice = "full"
        else:
            D = D2
            p = p2
            model_choice = "nested"

    print(f"Chosen model parameters for statistic: {model_choice}")
    if p < 0.05:
        print(f"Complex model is preferred.")
    else:
        print(f"Simple model is preferred.")

    return D, p, model_choice

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