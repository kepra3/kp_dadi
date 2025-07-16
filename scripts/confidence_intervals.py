#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
@author: Katharine Prata
@date created: 2/6/21
@description: obtains GIM or FIM uncertainty from optimised parameters and bootstraps.

Script modified from YRI_CEU.py if use Godambe have to cite Coffman et al. (2016).

Compatible with python 3.6.11 and dadi 2.1.1
"""

from argparse import Namespace
from dadi import Numerics, Spectrum, Inference, Godambe
import demo_models_kp
import numpy as np
import argparse
import plot_fs
import SETTINGS
import os
from glob import glob


def main(filepath, bootpath, function, model, eps, opt, PTS):
    """
    eps: Fractional stepsize to use when taking finite-difference derivatives.
        Note that if eps*param is < 1e-6, then the step size for that parameter
        will simply be eps, to avoid numerical issues with small parameter
        perturbations.
    """
    # Extract SNPs name from the file path
    snps = os.path.splitext(os.path.basename(filepath))[0]

    # Import spectrum
    fs = Spectrum.from_file(filepath)
    ns = fs.sample_sizes
    print(f"Loaded spectrum with sample sizes: {ns}")

    # The parameter confidence interval txt file
    out_name = "../results/{}_{}_confidence_intervals.txt".format(snps, model)

    # Calculate statistics
    statistic_name, statistic = plot_fs.calculate_statistic(fs)
    print(statistic_name, statistic)

    # Need to manually alter the upper and lower parameter limits in settings files,
    # model functions are defined in demo_models_kp.py
    # Use nicknames for models, e.g., "snm" instead of model function name, e.g., "no_divergence"
    model_fun, num, p_labels, upper, lower = SETTINGS.get_settings(model, ALL=True)

    # Model function
    func_ex = Numerics.make_extrap_log_func(model_fun)
    sim_model = func_ex(opt, ns, PTS)
    ll_model = Inference.ll_multinom(sim_model, fs)

    # The log likelihood of the model
    print('Maximum log composite likelihood: {0}'.format(ll_model))

    # The optimal value of theta given the model
    theta = Inference.optimal_sfs_scaling(sim_model, fs)
    print('Optimal value of theta: {0}'.format(theta))

    # Import all bootstraps from the specified directory
    bootstrap_files = sorted(glob(os.path.join(bootpath, "*.fs")))
    all_boot = [Spectrum.from_file(f) for f in bootstrap_files]
    print(f"Loaded {len(all_boot)} bootstrap spectra from {bootpath}")

    # Godambe uncertainties
    # param_confidence_intervals contains the estimated standard deviations of each parameter,
    # with theta as the final entry in the list.
    if function == "GIM":
        param_confidence_intervals = Godambe.GIM_uncert(func_ex, PTS, all_boot, opt, fs,
                                                        multinom=True, eps=eps, log=True)
    elif function == "FIM":
        # if want to do FIM
        # For comparison, we can estimate uncertainties with the Fisher Information
        # Matrix, which doesn't account for linkage in the data and thus underestimates
        # uncertainty. (Although it's a fine approach if you think your data is truly
        # unlinked.)
        param_confidence_intervals = Godambe.FIM_uncert(func_ex, PTS, opt, fs,
                                                        multinom=True, eps=eps, log=True)
    else:
        print("Choose uncertainty function")

    print('Estimated parameter standard deviations from {0}: {1}'.format(function, param_confidence_intervals))
    # Add optimal theta to parameter list for comparison
    opt.append(theta)
    print('Optimised parameters: {0}'.format(opt))

    # The lower bound confidence intervals
    low = np.subtract(opt, param_confidence_intervals)
    print('Estimated parameter lower from {0}: {1}'.format(function, low))
    for i in range(len(low)):
        if low[i] < 0:
            low[i] = 0
        else:
            print("not negative")
    low = np.around(low, 4)
    print('Adjusted estimated parameter lower from {0}: {1}'.format(function, low))

    # The upper bound confidence intervals
    upp = opt + param_confidence_intervals
    upp = np.around(upp, 4)
    print('Estimated parameter upper from {0}: {1}'.format(function, upp))

    p_labels = [x.strip() for x in p_labels.split(",")]
    
    # Write results to file
    if os.path.exists(out_name):
        with open(out_name, "a") as out:
            for i, label in enumerate(p_labels + ["theta"]):
                out.write(f"{label}\t{opt[i]}\t{low[i]}\t{upp[i]}\t{eps}\n")
        print(f"Results written to {out_name}")
    else: 
        with open(out_name, "w") as out:
            out.write("Parameter\tOptimised\tLower_CI\tUpper_CI\teps\n")
            for i, label in enumerate(p_labels + ["theta"]):
                out.write(f"{label}\t{opt[i]}\t{low[i]}\t{upp[i]}\t{eps}\n")
        print(f"Results written to {out_name}")
        

if __name__ == "__main__":
    # Arguments
    parser = argparse.ArgumentParser(prog="Parameter uncertainty",
    description="A script for obtaining parameter uncertainty from optimised parameters and bootstraps.",
    usage="%(prog)s [options] <filepath> <bootpath> <function> <model> <sims> <eps> [opt_params]")
    parser.add_argument(
        "filepath",
        help="Path to the data .fs file."
    )
    parser.add_argument(
        "bootpath",
        help="Directory for bootstraps, e.g., '../results/bootstraps/'"
    )
    parser.add_argument(
        "function",
        help="Function to use for uncertainty estimation (e.g., 'GIM' or 'FIM')."
    )
    parser.add_argument(
        "model",
        help="Model to use for the analysis from kp_dadi."
    )
    parser.add_argument(
        "eps",
        help="eps setting (e.g., 0.01, 0.001, 0.0001)",
        type=float
    )
    parser.add_argument(
        "-o", "--opt_params",
        help="optimised paramaters for specific model",
        nargs="+", type=float
        )
    args: Namespace = parser.parse_args()

    # Need to manually define in SETTINGS.py
    # Define optimisation bounds
    PTS = SETTINGS.SET_PTS
    print("PTS is {}".format(PTS))

    main(args.filepath, args.bootpath, args.function, args.model, args.eps, args.opt_params, PTS)
