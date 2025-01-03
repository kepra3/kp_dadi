#!/usr/anaconda3/env/dadi211/bin/python
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


def main(function, snps, method, model, sims, eps, opt, PTS):
    """
    eps: Fractional stepsize to use when taking finite-difference derivatives.
        Note that if eps*param is < 1e-6, then the step size for that parameter
        will simply be eps, to avoid numerical issues with small parameter
        perturbations.
    """
    # Import spectrum
    if method == "subsample":
        fs_path = "../data/fs/{}_subsampled.fs".format(snps)
    elif method == "projection":
        fs_path = "../data/fs/{}_projected.fs".format(snps)
    else:
        fs_path = "../data/fs/{}.fs".format(snps)
    fs = Spectrum.from_file(fs_path)
    ns = fs.sample_sizes

    # The parameter confidence interval txt file which will be output-1-1timeperiod-lowmask
    out_name = "../results/{}_{}_confidence_intervals.txt".format(snps, model)

    # Mask singletons and doubletons
    fs.mask[1, 0] = True
    fs.mask[0, 1] = True
    fs.mask[2, 0] = True
    fs.mask[0, 2] = True
    fs.mask[1, 1] = True

    # Print information about the spectrum
    print("\n* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n")
    print("Data for site frequency spectrum:\n")
    print("Sample sizes: {}".format(fs.sample_sizes))
    print("Sum of SFS: {}".format(np.around(fs.S(), 2)))
    print("FST of SFS: {}".format(np.around(fs.Fst(), 2)))
    print("\n* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n")

    # Models of choice
    if model == "iso_inbred":
        model_fun = demo_models_kp.iso_inbreeding
    elif model == "mig_inbred":
        model_fun = demo_models_kp.mig_inbreeding
    elif model == "anc_mig":
        model_fun = demo_models_kp.anc_sym_mig_inbred
    elif model == "sec_cont":
        model_fun = demo_models_kp.sec_contact_sym_mig_inbred
    else:
        model_fun = False
        print("Choose a correct model")

    # Model function
    func_ex = Numerics.make_extrap_log_func(model_fun)
    sim_model = func_ex(opt, ns, PTS)
    ll_model = Inference.ll_multinom(sim_model, fs)

    # The log likelihood of the model
    print('Maximum log composite likelihood: {0}'.format(ll_model))

    # The optimal value of theta given the model
    theta = Inference.optimal_sfs_scaling(sim_model, fs)
    print('Optimal value of theta: {0}'.format(theta))

    # Import all bootstraps
    all_boot = [Spectrum.from_file('../results/bootstraps/{}_bootstrap_vcf_{}.fs'.format(snps, i))
                for i in range(0, sims)]

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

    # If statements are to adjust the number of columns with the varying number of parameters within each of my models
    if model == "iso_inbred":
        with open(out_name, 'a') as file_out:
            file_out.write('Pops\tModel\ttheta\tnu1\tnu2\tF1\tF2\tT\teps\n')
        with open(out_name, 'a') as file_out:
            file_out.write(f'{snps}\t{model}\t{low[-1]}\t{low[0]}\t{low[1]}\t{low[2]}\t{low[3]}\t{low[4]}\t{eps}\n')
        with open(out_name, 'a') as file_out:
            file_out.write(f'{snps}\t{model}\t{upp[-1]}\t{upp[0]}\t{upp[1]}\t{upp[2]}\t{upp[3]}\t{upp[4]}\t{eps}\n')
    elif model == "mig_inbred":
        with open(out_name, 'a') as file_out:
            file_out.write('Pops\tModel\ttheta\tnu1\tnu2\tF1\tF2\tm\tT\teps\n')
        with open(out_name, 'a') as fh_out:
            fh_out.write(
                f'{snps}\t{model}\t{low[-1]}\t{low[0]}\t{low[1]}\t{low[2]}\t{low[3]}\t{low[4]}\t{low[5]}\t{eps}\n')
        with open(out_name, 'a') as fh_out:
            fh_out.write(
                f'{snps}\t{model}\t{upp[-1]}\t{upp[0]}\t{upp[1]}\t{upp[2]}\t{upp[3]}\t{upp[4]}\t{upp[5]}\t{eps}\n')
    else:
        with open('a', out_name) as file_out:
            file_out.write('Pops\tModel\ttheta\tnu1\tnu2\tF1\tF2\tm\tT1\tT2\teps\n')
        with open(out_name, 'a') as fh_out:
            fh_out.write(
                f'{snps}\t{model}\t{low[-1]}\t{low[0]}\t{low[1]}\t{low[2]}\t{low[3]}\t{low[4]}\t{low[5]}'
                f'\t{low[6]}\t{eps}\n')
        with open(out_name, 'a') as file_out:
            file_out.write(f'{snps}\t{model}\t{upp[-1]}\t{upp[0]}\t{upp[1]}\t{upp[2]}\t{upp[3]}\t{upp[4]}\t{upp[5]}'
                           f'\t{upp[6]}\t{eps}\n')


if __name__ == "__main__":
    # Arguments
    parser = argparse.ArgumentParser(prog="Parameter uncertainty")
    parser.add_argument("function")
    parser.add_argument("snps")
    parser.add_argument("method")
    parser.add_argument("model")
    parser.add_argument("sims", type=int)
    parser.add_argument("eps", type=float)
    parser.add_argument("-o", "--opt_params", nargs="+", type=float)
    args: Namespace = parser.parse_args()

    # Setting variables
    function = args.snps
    snps = args.snps
    method = args.method
    model = args.model
    sims = args.sims
    eps = args.eps
    opt = args.opt_params

    # Extrapolation of grid size
    PTS = [50, 60, 70]

    main(function, snps, method, model, sims, eps, opt, PTS)
