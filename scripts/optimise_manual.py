#!/usr/anaconda3/env/dadi211/bin/python
# -*- coding: utf-8 -*-

"""
@author: Katharine Prata
@date created: 5/5/21
@description: Optimising sample fit to demographic models with real time optimisation monitoring.

Optimises using Nelder-mead simplex (can change log_fmin within script if you require a different optimiser).

Inputs:
File: The fs with naming which includes the format of the fs, i,e., subsampled/projected/neither
Arguments:
fs = Pop1-Pop2
model = iso_inbred (refers to custom model module: demo_models_kp.py and used alias for model see script below)
masked = yes or no
method = subsample, projection or none
folds = any integer (best to use: 1, 2 or 3)

Run from script path directory and make sure your custom demographic module is within the directory.

Outputs: dadi_optimisation.txt
Contains the optimal parameters of a run of any model x pop combination

Compatible with python 3.6.11 and dadi 2.1.1
"""

import dadi
import demo_models_kp
import argparse
import numpy
import SETTINGS


def main(fs, model, masked, folds, int_params, PTS, method=None, path=None):
    # Import and define data constants
    if method == "subsample":
        data = dadi.Spectrum.from_file('../data/fs/{}_subsampled.fs'.format(fs))
    elif method == "projection":
        data = dadi.Spectrum.from_file('../data/fs/{}_projected.fs'.format(fs))
    else:
        data = dadi.Spectrum.from_file('../data/fs/{}.fs'.format(fs))
    pops = "{}".format(fs)
    pop_ids = pops.split("-")

    # Mask singletons and doubletons
    if masked == "low":
        data.mask[1, 0] = True
        data.mask[0, 1] = True
        data.mask[2, 0] = True
        data.mask[0, 2] = True
        data.mask[1, 1] = True
        print("Low frequencies masked")
    elif masked == "mid" and len(data.sample_sizes) == 1:
        mid = data.sample_sizes / 2
        mid = int(mid[0])
        data.mask[mid] = True
        print("Mid frequencies masked")
    else:
        print("No masking")

    # Print useful information about the sfs
    print("The datafile is {}".format(fs))
    print("\n* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n")
    print("Data for site frequency spectrum:")
    print("Sample sizes: {}".format(data.sample_sizes))
    print("Sum of SFS: {}".format(numpy.around(data.S(), 2)))
    if len(pop_ids) == 2:
        print("FST of SFS: {}".format(numpy.around(data.Fst(), 2)))
    if len(pop_ids) == 1:
        print("Tajima's D: {}".format(numpy.around(data.Tajima_D(), 2)))
    print("\n* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n")

    # Need to manually alter the upper and lower parameter limits in settings files,
    # model functions are defined in demo_models_kp.py
    # Use nicknames for models, e.g., "snm" instead of model function name, e.g., "no_divergence"
    num, p_labels, upper, lower, model_fun = SETTINGS.get_settings(model)

    # Create log file.
    out_name = path + "dadi_optimisation.txt"
    with open(out_name, 'a') as opt_out:
        if opt_out.tell() == 0:
            print('Creating a new file\n')
            opt_out.write(
                "Pop\tMask\tModel\tFolds\tlog-likelihood\tAIC\tchi-squared\ttheta\tinitial_params\toptimised_params\t"
                "optimised_params_labels\n")
        else:
            print('File exists, appending\n')

    # This is our initial guess for the parameters, which is somewhat arbitrary.
    if int_params is None:
        p0 = [1] * num
    else:
        p0 = int_params
    # Optional paste below your optimised params to start from a specified place.
    # For example:
    # if model == "iso_inbred":
    #    if fs == "AG1-AG2":
    #        p1 = [2.1498, 97.3976, 0.0374, 0.0231, 0.4141]
    #    elif fs == "AL1-AL2":
    #        p1 = [7.4065, 1.9955, 0.0302, 0.0148, 0.4264]
    #    elif fs == "AG1-AL1":
    #        p1 = [0.8382, 41.4364, 0.0371, 0.0229, 0.8367]
    #    elif fs == "AG1-AL2":
    #        p1 = [1.0488, 46.8787, 0.0234, 0.00001, 0.8585]
    #    elif fs == "AG2-AL1":
    #        p1 = [0.7615, 55.4268, 0.0261, 0.0712, 0.9558]
    #    elif fs == "AG2-AL2":
    #        p1 = [0.6783, 44.3377, 0.055, 0.0398, 0.7721]
    #    else:
    #        p1 = p0
    # if using custom starting parameter values then comment out the line below
    p1 = p0

    # Comments from software example (Gutenkunst et al., 2009)
    # Make the extrapolating version of our demographic model function.
    func_ex = dadi.Numerics.make_extrap_log_func(model_fun)

    # Perturb our parameters before optimisation. This does so by taking each
    # parameter a up to a factor of "folds" up or down.
    p1 = dadi.Misc.perturb_params(p1, fold=folds, upper_bound=upper, lower_bound=lower)
    # Do the optimization. By default we assume that theta is a free parameter,
    # since it's trivial to find given the other parameters. If you want to fix
    # theta, add a multinom=False to the call.
    # The maxiter argument restricts how long the optimiser will run. For real
    # runs, you will want to set this value higher (at least 10), to encourage
    # better convergence. You will also want to run optimization several times
    # using multiple sets of initial parameters, to be confident you've actually
    # found the true maximum likelihood parameters.

    print('\nInitial parameters are {}\n'.format(numpy.around(p1, 2)))
    print('* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *')
    print('* * * * * * * * * * * * * * * * * * Beginning  optimisation * * * * * * * * * * * * * * * * * *')
    print('* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *')
    param_opt = dadi.Inference.optimize_log_fmin(p1, data, func_ex, PTS,
                                                 lower_bound=lower,
                                                 upper_bound=upper,
                                                 verbose=1, maxiter=100)
    # The verbose argument controls how often progress of the optimizer should be
    # printed. It's useful to keep track of optimisation process.

    # Calculate sim model using parameters optimised (p_opt).
    sim_model = func_ex(param_opt, data.sample_sizes, PTS)

    # Calculate theta
    theta = dadi.Inference.optimal_sfs_scaling(sim_model, data)
    theta = numpy.around(theta, 4)

    # Calculate likelihood
    ll = dadi.Inference.ll_multinom(sim_model, data)
    ll = numpy.around(ll, 4)

    # Calculate AIC
    aic = (-2 * (float(ll))) + (2 * len(param_opt))

    # Calculate  Chi^2
    scaled_sim_model = sim_model * theta
    folded_sim_model = scaled_sim_model.fold()
    if data.folded == folded_sim_model.folded:
        chi2 = numpy.sum((folded_sim_model - data) ** 2 / folded_sim_model)
        chi2 = numpy.around(chi2, 4)
    else:
        chi2 = numpy.sum((scaled_sim_model - data) ** 2 / scaled_sim_model)
        chi2 = numpy.around(chi2, 4)

    # Store results with likelihoods, theta and parameters (max. likelihood)
    results = [ll, aic, chi2, theta, param_opt]

    # Export the results to file
    with open(out_name, "a") as opt_out:
        easy_p = ",".join([str(numpy.around(x, 4)) for x in results[4]])
        int_p = ",".join([str(numpy.around(x, 4)) for x in p0])
        opt_out.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\n".format(fs, masked, model, folds[0],
                                                                                        results[0], results[1],
                                                                                        results[2],
                                                                                        results[3], int_p, easy_p,
                                                                                        p_labels))
    print('* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *')
    print('* * * * * * * * * * * * * * * * * *  Finished optimisation  * * * * * * * * * * * * * * * * * *')
    print('* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *')


if __name__ == '__main__':
    # Arguments
    parser = argparse.ArgumentParser(
        prog="dadi Optimisation",
        description="A script for optimising dadi models. Specify datafile, model, and other parameters.",
        usage="%(prog)s [options] <fs> <model> <masked> <method> <folds> <out_path>"
    )

    # Required positional arguments
    parser.add_argument("fs", help="Name of the fs file.")
    parser.add_argument("model", help="Name of the model to optimise.")
    parser.add_argument("masked", help="Masking method (e.g., 'mid').")
    parser.add_argument("method", help="Whether fs is projected, subsampled or not (e.g., 'projection')")
    parser.add_argument("folds", type=int, help="Number of folds for cross-validation.", nargs=1)
    parser.add_argument("out_path", help="Output path for results.")

    # Optional arguments
    parser.add_argument(
        "-p", "--int_params",
        nargs="+",
        type=float,
        help="Initial parameters for optimisation (space-separated list)."
    )

    args = parser.parse_args()

    # Setting variables
    fs = args.fs
    model = args.model
    masked = args.masked
    method = args.method
    folds = args.folds
    int_params = args.int_params

    # Need to manually define in SETTINGS.py
    # Define optimisation bounds
    PTS = SETTINGS.SET_PTS
    print("PTS is {}".format(PTS))

    # If you are wanting to export data to a specific location,
    # uncomment the proceeding comment and argument parse,
    # then add path variable to main function.
    path = "{}".format(args.out_path)

    main(fs, model, masked, folds, int_params, PTS, method=method, path=path)
