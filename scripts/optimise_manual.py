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
snps = Pop1-Pop2
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


def main(snps, model, masked, folds, int_params, PTS, method=None):
    # Import and define data constants
    if method == "subsample":
        data = dadi.Spectrum.from_file('../data/fs/{}_subsampled.fs'.format(snps))
    elif method == "projection":
        data = dadi.Spectrum.from_file('../data/fs/{}_projected.fs'.format(snps))
    else:
        data = dadi.Spectrum.from_file('../data/fs/{}.fs'.format(snps))
    pops = "{}".format(snps)
    pop_ids = pops.split("-")

    # Mask singletons and doubletons
    if masked == "yes":
        data.mask[1, 0] = True
        data.mask[0, 1] = True
        data.mask[2, 0] = True
        data.mask[0, 2] = True
        data.mask[1, 1] = True
    else:
        print("Singletons and Doubletons are not masked")

    # Print useful information about the sfs
    print("The datafile is {}".format(snps))
    print("\n* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n")
    print("Data for site frequency spectrum:")
    print("Sample sizes: {}".format(data.sample_sizes))
    print("Sum of SFS: {}".format(numpy.around(data.S(), 2)))
    if len(pop_ids) == 2:
        print("FST of SFS: {}".format(numpy.around(data.Fst(), 2)))
    print("\n* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n")

    # Need to manually alter the upper and lower parameter limits, model functions are defined in demo_models_kp.py
    # Use nicknames for models, e.g., "snm" instead of model function name, e.g., "no_divergence"
    # Define metadata for models
    if model == "snm.1d":
        # standard neutral model 1d
        num = 1
        p_labels = "nu"
        upper = [200]
        lower = [0.001]
        model_fun = demo_models_kp.no_divergence_1d
    elif model == "size_change":
        num = 2
        p_labels = "nu, T"
        upper = [200, 200]
        lower = [0.001, 0.001]
        model_fun = demo_models_kp.instant_change
    elif model == "bottle":
        num = 3
        p_labels = "nuB, nuF, T"
        upper = [200, 200, 200]
        lower = [0.001, 0.001, 0.001]
        model_fun = demo_models_kp.bottlegrowth
    elif model == "bottle_neck":
        num = 4
        p_labels = "nuB, nuF, TB, TF"
        upper = [200, 200, 200, 200]
        lower = [0.001, 0.001, 0.001, 0.001]
        model_fun = demo_models_kp.bottleneck
    elif model == "no_mig":
        # divergence with no migration
        num = 3
        p_labels = "nu1, nu2, T"
        upper = [150, 150, 15]
        lower = [0.001, 0.001, 0.001]
        model_fun = demo_models_kp.no_migration
    elif model == "snm":
        # standard neutral model, no divergence
        num = 1
        p_labels = "nu"
        upper = [150]
        lower = [0.001]
        model_fun = demo_models_kp.no_divergence
    elif model == "sym_mig":
        # divergence with symmetrical migration
        num = 4
        p_labels = "nu1, nu2, m, T"
        upper = [150, 150, 10, 15]
        lower = [0.001, 0.001, 0.001, 0.001]
        model_fun = demo_models_kp.sym_migration
    elif model == "asym_mig":
        # divergence with asymmetrical migration
        num = 5
        p_labels = "nu1, nu2, m12, m21, T"
        upper = [150, 150, 10, 10, 15]
        lower = [0.001, 0.001, 0.001, 0.001, 0.001]
        model_fun = demo_models_kp.asym_migration
    elif model == "anc_sym_mig":
        num = 5
        p_labels = "nu1, nu2, m, T1, T2"
        upper = [150, 150, 10, 15, 15]
        lower = [0.001, 0.001, 0.001, 0.001, 0.001]
        model_fun = demo_models_kp.anc_sym_migration
    elif model == "anc_asym_mig":
        num = 6
        p_labels = "nu1, nu2, m12, m21, T1, T2"
        upper = [150, 150, 10, 10, 15, 15]
        lower = [0.001, 0.001, 0.001, 0.001, 0.001, 0.001]
        model_fun = demo_models_kp.anc_asym_migration
    elif model == "sec_cont_sym_mig":
        num = 5
        p_labels = "nu1, nu2, m, T1, T2"
        upper = [150, 150, 10, 15, 15]
        lower = [0.001, 0.001, 0.001, 0.001, 0.001]
        model_fun = demo_models_kp.sec_contact_sym_migration
    elif model == "sec_cont_asym_mig":
        num = 6
        p_labels = "nu1, nu2, m12, m21, T1, T2"
        upper = [150, 150, 10, 10, 15, 15]
        lower = [0.001, 0.001, 0.001, 0.001, 0.001, 0.001]
        model_fun = demo_models_kp.sec_contact_asym_migration
    elif model == "iso_inbred":
        # divergence with inbreeding
        num = 5
        p_labels = "nu1, nu2, F1, F2, T"
        upper = [150, 150, 1, 1, 15]
        lower = [0.001, 0.001, 0.00001, 0.00001, 0.001]
        model_fun = demo_models_kp.iso_inbreeding
    elif model == "mig_inbred":
        # divergence with migration with inbreeding
        num = 6
        p_labels = "nu1, nu2, F1, F2, m, T"
        upper = [150, 150, 1, 1, 10, 15]
        lower = [0.001, 0.001, 0.00001, 0.00001, 0.001, 0.001]
        model_fun = demo_models_kp.mig_inbreeding
    elif model == "anc_mig_inbred":
        # ancient migration with inbreeding
        num = 7
        p_labels = "nu1, nu2, F1, F2, m, T1, T2"
        upper = [150, 150, 1, 1, 10, 15, 15]
        lower = [0.001, 0.001, 0.00001, 0.00001, 0.001, 0.001, 0.001]
        model_fun = demo_models_kp.anc_sym_mig_inbred
    elif model == "sec_cont_inbred":
        # secondary contact with inbreeding
        num = 7
        p_labels = "nu1, nu2, F1, F2, m, T1, T2"
        upper = [150, 150, 1, 1, 10, 15, 15]
        lower = [0.001, 0.001, 0.00001, 0.00001, 0.001, 0.001, 0.001]
        model_fun = demo_models_kp.sec_contact_sym_mig_inbred
    elif model == "mig_be_inbred":
        # population size changes with divergence with migration and inbreeding
        num = 10
        p_labels = "nu1, nu2, nu1a, nu2a, F1, F2, m1, m2,  T1, T2"
        upper = [150, 150, 150, 150, 1, 1, 10, 10, 15, 15]
        lower = [0.001, 0.001, 0.001, 0.001, 0.00001, 0.00001, 0.001, 0.001, 0.001, 0.001]
        model_fun = demo_models_kp.mig_be_inbred
    else:
        print("model nickname undefined please check you are using the correct model nickname!")

    # Create log file.
    out_name = "../results/dadi_optimisation_ahya2.txt"
    with open(out_name, 'a') as opt_out:
        if opt_out.tell() == 0:
            print('Creating a new file\n')
            opt_out.write(
                "Pop\tModel\tFolds\tlog-likelihood\tAIC\tchi-squared\ttheta\tinitial_params\toptimised_params\t"
                "optimised_params_labels\n")
        else:
            print('File exists, appending\n')

    # This is our initial guess for the parameters, which is somewhat arbitrary.
    if int_params is None:
        p0 = [1] * num
    else:
        p0 = int_params
    # Paste below your optimised params to start from a specified place.
    # For example:
    # if model == "iso_inbred":
    #    if snps == "AG1-AG2":
    #        p1 = [2.1498, 97.3976, 0.0374, 0.0231, 0.4141]
    #    elif snps == "AL1-AL2":
    #        p1 = [7.4065, 1.9955, 0.0302, 0.0148, 0.4264]
    #    elif snps == "AG1-AL1":
    #        p1 = [0.8382, 41.4364, 0.0371, 0.0229, 0.8367]
    #    elif snps == "AG1-AL2":
    #        p1 = [1.0488, 46.8787, 0.0234, 0.00001, 0.8585]
    #    elif snps == "AG2-AL1":
    #        p1 = [0.7615, 55.4268, 0.0261, 0.0712, 0.9558]
    #    elif snps == "AG2-AL2":
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
    theta = numpy.around(theta, 2)

    # Calculate likelihood
    ll = dadi.Inference.ll_multinom(sim_model, data)
    ll = numpy.around(ll, 2)

    # Calculate AIC
    aic = (-2 * (float(ll))) + (2 * len(param_opt))

    # Calculate  Chi^2
    scaled_sim_model = sim_model * theta
    folded_sim_model = scaled_sim_model.fold()
    chi2 = numpy.sum((folded_sim_model - data) ** 2 / folded_sim_model)
    chi2 = numpy.around(chi2, 2)

    # Store results with likelihoods, theta and parameters (max. likelihood)
    results = [ll, aic, chi2, theta, param_opt]

    # Export the results to file
    with open(out_name, "a") as opt_out:
        easy_p = ",".join([str(numpy.around(x, 4)) for x in results[4]])
        int_p = ",".join([str(numpy.around(x, 4)) for x in p0])
        opt_out.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\n".format(snps, model, folds,
                                                                                  results[0], results[1], results[2],
                                                                                  results[3], int_p, easy_p, p_labels))
    print('* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *')
    print('* * * * * * * * * * * * * * * * * *  Finished optimisation  * * * * * * * * * * * * * * * * * *')
    print('* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *')


if __name__ == '__main__':
    # Arguments
    parser = argparse.ArgumentParser(prog="dadi Optimisation", usage='[options]')
    parser.add_argument("snps")
    parser.add_argument("model")
    parser.add_argument("masked")
    parser.add_argument("method")
    parser.add_argument("folds", type=int, nargs=1)
    parser.add_argument("-p", "--int_params", nargs="+", type=float)
    # parser.add_argument("out_path")
    args = parser.parse_args()

    # Setting variables
    snps = args.snps
    model = args.model
    masked = args.masked
    method = args.method
    folds = args.folds
    int_params = args.int_params

    # Need to manually define!
    # Define optimisation bounds
    PTS = [100, 120, 130]

    # If you are wanting to export data to a specific location,
    # uncomment the proceeding comment and argument parse,
    # then add path variable to main function.
    # path = "{}".format(args.out_path)

    main(snps, model, masked, folds, int_params, PTS, method=method)
