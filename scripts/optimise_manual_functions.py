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
method = subsample projection or none
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

parser = argparse.ArgumentParser(prog="dadi Optimisation", usage='[options]')
parser.add_argument("snps")
parser.add_argument("model")
parser.add_argument("masked")
parser.add_argument("method")
parser.add_argument("folds", type=int)
# parser.add_argument("out_path")
args = parser.parse_args()

# Import and define data constants
if args.method == "subsample":
    data = dadi.Spectrum.from_file('../data/{}_subsampled.fs'.format(args.snps))
elif args.method == "projection":
    data = dadi.Spectrum.from_file('../data/{}_projected.fs'.format(args.snps))
else:
    data = dadi.Spectrum.from_file('../data/{}.fs'.format(args.snps))
pops = "{}".format(args.snps)
pop_ids = pops.split("-")

# If you are wanting to export data to a specific location,
# uncomment the proceeding comment.
# path = "{}".format(args.out_path)

# Define optimisation bounds
PTS = [50, 60, 70]

# Mask singletons and doubletons
if args.masked == "yes":
    data.mask[1, 0] = True
    data.mask[0, 1] = True
    data.mask[2, 0] = True
    data.mask[0, 2] = True
else:
    print("Singletons and Doubletons are not masked")

# Print useful information about the sfs
print("The datafile is {}".format(args.snps))
print("\n* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n")
print("Data for site frequency spectrum:")
print("Sample sizes: {}".format(data.sample_sizes))
print("Sum of SFS: {}".format(numpy.around(data.S(), 2)))
print("FST of SFS: {}".format(numpy.around(data.Fst(), 2)))
print("\n* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n")

# Define model from argument
model = args.model

# Define metadata for models
if model == "iso_inbred":
    num = 5
    p_labels = "nu1, nu2, F1, F2, T"
    upper = [150, 150, 1, 1, 15]
    lower = [0.001, 0.001, 0.00001, 0.00001, 0.001]
    model = demo_models_kp.iso_inbreeding
elif model == "mig_inbred":
    num = 6
    p_labels = "nu1, nu2, F1, F2, m, T"
    upper = [150, 150, 1, 1, 10, 15]
    lower = [0.001, 0.001, 0.00001, 0.00001, 0.001, 0.001]
    model = demo_models_kp.mig_inbreeding
elif model == "anc_mig":
    num = 7
    p_labels = "nu1, nu2, F1, F2, m, T1, T2"
    upper = [150, 150, 1, 1, 10, 15, 15]
    lower = [0.001, 0.001, 0.00001, 0.00001, 0.001, 0.001, 0.001]
    model = demo_models_kp.anc_sym_mig_inbred
elif model == "snm":
    num = 1
    p_labels = "nu"
    upper = [150]
    lower = [0.001]
    model = demo_models_kp.no_divergence
elif model == "mig_be":
    num = 10
    p_labels = "nu1, nu2, nu1a, nu2a, F1, F2, m1, m2,  T1, T2"
    upper = [150, 150, 150, 150, 1, 1, 10, 10, 15, 15]
    lower = [0.001, 0.001, 0.001, 0.001, 0.00001, 0.00001, 0.001, 0.001, 0.001, 0.001]
    model = demo_models_kp.mig_be
else:
    num = 7
    p_labels = "nu1, nu2, F1, F2, m, T1, T2"
    upper = [150, 150, 1, 1, 10, 15, 15]
    lower = [0.001, 0.001, 0.00001, 0.00001, 0.001, 0.001, 0.001]
    model = demo_models_kp.sec_contact_sym_mig_inbred

# Create log file.
out_name = "../results/dadi_optimisation.txt"
with open(out_name, 'a') as opt_out:
    if opt_out.tell() == 0:
        print('Creating a new file\n')
        opt_out.write(
            "Pop\tModel\tFolds\tlog-likelihood\tAIC\tchi-squared\ttheta\toptimised_params({})\n".format(p_labels))
    else:
        print('File exists, appending\n')

# This is our initial guess for the parameters, which is somewhat arbitrary.
p0 = [1] * num
# Paste below your optimised params to start from a specified place.
# For example:
# if args.model == "iso_inbred":
#    if args.snps == "AG1-AG2":
#        p1 = [2.1498, 97.3976, 0.0374, 0.0231, 0.4141]
#    elif args.snps == "AL1-AL2":
#        p1 = [7.4065, 1.9955, 0.0302, 0.0148, 0.4264]
#    elif args.snps == "AG1-AL1":
#        p1 = [0.8382, 41.4364, 0.0371, 0.0229, 0.8367]
#    elif args.snps == "AG1-AL2":
#        p1 = [1.0488, 46.8787, 0.0234, 0.00001, 0.8585]
#    elif args.snps == "AG2-AL1":
#        p1 = [0.7615, 55.4268, 0.0261, 0.0712, 0.9558]
#    elif args.snps == "AG2-AL2":
#        p1 = [0.6783, 44.3377, 0.055, 0.0398, 0.7721]
#    else:
#        p1 = p0
# if using custom starting parameter values then comment the line below
p1 = p0

# Comments from software example (Gutenkunst et al., 2009)
# Make the extrapolating version of our demographic model function.
func_ex = dadi.Numerics.make_extrap_log_func(model)

# Comments taken from Gutenkunst et al. (2009)
# Perturb our parameters before optimisation. This does so by taking each
# parameter a up to a factor of "folds" up or down.
p1 = dadi.Misc.perturb_params(p1, fold=args.folds, upper_bound=upper, lower_bound=lower)
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
                                             verbose=1, maxiter=1000)
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
    opt_out.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n".format(args.snps, args.model, args.folds,
                                                                    results[0], results[1], results[2],
                                                                    results[3], easy_p))
print('* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *')
print('* * * * * * * * * * * * * * * * * *  Finished optimisation  * * * * * * * * * * * * * * * * * *')
print('* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *')
