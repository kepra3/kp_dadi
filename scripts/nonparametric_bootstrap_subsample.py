#!/usr/anaconda3/bin/python
# -*- coding: utf-8 -*-
# remember to edit python environment if required.

"""
@author: kprata
@date created: 2/6/21
@description: Create non-parametric bootstraps using the subsample option.
"""

import argparse
import numpy
import os.path
import demo_models_kp
from dadi import Misc, Numerics, Inference

parser = argparse.ArgumentParser(prog="vcf bootstraps")
parser.add_argument("snps")
parser.add_argument("model")
parser.add_argument("sims", type=int)
parser.add_argument("-g", "--genotypes", nargs="+", type=int)
parser.add_argument("-o", "--opt_params", nargs="+", type=float)
args: argparse.Namespace = parser.parse_args()

# import the spectrum and popfile from data/vcf and data/popfile
snp_path = "../data/vcf/" + args.snps + ".vcf"
pop_path = "../data/popfile/pop_" + args.snps + ".txt"
pops = "{}".format(args.snps)
pop_ids = pops.split("-")

# make a file with statistics about your bootstraps
out_name = "../results/bootstraps/{}_{}_{}_nonparametric_bootstraps_vcf.txt".format(args.snps, args.model, args.sims)
with open(out_name, 'a') as stat_out:
    stat_out.write("Pops\tModel\tLikelihood\tAIC\tTheta\tsfs_sum\tchi-squared\tbootstrap\n")

# Configuring haplotypes and genotypes
if len(pop_ids) == 1:
    proj = [args.genotypes[0] * 2]
    subsample = {pop_ids[0]: args.genotypes[1]}
elif len(pop_ids) == 2:
    proj = [args.genotypes[0] * 2, args.genotypes[1] * 2]
    subsample = {pop_ids[0]: args.genotypes[0], pop_ids[1]: args.genotypes[1]}
else:
    proj = []
    subsample = {}
    for i in range(len(pop_ids)):
        proj.append(20)
        subsample["Pop{}".format(i)] = 10

# Making subsampled bootstraps
# May need to alter chunk_size here depending on how long your contigs are
boots_subsample = Misc.bootstraps_subsample_vcf(snp_path, pop_path, chunk_size=100, Nboot=args.sims,
                                                subsample=subsample,
                                                pop_ids=pop_ids, polarized=False)
# Saving your bootstraps to file
sims = args.sims
for i in range(0, sims):
    bootstrap_out_name = "../results/bootstraps/{}_bootstrap_vcf_{}.fs".format(args.snps, i)
    if os.path.isfile(bootstrap_out_name):
        print("bootstrap exists")
    else:
        boots_subsample[i].to_file("../results/bootstraps/{}_bootstrap_vcf_{}.fs".format(args.snps, i))

# Place your models of choice here
model = args.model
if model == "iso_inbred":
    model = demo_models_kp.iso_inbreeding
elif model == "mig_inbred":
    model = demo_models_kp.mig_inbreeding
elif model == "anc_mig":
    model = demo_models_kp.anc_sym_mig_inbred
elif model == "sec_cont":
    model = demo_models_kp.sec_contact_sym_mig_inbred
else:
    model = False
    print("Please specify the correct model: iso_inbred, mig_inbred, anc_mig or sec_cont.")

# Calculate how well your bootstraps fit your optimised parameters
for i in range(0, args.sims):
    # Define optimisation bounds - use the same ones as before
    PTS = [50, 60, 70]
    # Uncomment below if not masking singletons and doubletons
    boots_subsample[i].mask[1, 0] = True
    boots_subsample[i].mask[0, 1] = True
    boots_subsample[i].mask[2, 0] = True
    boots_subsample[i].mask[0, 2] = True
    opt = args.opt_params
    func_exec = Numerics.make_extrap_func(model)
    sim_model = func_exec(opt, boots_subsample[i].sample_sizes, PTS)
    ll = Inference.ll_multinom(sim_model, boots_subsample[i])
    ll = numpy.around(ll, 2)
    aic = (-2 * (float(ll))) + (2 * len(opt))
    theta = Inference.optimal_sfs_scaling(sim_model, boots_subsample[i])
    scaled_sim_model = sim_model * theta
    folded_sim_model = scaled_sim_model.fold()
    chisq = numpy.sum((folded_sim_model - boots_subsample[i]) ** 2 / folded_sim_model)
    chisq = numpy.around(chisq, 2)
    print("\n* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n")
    print("{}".format(i))
    print("\nData for site frequency spectrum:\n")
    print("Sample sizes: {}".format(boots_subsample[i].sample_sizes))
    print("Sum of SFS: {}".format(numpy.around(boots_subsample[i].S(), 2)))
    print("FST of SFS: {}".format(numpy.around(boots_subsample[i].Fst(), 2)))
    print("\n* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n")
    print("\nAnalysis sample size statistics:\n")
    print("Likelihood = {}".format(ll))
    print("AIC = {}".format(aic))
    print("Theta = {}".format(theta))
    print("Chi squared = {}".format(chisq))
    print("\n* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n")
    # make a file to get the likelihoods of all the nonparametric bootstraps
    with open(out_name, 'a') as stat_out:
        stat_out.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n".format(args.snps, args.model, ll, aic,
                                                                         theta, boots_subsample[i].S(), chisq, i))
