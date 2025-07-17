#!/usr/anaconda3/env/dadi211/bin/python
# -*- coding: utf-8 -*-

"""
@author: Katharine Prata
@date created: 2/6/21
@description: Create non-parametric bootstraps using the subsample option.

script modified from YRI_CEU.py
"""

import argparse
import numpy
import os.path
import demo_models_kp
from dadi import Misc, Numerics, Inference
import SETTINGS
import plot_fs


def main(snps, model, sims, genotypes, chunk_size, opt, PTS, mask_type):
    # import the spectrum and popfile from data/vcf and data/popfile
    snp_path = "../data/vcf/" + snps + ".vcf"
    pop_path = "../data/popfile/pop_" + snps + ".txt"
    pops = "{}".format(snps)
    pop_ids = pops.split("-")

    # make a file with statistics about your bootstraps
    out_name = "../results/bootstraps/{}_{}_{}_nonparametric_bootstraps_vcf.txt".format(args.snps, args.model,
                                                                                        args.sims)
    with open(out_name, 'a') as stat_out:
        stat_out.write("Pops\tModel\tLikelihood\tAIC\tTheta\tsfs_sum\tchi-squared\tbootstrap\n")

    # Configuring haplotypes and genotypes
    if len(pop_ids) == 1:
        proj = [genotypes[0] * 2]
        subsample = {pop_ids[0]: genotypes[1]}
    elif len(pop_ids) == 2:
        proj = [genotypes[0] * 2, genotypes[1] * 2]
        subsample = {pop_ids[0]: genotypes[0], pop_ids[1]: genotypes[1]}
    else:
        proj = []
        subsample = {}
        for i in range(len(pop_ids)):
            proj.append(20)
            subsample["Pop{}".format(i)] = 10

        # Making subsampled bootstraps
        # Need to alter chunk_size here depending on how long your contigs are
    boots_subsample = Misc.bootstraps_subsample_vcf(snp_path, pop_path, chunk_size=chunk_size, Nboot=sims,
                                                    subsample=subsample,
                                                    pop_ids=pop_ids, polarized=False)
    # Saving your bootstraps to file
    for i in range(0, sims):
        bootstrap_out_name = "../results/bootstraps/{}_bootstrap_vcf_{}.fs".format(snps, i)
        if os.path.isfile(bootstrap_out_name):
            print("bootstrap exists")
        else:
            boots_subsample[i].to_file(bootstrap_out_name)

    # Place your models of choice here
    if model == "iso_inbred":
        model_fun = demo_models_kp.iso_inbreeding
    elif model == "mig_inbred":
        model_fun = demo_models_kp.mig_inbreeding
    elif model == "anc_mig":
        model_fun = demo_models_kp.anc_sym_mig_inbred
    elif model == "anc_asym_mig":
        model_fun = demo_models_kp.anc_asym_migration
    elif model == "sec_cont":
        model_fun = demo_models_kp.sec_contact_sym_mig_inbred
    else:
        model_fun = False
        print("Please specify the correct model: iso_inbred, mig_inbred, anc_mig or sec_cont.")

    # Calculate how well your bootstraps fit your optimised parameters
    for i in range(0, sims):
        # Use mask function from plot_fs.py
        plot_fs.apply_mask(boots_subsample[i], mask_type)
        func_exec = Numerics.make_extrap_func(model_fun)
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
            stat_out.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n".format(snps, model, ll, aic,
                                                                             theta, boots_subsample[i].S(), chisq, i))


if __name__ == '__main__':
    # Arguments
    parser = argparse.ArgumentParser(
        prog="vcf bootstraps",
        description="Create non-parametric bootstraps using the subsample option.",
        usage="%(prog)s [options] <snps> <model> <sims> <chunksize> [-g GENOTYPES] [-o OPT_PARAMS] [-m MASK]"
    )
    parser.add_argument(
        "snps",
        help="SNPs name (without extension)."
    )
    parser.add_argument(
        "model",
        help="Model to use for the analysis from kp_dadi."
    )
    parser.add_argument(
        "sims",
        type=int,
        help="Number of bootstraps to generate."
    )
    parser.add_argument(
        "chunksize",
        type=int,
        help="Chunk size for bootstrapping."
    )
    parser.add_argument(
        "-g", "--genotypes",
        nargs="+", type=int,
        help="Genotype counts for each population."
    )
    parser.add_argument(
        "-o", "--opt_params",
        nargs="+", type=float,
        help="Optimised parameters for the model."
    )
    parser.add_argument(
        "-m", "--mask",
        default="low",
        help="Type of masking to use (e.g., 'mid', 'low', or 'no'). Default is 'low'."
    )
    args: argparse.Namespace = parser.parse_args()

    # Setting variables
    snps = args.snps
    model = args.model
    sims = args.sims
    genotypes = args.genotypes
    chunk_size = args.chunksize
    opt = args.opt_params
    mask_type = args.mask

    # Import PTS from SETTINGS
    PTS = SETTINGS.SET_PTS

    main(snps, model, sims, genotypes, chunk_size, opt, PTS, mask_type)


