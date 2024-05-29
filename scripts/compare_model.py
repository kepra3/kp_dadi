#!/usr/anaconda3/env/dadi211/bin/python
# -*- coding: utf-8 -*-

"""
@author: Katharine Prata
@date created: 2/6/21
@description: Plotting the data, model and residuals for comparison.

Compatible with python 3.6.11 and dadi 2.1.1
"""

import argparse
import numpy
import pylab
import demo_models_kp
import os
from dadi import Plotting, Spectrum, Numerics, Inference
import matplotlib.pyplot as plt


def main(snps, model, mask, fold, vmin, opt, PTS, figsize, figsize2):
    # Parse the data file to generate the sfs
    snp_path = '../data/fs/' + snps + ".fs"
    pops = "{}".format(snps)
    fs = Spectrum.from_file(snp_path)

    # Name of plots
    out_name = "../plots/" + snps

    # Plot font size
    plt.rcParams.update({'font.size': 8})

    # Masking the spectrum
    if mask == "low" and len(fs.sample_sizes) == 2:
        fs.mask[1, 0] = True
        fs.mask[0, 1] = True
        fs.mask[2, 0] = True
        fs.mask[0, 2] = True
        fs.mask[1, 1] = True
        print("Singletons and doubletons masked")
    elif mask == "mid" and len(fs.sample_sizes) == 1:
        mid = fs.sample_sizes/2
        mid = int(mid[0])
        fs.mask[mid] = True
        print("Mid frequencies masked")
    else:
        print("Spectrum not masked")

    # Choose the model
    if model == "iso_inbred":
        model_fun = demo_models_kp.iso_inbreeding
        dim = "2d"
    elif model == "mig_inbred":
        model_fun = demo_models_kp.mig_inbreeding
        dim = "2d"
    elif model == "no_mig":
        model_fun = demo_models_kp.no_migration
        dim = "2d"
    elif model == "asym_mig":
        model_fun = demo_models_kp.asym_migration
        dim = "2d"
    elif model == "anc_mig":
        model_fun = demo_models_kp.anc_sym_mig_inbred
        dim = "2d"
    elif model == "sec_cont":
        model_fun = demo_models_kp.sec_contact_sym_mig_inbred
        dim = "2d"
    elif model == "anc_asym_mig":
        model_fun = demo_models_kp.anc_asym_migration
        dim = "2d"
    elif model == "sec_cont_asym_mig":
        model_fun = demo_models_kp.sec_contact_asym_migration
        dim = "2d"
    elif model == "split_bottle":
        model_fun = demo_models_kp.split_bottlegrowth
        dim = "2d"
    elif model == "split_bottle_asym_mig":
        model_fun = demo_models_kp.split_bottlegrowth_asym_mig
        dim = "2d"
    elif model == "split_bottle_anc_asym_mig":
        model_fun = demo_models_kp.split_bottlegrowth_ancient_asym_mig
        dim = "2d"
    elif model == "split_bottle_sec_asym_mig":
        model_fun = demo_models_kp.split_bottlegrowth_second_asym_mig
        dim = "2d"
    elif model == "split_sizechange":
        model_fun = demo_models_kp.split_sizechange
        dim = "2d"
    elif model == "split_sizechange_asym_mig":
        model_fun = demo_models_kp.split_sizechange_asym_mig
        dim = "2d"
    elif model == "split_sizechange_anc_asym_mig":
        model_fun = demo_models_kp.split_sizechange_ancient_asym_mig
        dim = "2d"
    elif model == "split_sizechange_sec_asym_mig":
        model_fun = demo_models_kp.split_sizechange_second_asym_mig
        dim = "2d"
    elif model == "snm.1d":
        model_fun = demo_models_kp.no_divergence_1d
        dim = "1d"
    elif model == "size_change":
        model_fun = demo_models_kp.instant_change
        dim = "1d"
    elif model == "bottle":
        model_fun = demo_models_kp.bottlegrowth
        dim = "1d"
    elif model == "bottle_neck":
        model_fun = demo_models_kp.bottleneck
        dim = "1d"
    else:
        model_fun = False
        print("Choose correct model name")

    # Print spectrum statistics
    print("\n* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n")
    print("Sum of SFS: {}".format(numpy.around(fs.S(), 2)))
    if dim == "2d":
        print("FST: {}".format(numpy.around(fs.Fst(), 2)))
    elif dim == "1d":
        print("Tajima's D: {}".format(numpy.around(fs.Tajima_D()), 2))
        print("Nucleotide diversity: {}".format(numpy.around(fs.pi()), 2))
    print("\n* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n")

    # Create an extrapolating function
    func_exec = Numerics.make_extrap_log_func(model_fun)

    # Simulate the model with the optimised parameters
    sim_model = func_exec(opt, fs.sample_sizes, PTS)
    theta = Inference.optimal_sfs_scaling(sim_model, fs)
    scaled_sim_model = sim_model * theta
    if fold == "yes":
        folded_sim_model = scaled_sim_model.fold()
        # Calculate the residuals
        resid = Inference.Anscombe_Poisson_residual(folded_sim_model, fs)
    elif fold == "no":
        resid = Inference.Anscombe_Poisson_residual(scaled_sim_model, fs)
    else:
        print("choose yes or no for folding")

    if dim == "2d":
        # Plot the data fs
        fig1 = pylab.figure(figsize=figsize)
        Plotting.plot_single_2d_sfs(fs, vmin=vmin, cmap=pylab.cm.hsv)
        fig1.tight_layout()
        if os.path.isfile(out_name + "_data.pdf"):
            print("data file exists")
        else:
            fig1.savefig(out_name + "_data.pdf", dpi=300)
        # Plot figure to (residuals)
        fig2 = pylab.figure(figsize=figsize)
        Plotting.plot_2d_resid(resid, resid_range=3)
        fig2.tight_layout()
        fig2.savefig(out_name + "_" + model + "_residual.pdf", dpi=300)

        # Plot figure 3 (the simulated model)
        fig3 = pylab.figure(figsize=figsize)
        if fold == "yes":
            Plotting.plot_single_2d_sfs(folded_sim_model, vmin=vmin, cmap=pylab.cm.hsv)
        elif fold == "no":
            Plotting.plot_single_2d_sfs(scaled_sim_model, vmin=vmin, cmap=pylab.cm.hsv)
        else:
            print("choose folding")
        fig3.tight_layout()
        fig3.savefig(out_name + "_" + model + "_model.pdf", dpi=300)

        # Plot figure 4 (all together and distribution of residuals)
        fig4 = pylab.figure(figsize=figsize2)
        if fold == "yes":
            Plotting.plot_2d_comp_multinom(folded_sim_model, fs, resid_range=3, vmin=vmin)
        elif fold == "no":
            Plotting.plot_2d_comp_multinom(scaled_sim_model, fs, resid_range=3, vmin=vmin)
        else:
            print("choose folding")
        fig4.tight_layout()
        fig4.savefig(out_name + model + "_all4.pdf", dpi=300)
    elif dim == "1d":
        # Plot the data fs
        fig1 = pylab.figure(figsize=figsize)
        Plotting.plot_1d_fs(fs)
        fig1.tight_layout()
        if os.path.isfile(out_name + "_data.pdf"):
            print("data file exists")
        else:
            fig1.savefig(out_name + "_data.pdf", dpi=300)
        # Fig 2 Poisson
        fig2 = pylab.figure(figsize=figsize)
        if fold == "yes":
            Plotting.plot_1d_comp_Poisson(folded_sim_model, fs)
        elif fold == "no":
            Plotting.plot_1d_comp_Poisson(scaled_sim_model, fs)
        else:
            print('choose folding')
        fig2.tight_layout()
        fig2.savefig(out_name + "_" + model + "_poisson_residual.pdf", dpi=300)

        # Fig 3 multimnom
        fig3 = pylab.figure(figsize=figsize)
        if fold == "yes":
            Plotting.plot_1d_comp_multinom(folded_sim_model, fs)
        elif fold == "no":
            Plotting.plot_1d_comp_multinom(scaled_sim_model, fs)
        else:
            print('choose folding')
        fig3.tight_layout()
        fig3.savefig(out_name + "_" + model + "_mulitnom_residual.pdf", dpi=300)


if __name__ == '__main__':
    # Arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("snps")
    parser.add_argument("model")
    parser.add_argument("mask")
    parser.add_argument("fold")
    parser.add_argument("-o", "--opt", nargs='+', type=float)
    args = parser.parse_args()

    # Setting variables
    snps = args.snps
    model = args.model
    mask = args.mask
    fold = args.fold
    opt = args.opt

    # Extrapolating grid size
    PTS = [100, 120, 130]

    # Figure sizes
    figsize = (2.5, 2)
    figsize2 = (5, 4)
    vmin = 0.05

    main(snps, model, mask, fold, vmin, opt, PTS, figsize, figsize2)
