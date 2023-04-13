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


def main(snps, model, vmin, opt, PTS, figsize, figsize2):
    # Parse the data file to generate the sfs
    snp_path = '../data/fs/' + snps + ".fs"
    pops = "{}".format(snps)
    fs = Spectrum.from_file(snp_path)

    # Name of plots
    out_name = "../plots/" + snps

    # Plot font size
    plt.rcParams.update({'font.size': 8})

    # Masking the spectrum
    fs.mask[1, 0] = True
    fs.mask[0, 1] = True
    fs.mask[2, 0] = True
    fs.mask[0, 2] = True
    fs.mask[1, 1] = True

    # Print spectrum statistics
    print("\n* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n")
    print("Sum of SFS: {}".format(numpy.around(fs.S(), 2)))
    print("FST: {}".format(numpy.around(fs.Fst(), 2)))
    print("\n* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n")

    # Plot the data fs
    fig1 = pylab.figure(figsize=figsize)
    Plotting.plot_single_2d_sfs(fs, vmin=vmin, cmap=pylab.cm.hsv)
    fig1.tight_layout()
    if os.path.isfile(out_name + "_data.pdf"):
        print("data file exists")
    else:
        fig1.savefig(out_name + "_data.pdf", dpi=300)

    # Choose the model
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
        print("Choose correct model name")

    # Create an extrapolating function
    func_exec = Numerics.make_extrap_log_func(model_fun)

    # Simulate the model with the optimised parameters
    sim_model = func_exec(opt, fs.sample_sizes, PTS)
    theta = Inference.optimal_sfs_scaling(sim_model, fs)
    scaled_sim_model = sim_model * theta
    folded_sim_model = scaled_sim_model.fold()

    # Calculate the residuals
    resid = Inference.Anscombe_Poisson_residual(folded_sim_model, fs)

    # Plot figure to (residuals)
    fig2 = pylab.figure(figsize=figsize)
    Plotting.plot_2d_resid(resid, resid_range=3)
    fig2.tight_layout()
    fig2.savefig(out_name + model + "_residual.pdf", dpi=300)

    # Plot figure 3 (the simulated model)
    fig3 = pylab.figure(figsize=figsize)
    Plotting.plot_single_2d_sfs(folded_sim_model, vmin=vmin, cmap=pylab.cm.hsv)
    fig3.tight_layout()
    fig3.savefig(out_name + args.model + "_model.pdf", dpi=300)

    # Plot figure 4 (all together and distribution of residuals)
    fig4 = pylab.figure(figsize=figsize2)
    Plotting.plot_2d_comp_multinom(folded_sim_model, fs, resid_range=3, vmin=vmin)
    fig4.tight_layout()
    fig4.savefig(out_name + model + "_all4.pdf", dpi=300)


if __name__ == '__main__':
    # Arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("snps")
    parser.add_argument("model")
    parser.add_argument("vmin", type=float)
    parser.add_argument("-o", "--opt", nargs='+', type=float)
    args = parser.parse_args()

    # Setting variables
    snps = args.snps
    model = args.model
    vmin = args.vmin
    opt = args.opt

    # Extrapolating grid size
    PTS = [50, 60, 70]

    # Figure sizes
    figsize = (2.5, 2)
    figsize2 = (5, 4)

    main(snps, model, vmin, opt, PTS, figsize, figsize2)
