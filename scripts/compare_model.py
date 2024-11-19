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
import SETTINGS


def main(data, model, mask, fold, vmin, opt, PTS, figsize, figsize2, resid_range):
    # Parse the data file to generate the sfs
    snp_path = '../data/fs/' + data + ".fs"
    fs = Spectrum.from_file(snp_path)

    # Name of plots
    out_name = "../plots/models/" + data

    # Plot font size
    plt.rcParams.update({'font.size': 8})

    extras = ""

    # Masking the spectrum
    if mask == "low" and len(fs.sample_sizes) == 2:
        fs.mask[1, 0] = True
        fs.mask[0, 1] = True
        fs.mask[2, 0] = True
        fs.mask[0, 2] = True
        fs.mask[1, 1] = True
        print("Singletons and doubletons masked")
        extras += "mask{}".format(mask)
    elif mask == "low" and len(fs.sample_sizes) == 1:
        fs.mask[1] = True
        fs.mask[2] = True
        extras += "mask{}".format(mask)
    elif mask == "mid" and len(fs.sample_sizes) == 2:
        mid0 = int(fs.sample_sizes[0]/2)
        mid1 = int(fs.sample_sizes[1]/2)
        fs.mask[mid0, :] = True
        fs.mask[:, mid1] = True
        print("Mid frequencies masked")
        extras += "mask{}".format(mask)
    elif mask == "mid" and len(fs.sample_sizes) == 1:
        mid = fs.sample_sizes/2
        mid = int(mid[0])
        fs.mask[mid] = True
        print("Mid frequencies masked")
        extras += "mask{}".format(mask)
    elif mask == "both" and len(fs.sample_sizes) == 2:
        mid0 = int(fs.sample_sizes[0]/2)
        mid1 = int(fs.sample_sizes[1]/2)
        fs.mask[mid0, :] = True
        fs.mask[:, mid1] = True
        fs.mask[1, 0] = True
        fs.mask[0, 1] = True
        fs.mask[2, 0] = True
        fs.mask[0, 2] = True
        fs.mask[1, 1] = True
        extras += "maskmidlow"
        print("Mid frequencies and singletons and doubletons masked")
    else:
        print("Spectrum not masked")

    model_fun = SETTINGS.get_settings(model)

    # Print spectrum statistics
    print("\n* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n")
    print("Sum of SFS: {}".format(numpy.around(fs.S(), 2)))
    if len(fs.sample_sizes) == 2:
        print("FST: {}".format(numpy.around(fs.Fst(), 2)))
    elif len(fs.sample_sizes) == 1:
        print("Tajima's D: {}".format(numpy.around(fs.Tajima_D()), 2))
        print("Nucleotide diversity: {}".format(numpy.around(fs.pi()), 2))
    print("\n* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n")

    # Create an extrapolating function
    func_exec = Numerics.make_extrap_log_func(model_fun)

    # Simulate the model with the optimised parameters
    sim_model = func_exec(opt, fs.sample_sizes, PTS)
    theta = Inference.optimal_sfs_scaling(sim_model, fs)
    scaled_sim_model = sim_model * theta
    if fold == "folded":
        folded_sim_model = scaled_sim_model.fold()
        # Calculate the residuals
        resid = Inference.Anscombe_Poisson_residual(folded_sim_model, fs)
        extras += "folded"
    elif fold == "unfolded":
        resid = Inference.Anscombe_Poisson_residual(scaled_sim_model, fs)
        extras += "unfolded"
    else:
        print("choose folded or unfolded for folding")

    if len(fs.sample_sizes) == 2:
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
        Plotting.plot_2d_resid(resid, resid_range=resid_range)
        fig2.tight_layout()
        fig2.savefig(out_name + "_" + model + "_" + extras + "_residual.pdf", dpi=300)

        # Plot figure 3 (the simulated model)
        fig3 = pylab.figure(figsize=figsize)
        if fold == "folded":
            Plotting.plot_single_2d_sfs(folded_sim_model, vmin=vmin, cmap=pylab.cm.hsv)
        elif fold == "unfolded":
            Plotting.plot_single_2d_sfs(scaled_sim_model, vmin=vmin, cmap=pylab.cm.hsv)
        else:
            raise ValueError("choose folding")
        fig3.tight_layout()
        fig3.savefig(out_name + "_" + model + "_" + extras + "_model.pdf", dpi=300)

        # Plot figure 4 (all together and distribution of residuals)
        fig4 = pylab.figure(figsize=figsize2)
        if fold == "folded":
            Plotting.plot_2d_comp_multinom(folded_sim_model, fs, vmin=vmin)
        elif fold == "unfolded":
            Plotting.plot_2d_comp_multinom(scaled_sim_model, fs, vmin=vmin)
        else:
            raise ValueError("choose folding")
        fig4.tight_layout()
        fig4.savefig(out_name + model + "_" + extras + "_all4.pdf", dpi=300)
    elif len(fs.sample_sizes) == 1:
        # Plot the data fs
        fig1 = pylab.figure(figsize=figsize)
        Plotting.plot_1d_fs(fs)
        fig1.tight_layout()
        if os.path.isfile(out_name + "_" + extras + "_data.pdf"):
            print("data file exists")
        else:
            fig1.savefig(out_name + "_" + extras + "_data.pdf", dpi=300)
        # Fig 2 Poisson
        fig2 = pylab.figure(figsize=figsize)
        if fold == "folded":
            Plotting.plot_1d_comp_Poisson(folded_sim_model, fs)
        elif fold == "unfolded":
            Plotting.plot_1d_comp_Poisson(scaled_sim_model, fs)
        else:
            raise ValueError("choose folding")
        fig2.tight_layout()
        fig2.savefig(out_name + "_" + model + "_" + extras + "_poisson_residual.pdf", dpi=300)

        # Fig 3 multimnom
        fig3 = pylab.figure(figsize=figsize)
        if fold == "folded":
            Plotting.plot_1d_comp_multinom(folded_sim_model, fs)
        elif fold == "unfolded":
            Plotting.plot_1d_comp_multinom(scaled_sim_model, fs)
        else:
            raise ValueError("choose folding")
        fig3.tight_layout()
        fig3.savefig(out_name + "_" + model + "_" + extras + "_mulitnom_residual.pdf", dpi=300)


if __name__ == '__main__':
    # Arguments
    parser = argparse.ArgumentParser(
        prog="compare dadi model to data",
        description="A script for comparing the optimised dadi model compared to the data",
        usage="%(prog)s [options] <data> <model> <mask> <vmin> <fold> <-o opt>"
    )

    # Required positional arguments
    parser.add_argument(
        "data",
        help="Name of the fs file."
    )
    parser.add_argument(
        "model",
        help="The name of the model used for optimisation."
    )
    parser.add_argument(
        "mask",
        help="Type of masking to apply (e.g., 'low', 'mid')."
    )
    parser.add_argument(
        "fold",
        help="Whether the spectra is folded or unfolded (e.g., 'folded' or 'unfolded')."
    )
    parser.add_argument(
        "vmin",
        type=float,
        help="The scale bar for the SNPs"
    )
    parser.add_argument(
        "-o", "--opt",
        nargs='+',
        type=float,
        help="Optimised parameters provided as space-separated floats (e.g., -o 1.0 2.0 3.0)."
    )

    # Parse arguments
    args = parser.parse_args()

    # Setting variables
    data = args.data
    model = args.model
    mask = args.mask
    fold = args.fold
    vmin = args.vmin
    opt = args.opt

    # Example usage of variables
    print(f"Data file: {data}")
    print(f"Model: {model}")
    print(f"Mask type: {mask}")
    print(f"Fold type: {fold}")
    print(f"Scale bar: {vmin}")
    print(f"Optimised parameters: {opt}")

    # Extrapolating grid size
    PTS = SETTINGS.SET_PTS

    # Figure sizes
    # Custom setting here..
    figsize = (2.5, 2)  # individual plots
    figsize2 = (5, 4)  # larger 4-way plot
    resid_range = 3  # custom residual range

    main(data, model, mask, fold, vmin, opt, PTS, figsize, figsize2, resid_range)
