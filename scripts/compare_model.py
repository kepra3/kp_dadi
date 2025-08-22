#!/usr/anaconda3/env/dadi211/bin/python
# -*- coding: utf-8 -*-

"""
@author: Katharine Prata
@date created: 2/6/21
@description: Plotting the data, model, and residuals for comparison.

Compatible with Python 3.6.11 and dadi 2.1.1
"""

import argparse
import numpy as np
import matplotlib.pyplot as plt
import os
from dadi import Plotting, Spectrum, Numerics, Inference
import SETTINGS


def apply_mask(fs, mask_type):
    """Apply the specified mask to the frequency spectrum."""
    if mask_type == "low":
        low_mask = slice(1, 3)
        if len(fs.sample_sizes) == 2:
            fs.mask[low_mask, :] = True
            fs.mask[:, low_mask] = True
        elif len(fs.sample_sizes) == 1:
            fs.mask[low_mask] = True
        print("Low frequencies (singletons and doubletons) masked.")
    elif mask_type == "mid":
        mid_idx = [int(n / 2) for n in fs.sample_sizes]
        for idx, size in enumerate(fs.sample_sizes):
            if len(fs.sample_sizes) == 2:
                fs.mask[mid_idx[0], :] = True
                fs.mask[:, mid_idx[1]] = True
            else:
                fs.mask[mid_idx[0]] = True
        print("Mid frequencies masked.")
    elif mask_type == "both":
        apply_mask(fs, "low")
        apply_mask(fs, "mid")
        print("Low and mid frequencies masked.")
    else:
        print("No masking applied.")
    return mask_type


def plot_and_save(fig, filename, tight=True):
    """Save the figure to the specified filename."""
    if tight:
        fig.tight_layout()
    fig.savefig(filename, dpi=300)
    plt.close(fig)
    print(f"Saved plot to {filename}")


def plot_1d(fs, scaled_model, fold, out_name, extras):
    """Plot 1D spectra and residuals."""
    fig = plt.figure(figsize=(2.5, 2))
    Plotting.plot_1d_fs(fs)
    plot_and_save(fig, f"{out_name}_{extras}_data.pdf")

    fig = plt.figure(figsize=(2.5, 2))
    plot_func = Plotting.plot_1d_comp_Poisson if fold == "folded" else Plotting.plot_1d_comp_multinom
    plot_func(scaled_model, fs)
    plot_and_save(fig, f"{out_name}_{extras}_residual.pdf")


def plot_2d(fs, scaled_model, resid, fold, out_name, model, extras, resid_range):
    """Plot 2D spectra and residuals."""
    fig = plt.figure(figsize=(2.5, 2))
    Plotting.plot_single_2d_sfs(fs, vmin=1, cmap=plt.cm.hsv)
    plot_and_save(fig, f"{out_name}_data.pdf")

    fig = plt.figure(figsize=(2.5, 2))
    Plotting.plot_2d_resid(resid, resid_range=resid_range)
    plot_and_save(fig, f"{out_name}_{model}_{extras}_residual.pdf")

    fig = plt.figure(figsize=(5, 4))
    Plotting.plot_2d_comp_multinom(scaled_model, fs, vmin=1)
    plot_and_save(fig, f"{out_name}_{model}_{extras}_all.pdf")


def main(data, model, mask, fold, opt, PTS, resid_range):
    """Main execution logic."""
    fs = Spectrum.from_file(f"../data/fs/{data}.fs")
    out_name = f"../plots/models/{data}"
    plt.rcParams.update({'font.size': 8})

    extras = apply_mask(fs, mask)

    # Print statistics
    print("\nSpectrum Statistics:")
    print(f"Sum of SFS: {np.around(fs.S(), 2)}")
    if len(fs.sample_sizes) == 2:
        print(f"FST: {np.around(fs.Fst(), 2)}")
    elif len(fs.sample_sizes) == 1:
        print(f"Tajima's D: {np.around(fs.Tajima_D(), 2)}")
        print(f"Nucleotide diversity: {np.around(fs.pi(), 2)}")
    print()

    # Simulate the model
    model_func = Numerics.make_extrap_log_func(SETTINGS.get_settings(model))
    sim_model = model_func(opt, fs.sample_sizes, PTS)
    theta = Inference.optimal_sfs_scaling(sim_model, fs)
    scaled_model = sim_model * theta
    resid = Inference.Anscombe_Poisson_residual(scaled_model.fold() if fold == "folded" else scaled_model, fs)

    # Plot based on sample sizes
    if len(fs.sample_sizes) == 2:
        plot_2d(fs, scaled_model, resid, fold, out_name, model, extras, resid_range)
    elif len(fs.sample_sizes) == 1:
        plot_1d(fs, scaled_model, fold, out_name, extras)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compare dadi model to data.")
    parser.add_argument("data", help="Name of the fs file.")
    parser.add_argument("model", help="Optimized model name.")
    parser.add_argument("mask", help="Type of masking to apply ('low', 'mid', 'both', or 'none').")
    parser.add_argument("fold", choices=["folded", "unfolded"], help="Fold type.")
    parser.add_argument("resid_range", type=float, help="Scale bar for the SNPs.")
    parser.add_argument("-o", "--opt", nargs='+', type=float, required=True, help="Optimized parameters.")

    args = parser.parse_args()

    PTS = SETTINGS.SET_PTS
    main(args.data, args.model, args.mask, args.fold, args.opt, PTS, args.resid_range)
