#!/usr/anaconda3/env/dadi211/bin/python
# -*- coding: utf-8 -*-

"""
@author: Katharine Prata
@date created: 23/4/21
@description: vcf file is imported through make_data_dict_vcf and an fs is created.

Input:
File: .vcf file named with population names.
Arguments:
snps = vcf file named after your two populations broken by a hyphen (-)
fold = unfolded or folded
masked = yes or no
method = subsample or projection
genotypes = the number you want to subsample by in units of genotypes
run from the scripts directory or change path variables

Output: a subsampled fs for input into dadi analysis, a text file with statistics and plots of fs produced.

Compatible with python 3.6.11 and dadi 2.1.1
"""

import argparse
import copy
from argparse import Namespace
import dadi
import numpy as np
import pylab
import matplotlib as plt


def main(snps, fold, masked, method, genotypes):
    # Import the spectrum and popfile from data/vcf/ and data/popfile/
    snp_path = "../data/vcf/" + snps + ".vcf"
    pop_path = "../data/popfile/pop_" + snps + ".txt"
    pops = "{}".format(snps)
    pop_ids = pops.split("-")

    # make a file with statistics about your sfs
    stats_out_name = "../results/sfs_stats-2.txt"
    with open(stats_out_name, 'a') as stats_out:
        if stats_out.tell() == 0:
            print('Creating a new file')
            stats_out.write("Pop\tSample sizes\tFold\tSum of SFS\tStat\n")
        else:
            print('file exists, appending')

    # Configuring haplotypes and genotypes
    if len(pop_ids) == 1:
        print("One dimensional sfs")
        proj = [genotypes[0] * 2]
        subsample = {pop_ids[0]: genotypes[0]}
    elif len(pop_ids) == 2:
        print("Two dimensional sfs")
        proj = [genotypes[0] * 2, genotypes[1] * 2]
        subsample = {pop_ids[0]: genotypes[0], pop_ids[1]: genotypes[1]}
    else:
        proj = []
        subsample = {}
        for i in range(len(pop_ids)):
            proj.append(20)
            subsample["Pop{}".format(i)] = 10

    if fold == "folded":
        if method == "subsample":
            dd = dadi.Misc.make_data_dict_vcf(snp_path, pop_path, subsample=subsample)
            fs = dadi.Spectrum.from_data_dict(dd, pop_ids=pop_ids, projections=proj, polarized=False)
            fs.to_file("../data/fs/{}_subsampled.fs".format(snps))
        elif method == "projection":
            dd = dadi.Misc.make_data_dict_vcf(snp_path, pop_path)
            fs = dadi.Spectrum.from_data_dict(dd, pop_ids=pop_ids, projections=proj, polarized=False)
            fs.to_file("../data/fs/{}_projected.fs".format(snps))
        elif method == "no":
            dd = dadi.Misc.make_data_dict_vcf(snp_path, pop_path)
            fs = dadi.Spectrum.from_data_dict(dd, pop_ids=pop_ids, projections=proj, polarized=False)
            fs.to_file("../data/fs/{}_projected.fs".format(snps))
        else:
            raise ValueError("Choose projection, subsample or no")
    elif fold == "unfolded":
        if method == "subsample":
            dd = dadi.Misc.make_data_dict_vcf(snp_path, pop_path, subsample=subsample)
            fs = dadi.Spectrum.from_data_dict(dd, pop_ids=pop_ids, projections=proj, polarized=True)
            fs.to_file("../data/fs/{}_unfolded_subsampled.fs".format(snps))
        elif method == "projection":
            dd = dadi.Misc.make_data_dict_vcf(snp_path, pop_path)
            fs = dadi.Spectrum.from_data_dict(dd, pop_ids=pop_ids, projections=proj, polarized=True)
            fs.to_file("../data/fs/{}_unfolded_projected.fs".format(snps))
    else:
        raise ValueError("Need to choose whether folded or unfolded spectra")

    # Printing out stats for the fs
    print("The datafile will be named {}".format(snps))
    print("\n* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n")
    print("Data for site frequency spectrum:\n")
    print("Sample sizes: {}".format(fs.sample_sizes))
    print("Sum of SFS: {}".format(np.around(fs.S(), 2)))
    if len(pop_ids) == 2:
        print("FST of SFS: {}".format(np.around(fs.Fst(), 2)))
        stat = np.around(fs.Fst(), 2)
    elif len(pop_ids) == 1:
        print("Tajima's D of SFS: {}".format(np.around(fs.Tajima_D())))
        stat = np.around(fs.Tajima_D(), 2)
    print("\n* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n")

    with open(stats_out_name, "a") as stats_out:
        stats_out.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(snps, fs.sample_sizes, fold, np.around(fs.S(), 2),
                                                           stat))

    # Plotting sfs and masked sfs
    # Not using rainbow palette anymore due to error
    colour_map = copy.copy(pylab.cm.get_cmap("hsv"))
    # within figure function, cmap=colour_map
    v_min = 0.05
    fig_size = (2.5, 2)
    fig1 = pylab.figure(figsize=fig_size)
    plt.rcParams.update({'font.size': 8})

    plot_out_name = args.snps + "_" + "-".join(map(str, proj))

    if len(pop_ids) == 1:
        if masked == "low":
            fs.mask[1] = True
            fs.mask[2] = True
            fig2 = pylab.figure(figsize=fig_size, cmap=colour_map)
            dadi.Plotting.plot_1d_fs(fs)
            fig2.tight_layout()
            fig2.savefig("../plots/spectra/" + plot_out_name + "_masked.png", dpi=300)
        elif masked == "no":
            dadi.Plotting.plot_1d_fs(fs)
            fig1.tight_layout()
            fig1.savefig("../plots/spectra/" + plot_out_name + ".png", dpi=300)
        else:
            print(":(")
    elif len(pop_ids) == 2:
        if masked == "low":
            fs.mask[0, 1] = True
            fs.mask[1, 0] = True
            fs.mask[2, 0] = True
            fs.mask[0, 2] = True
            fs.mask[1, 1] = True
            fig2 = pylab.figure(figsize=fig_size)
            dadi.Plotting.plot_single_2d_sfs(fs, vmin=v_min, cmap=colour_map)
            fig2.tight_layout()
            fig2.savefig("../plots/spectra/" + plot_out_name + "_2D_masked.png", dpi=300)
        elif masked == "no":
            dadi.Plotting.plot_single_2d_sfs(fs, vmin=v_min, cmap=colour_map)
            #dadi.Plotting.plot_single_2d_sfs(fs, vmin=v_min)
            fig1.tight_layout()
            fig1.savefig("../plots/spectra/" + plot_out_name + "_2D.png", dpi=300)
        else:
            raise ValueError("choose appropriate masking!")
    elif len(pop_ids) == 3:
        if masked == "low":
            print('Masking singletons and doubletons')
            fs.mask[0, 0, 1] = True
            fs.mask[0, 1, 0] = True
            fs.mask[1, 0, 0] = True
            fs.mask[0, 1, 1] = True
            fs.mask[1, 0, 1] = True
            fs.mask[1, 1, 0] = True
            fs.mask[2, 0, 0] = True
            fs.mask[0, 2, 0] = True
            fs.mask[0, 0, 2] = True
            dadi.Plotting.plot_3d_spectrum(fs, vmin=v_min, cmap=colour_map)
            fig1.tight_layout()
            fig1.savefig("../plots/spectra/" + plot_out_name + "_3D_masked.png", dpi=300)
        elif masked == "no":
            dadi.Plotting.plot_3d_spectrum(fs, vmin=v_min, cmap=colour_map)
            fig1.tight_layout()
            fig1.savefig("../plots/spectra/" + plot_out_name + "_3D.png", dpi=300)
        else:
            print(":(")
    else:
        raise ValueError("Pop IDs not configured appropriately - check input :)")

    if masked == "low":
        print("\n...Masking spectra...\n")
        fs.mask[0, 1] = True
        fs.mask[1, 0] = True
        fs.mask[2, 0] = True
        fs.mask[0, 2] = True
        fs.mask[1, 1] = True

        print("\n* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n")
        print("Data for masked site frequency spectrum:\n")
        print("Sample sizes: {}".format(fs.sample_sizes))
        print("Sum of SFS: {}".format(np.around(fs.S(), 2)))
        if len(pop_ids) == 2:
            print("FST of SFS: {}".format(np.around(fs.Fst(), 2)))
            stat = np.around(fs.Fst(), 2)
        elif len(pop_ids) == 1:
            print("Tajima's D of SFS: {}".format(np.around(fs.Tajima_D())))
            stat = np.around(fs.Tajima_D(), 2)
        print("\n* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n")

        with open(stats_out_name, "a") as stats_out:
            stats_out.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(snps, fs.sample_sizes, fold, np.around(fs.S(), 2),
                                                          stat))
    elif masked == "no":
        print("low frequency SNP masking not performed |-O-O-|")
    else:
        print("choose an appropriate option for masking - \"low\" \"mid\" or \"no\"")


if __name__ == '__main__':
    # Arguments
    parser = argparse.ArgumentParser(
        prog="dadi fs",
        description="A script for generating frequency spectra (fs) using dadi. Specify SNPs and other parameters.",
        usage="%(prog)s [options] <snps> <fold> <masked> <method> <genotypes>"
    )

    # Required positional arguments
    parser.add_argument(
        "snps",
        help="Name of the vcf file (e.g., AG1-AG2)."
    )
    parser.add_argument(
        "fold",
        help="Fold type (e.g., 'unfolded' or 'folded')."
    )
    parser.add_argument(
        "masked",
        help="Masking method (e.g., 'mid')."
    )
    parser.add_argument(
        "method",
        help="Where the fs is projected, subsampled or not (e.g., 'projection')."
    )
    parser.add_argument(
        "genotypes",
        type=int,
        nargs="+",
        help="List of integers specifying the genotype counts."
    )

    # Parse arguments
    args: Namespace = parser.parse_args()

    # Setting variables
    snps = args.snps
    fold = args.fold
    masked = args.masked
    method = args.method
    genotypes = args.genotypes

    # Call the main function with parsed arguments
    main(snps, fold, masked, method, genotypes)


