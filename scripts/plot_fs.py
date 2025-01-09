#!/usr/anaconda3/bin/python
# -*- coding: utf-8 -*-
# remember to edit python environment if required.

"""
@author: kprata
@date created: 25/1/24
@description: Plotting the fs spectrum for a premade .fs file
"""

import argparse
import os
import dadi
import numpy as np
import copy
import pylab
import matplotlib.pyplot as pyplot


def main(snps, fold, method, mask):
    # Import spectrum
    fs = dadi.Spectrum.from_file('../data/fs/{}.fs'.format(snps))
    pops = "{}".format(snps)
    pop_ids = pops.split("-")

    if mask == "low" and len(fs.sample_sizes) == 2:
        fs.mask[1, 0] = True
        fs.mask[0, 1] = True
        fs.mask[2, 0] = True
        fs.mask[0, 2] = True
        fs.mask[1, 1] = True
    elif mask == "mid" and len(fs.sample_sizes) == 1:
        mid = fs.sample_sizes / 2
        mid = int(mid[0])
        fs.mask[mid] = True
    else:
        print("No masking")

    if len(fs.sample_sizes) == 2:
        statistic_name = "FST"
        statistic = fs.Fst()
    elif len(fs.sample_sizes) == 1:
        statistic_name = "TajimasD"
        statistic = fs.Tajima_D()
    else:
        print("Choose dim")

    # make a file with statistics about your sfs
    stats_out_name = "../results/sfs_stats-2.txt"
    with open(stats_out_name, 'a') as stats_out:
        if stats_out.tell() == 0:
            print('Creating a new file')
            stats_out.write("Pop\tSample sizes\tSum of SFS\t{}\n".format(statistic_name))
        else:
            print('file exists, appending')

    # Printing out stats for the fs
    print("\nThe datafile will be named {}".format(snps))
    print("\n* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n")
    print("Data for site frequency spectrum:\n")
    print("Sample sizes: {}".format(fs.sample_sizes))
    print("Sum of SFS: {}".format(np.around(fs.S(), 2)))
    print("{} of SFS: {}".format(statistic_name, np.around(statistic, 2)))
    print("\n* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n")

    with open(stats_out_name, "a") as stats_out:
        stats_out.write("{0}\t{1}\t{2}\t{3}\n".format(snps, fs.sample_sizes, np.around(fs.S(), 2),
                                                      np.around(statistic, 2)))

    # Plotting settings
    fig_size = (5, 5)
    fig1 = pylab.figure(figsize=fig_size)
    colour_map = copy.copy(pylab.cm.get_cmap("hsv"))

    extras = ""

    # Folding
    if fold == "folded":
        fs = fs.fold()
        extras += "fold"
        fs.to_file("../data/fs/{}_folded.fs".format(snps))

    if mask != "no":
        extras += "mask{}".format(mask)

    print("Plotting spetrum for {}".format(snps))

    if len(fs.sample_sizes) == 1:
        if method == "projection":
            project = fs.sample_sizes * 0.8
            project_int = int(project)
            fs = fs.project([project_int])
            extras += "project"
            fs.to_file("../data/fs/{}_projected0.8{}.fs".format(snps, extras))
        dadi.Plotting.plot_1d_fs(fs)
    elif len(fs.sample_sizes == 2):
        # Projecting
        if method == "projection":
            project = fs.sample_sizes * 0.8
            project_int = [int(x) for x in project]
            fs = fs.project([project_int[0], project_int[1]])
            extras += "project"
            fs.to_file("../data/fs/{}_projected0.8{}.fs".format(snps, extras))
        dadi.Plotting.plot_single_2d_sfs(fs, vmin=1, cmap=colour_map)

    fig1.tight_layout()
    fig1.savefig("../plots/spectra/" + snps + extras + ".png", dpi=300)


if __name__ == '__main__':
    # Arguments
    parser = argparse.ArgumentParser(
        prog="Plotting fs",
        description="A script for plotting fs from .fs file.",
        usage="%(prog)s [options] <snps> <fold> <method> <mask>"
    )

    # Required positional arguments
    parser.add_argument(
        "snps",
        help="Name of the data .fs file."
    )
    parser.add_argument(
        "fold",
        help="If you want your spectrum to be folded (e.g., 'folded'). "
             "Only use on unfolded spectrum!"
    )
    parser.add_argument(
        "method",
        help="Method to use for projection or subsample (e.g., 'projection, "
             "'subsample' or 'no')"
    )
    parser.add_argument(
        "mask",
        help="Type of masking to use (e.g., 'mid', 'low', or 'no')."
    )

    # Parse arguments
    args = parser.parse_args()

    # Set variables from parsed arguments
    snps = args.snps
    fold = args.fold
    method = args.method
    mask = args.mask

    # Call the main function
    main(snps, fold, method, mask)
