#!/usr/anaconda3/bin/python
# -*- coding: utf-8 -*-
# remember to edit python environment if required.

"""
@author: kprata
@date created: 25/1/24
@description: TODO
"""

import argparse
import os
import dadi
import numpy as np
import copy
import pylab
import matplotlib.pyplot as pyplot

# Arguments
parser = argparse.ArgumentParser(prog="dadi Optimisation", usage='[options]')
parser.add_argument("snps")
parser.add_argument("dim")
parser.add_argument("fold")
parser.add_argument("proj")
parser.add_argument("mask")
args = parser.parse_args()

snps = args.snps
dim = args.dim
fold = args.fold
proj = args.proj
mask = args.mask
# snps = "ahya.unfold.het05.group1.group2_snp6966117.2d"
# dim = "2D"
# os.chdir("scripts")
# proj = "yes" # projecting not working??

# def main()
# Import spectrum
fs = dadi.Spectrum.from_file('../data/fs/{}.fs'.format(snps))
pops = "{}".format(snps)
pop_ids = pops.split("-")

if mask == "low" and dim == "2D":
    fs.mask[1, 0] = True
    fs.mask[0, 1] = True
    fs.mask[2, 0] = True
    fs.mask[0, 2] = True
    fs.mask[1, 1] = True
elif mask == "mid" and dim == "1D":
    mid = fs.sample_sizes / 2
    mid = int(mid[0])
    fs.mask[mid] = True
else:
    print("No masking")

if dim == "2D":
    statistic_name = "FST"
    statistic = fs.Fst()
elif dim == "1D":
    statistic_name = "TajimasD"
    statistic = fs.Tajima_D()
else:
    print("Choose dim")

# make a file with statistics about your sfs
stats_out_name = "../results/sfs_stats_ahya.txt"
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

# Plotting
fig_size = (5, 5)
fig1 = pylab.figure(figsize=fig_size)
colour_map = copy.copy(pylab.cm.get_cmap("hsv"))

extras = ""

# Folding
if fold == "yes":
    fs = fs.fold()
    extras += "fold"
    fs.to_file("../data/fs/{}_folded.fs".format(snps))

if mask != "no":
    extras += "mask{}".format(mask)

if dim == "1D":
    if proj == "yes":
        project = fs.sample_sizes * 0.8
        project_int = int(project)
        fs = fs.project([project_int])
        extras += "project"
        fs.to_file("../data/fs/{}_projected.fs".format(snps))
    dadi.Plotting.plot_1d_fs(fs)
elif dim == "2D":
    # Projecting
    if proj == "yes":
        project = fs.sample_sizes * 0.8
        project_int = [int(x) for x in project]
        fs = fs.project([project_int[0], project_int[1]])
        extras += "project"
        fs.to_file("../data/fs/{}_projected.fs".format(snps))
    dadi.Plotting.plot_single_2d_sfs(fs, vmin=1, cmap=colour_map)

fig1.tight_layout()
fig1.savefig("../plots/" + snps + dim + extras + ".png", dpi=300)
