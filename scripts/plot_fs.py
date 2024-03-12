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
import matplotlib as plt

# TODO: Make compatible for 2D
# TODO: Make various plotting functions

# Arguments
parser = argparse.ArgumentParser(prog="dadi Optimisation", usage='[options]')
parser.add_argument("snps")
parser.add_argument("dim")
args = parser.parse_args()

snps = args.snps
dim = args.dim
#snps = "out.unfold.folded.group1.saf.major"
#dim = "1D"

# def main()
# Import spectrum
fs = dadi.Spectrum.from_file('../data/fs/{}.fs'.format(snps))
pops = "{}".format(snps)
pop_ids = pops.split("-")

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
print("The datafile will be named {}".format(snps))
print("\n* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n")
print("Data for site frequency spectrum:\n")
print("Sample sizes: {}".format(fs.sample_sizes))
print("Sum of SFS: {}".format(np.around(fs.S(), 2)))
print("Tajima's D of SFS: {}".format(np.around(statistic, 2)))
print("Watterson's theta of SFS: {}".format(np.around(fs.Watterson_theta(), 2)))
print("\n* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n")

with open(stats_out_name, "a") as stats_out:
    stats_out.write("{0}\t{1}\t{2}\t{3}\n".format(snps, fs.sample_sizes, np.around(fs.S(), 2),
                                                  np.around(fs.Tajima_D(), 2),
                                                  np.around(fs.Watterson_theta(), 2)))

# Plotting
fig_size = (2.5, 2)
fig1 = pylab.figure(figsize=fig_size)
plt.rcParams.update({'font.size': 8})

if dim == "1D":
    dadi.Plotting.plot_1d_fs(fs)
elif dim == "2D":
    dadi.Plotting.plot_single_2d_sfs(fs)

fig1.tight_layout()
fig1.savefig("../plots/" + snps + dim + ".png", dpi=300)
