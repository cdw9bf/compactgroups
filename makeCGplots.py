"""
makeCGplots.py
Make informative plots for compact groups and galaxy members
v1.0.0 - Sophia Xiao - 6 February 2015
"""
vers = "v1.0.0"

# =============================================================================
#                    Import Modules and Parsers
# =============================================================================

import sys
import argparse
import numpy as np
import matplotlib
matplotlib.use("agg")
matplotlib.rcParams.update({'font.size': 16, 'text.usetex': True})
from matplotlib import pyplot

# Define parser variables
parser = argparse.ArgumentParser(
    description="Make Plots For Compact Groups And Members of Compact Groups",
    prog='makeCGplots.py',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('--version', action='version', version='%(prog)s '+vers)

parser.add_argument('--groupfile',type=str,
                     help="input file with group statistics",
                     default='groups.txt')

parser.add_argument('--memberfile',type=str,
                     help="input file with group member statistics",
                     default='members.txt')

args = parser.parse_args()

# =============================================================================
#                                Main
# =============================================================================

# read in data file
print "Reading in data file..."
groups = np.genfromtxt(args.groupfile, dtype=None, delimiter=',',
                       names=True, comments="#")
members = np.genfromtxt(args.memberfile, dtype=None, delimiter=',',
                       names=True, comments="#")
print "Found {0} groups in data file...".format(len(groups))
print "Fount {0} galaxies in data file...".format(len(members))

veldiff = []
for i in range(len(groups['group_id'])):
    for j in range(len(members['group_id'])):
        if groups['group_id'][i] == members['group_id'][j]:
            veldiff.append(members['velTot'][j] - groups['vel_median'][i])

# =============================================================================
#                          Make informative plots
# =============================================================================
print "Writing graphs..."
fig = pyplot.figure()

ax = fig.gca()
pyplot.hist(veldiff, bins=100, color='#19198D')
ax.grid()
ax.set_xlabel("$v-v_{med}$ ($km/s$)")
ax.set_ylabel("Number of Galaxies")
fig.savefig("hist1d_velocity.png")
fig.clf()

ax = fig.gca()
pyplot.hist(members['mag_r'], bins=100, color='#19198D')
ax.grid()
ax.set_xlabel("R Band Magnitude")
ax.set_ylabel("Number of Galaxies")
fig.savefig("hist1d_mag_r.png")
fig.clf()
