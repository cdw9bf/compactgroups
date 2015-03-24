"""
makeCGplots.py
Make informative plots for compact groups and galaxy members
v1.0.0 - Sophia Xiao - 6 February 2015
         created "single" mode: make 1d histograms of cg & galaxies parameters
v1.0.1 - Sophia Xiao - 17 February 2015
         added in "multi" mode: make money plots:
                                (1) cg number evolution in z
                                (2) and normalized galaxy number in z
         added in "filter" mode: make 2d histograms for parameters in z
"""
vers = "v1.0.1"

# =============================================================================
#                    Import Modules and Parsers
# =============================================================================

import sys, glob, re
import argparse
import numpy as np
import matplotlib
matplotlib.use("agg")
matplotlib.rcParams.update({'font.size': 16, 'text.usetex': True})
import matplotlib.pyplot as plt

# Define parser variables
parser = argparse.ArgumentParser(
    description="Make Plots For Compact Groups And Members of Compact Groups",
    prog='makeCGplots.py',
    formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument('--version', '-v', action='version', version='%(prog)s '+vers)

parser.add_argument('--mode',type=str,
                     help="choose plotting mode [MODE: single OR multi OR filter]\n"
                           "- single: use 1 group file and 1 member file to make histograms\n"
                           "- multi: use 1 directory of group files in different redshifts\n"
                           "to make money plots\n"
                           "- filter: use three num_gp files for 2d histograms plotting,\n"
                           "showing the effects of different constraints on the number of compact groups",
                     default="single")

parser.add_argument('--groupfile',type=str,metavar="FILE",
                     help="(mode: single) input file with group statistics",
                     default='groups.txt')

parser.add_argument('--memberfile',type=str,metavar="FILE",
                     help="(mode: single) input file with galaxy statistics",
                     default='members.txt')

parser.add_argument('--dirname',type=str,metavar="DIR",
                     help="(mode: multi) specify directory of redshift runs, ending with '/'",
                     default="cgdatafiles/")

parser.add_argument('--vel_filter_file',type=str,metavar="FILE",
                     help="(mode: filter) input file with velocity filter statistics",
                     default='vel_filter_snapnum.txt')

parser.add_argument('--bandwidth_file',type=str,metavar="FILE",
                     help="(mode: filter) input file with bandwidth statistics",
                     default='bandwidth_snapnum.txt')

parser.add_argument('--sep_ratio_file',type=str,metavar="FILE",
                     help="(mode: filter) input file with separation ratio statistics",
                     default='sep_ratio_snapnum.txt')

args = parser.parse_args()

# =============================================================================
#                         Single Run Plotting Mode
# =============================================================================

if args.mode == "single":
    # read in data file
    print "MODE: single\n"
    print "Reading in data files..."
    groups = np.genfromtxt(args.groupfile, dtype=None, delimiter=',',
                          names=True, comments="#")
    members = np.genfromtxt(args.memberfile, dtype=None, delimiter=',',
                            names=True, comments="#")
    print "Found {0} groups in data file...".format(len(groups))
    print "Fount {0} galaxies in data file...".format(len(members))

    # calculate velocity difference from group median
    veldiff = []
    for i in range(len(groups['group_id'])):
       for j in range(len(members['group_id'])):
            if groups['group_id'][i] == members['group_id'][j]:
                veldiff.append(members['velTot'][j] - groups['vel_median'][i])

    # make 1d histograms
    print "Writing 1-D histograms..."
    fig = plt.figure()

    # histograms for groups.txt
    # cg median velocity histogram
    ax = fig.gca()
    plt.hist(groups['vel_median'], bins=100, color='#19198D')
    ax.grid()
    ax.set_xlabel("$v_{med}$ of CGs [$km/s$]")
    ax.set_ylabel("Number of CGs")
    fig.savefig("hist1d_cg_medianVelocity.png")
    fig.clf()

    # velocity dispersion plot
    ax = fig.gca()
    plt.hist(veldiff, bins=100, color='#19198D')
    ax.grid()
    ax.set_xlabel("$v-v_{med}$ [$km/s$]")
    ax.set_ylabel("Number of CGs")
    fig.savefig("hist1d_cg_velDispersion.png")
    fig.clf()

    # radius histogram
    ax = fig.gca()
    plt.hist(groups['radius'], bins=100, color='#19198D')
    ax.grid()
    ax.set_xlabel("CG Radius [Mpc/h]")
    ax.set_ylabel("Number of CGs")
    fig.savefig("hist1d_cg_radius.png")
    fig.clf()

    # num_members histogram
    ax = fig.gca()
    plt.hist(groups['num_members'], bins=len(set(groups['num_members']))-1, color='#19198D')
    ax.grid()
    ax.set_xlabel("Number of Galaxy Members")
    ax.set_ylabel("Number of CGs")
    fig.savefig("hist1d_cg_numGalaxy.png")
    fig.clf()

    # histograms for members.txt
    # magnitude histogram
    ax = fig.gca()
    plt.hist(members['mag_r'], bins=100, color='#19198D')
    ax.grid()
    ax.set_xlabel("R Band Magnitude")
    ax.set_ylabel("Number of Galaxies")
    fig.savefig("hist1d_galaxy_mag_r.png")
    fig.clf()

    # galaxy velocity histogram
    ax = fig.gca()
    plt.hist(members['velTot'], bins=100, color='#19198D')
    ax.grid()
    ax.set_xlabel("$v$ [$km/s$]")
    ax.set_ylabel("Number of Galaxies")
    fig.savefig("hist1d_galaxy_vel.png")
    fig.clf()

# =============================================================================
#                         Multi Runs Plotting Mode
# =============================================================================

if args.mode == "multi":
    # read in files
    print "MODE: multi\n"
    snapnum = []
    num_gp = []
    num_mem = []
    num_totGal = []

    print "Reading in data files..."
    # group files
    for files in glob.glob(args.dirname+'groups_*.txt'):
        groups = np.genfromtxt(files, dtype=None, delimiter=',',
                               names=True, comments="#")
        num_gp.append(len(groups))

    # member files
    for files in glob.glob(args.dirname+'members_*.txt'):
        members = np.genfromtxt(files, dtype=None, delimiter=',',
                                names=True, comments="#")
        num_mem.append(len(members))
        snapnum.append(re.sub(args.dirname+'|members_|.txt','',files))

    # data files
    datalist = glob.glob(args.dirname+'datafile_*.txt')
    for ind, files in enumerate(datalist):
        sys.stdout.write("\r[{0:10s}] {1}%".\
                         format('#'*(10*(1+ind)/len(datalist)), 100*(1+ind)/len(datalist)))
        sys.stdout.flush()
        z1 = re.sub(args.dirname+'|datafile_|.txt','',files)
        for z2 in snapnum:
            if z1 == z2:
                data = np.genfromtxt(files, dtype=None, delimiter=',',
                                     names=True, comments="#")
                num_totGal.append(len(data))
    print
    num_totGal = [float(i) for i in num_totGal]

    # make money plot
    print "Writing money plot..."
    fig = plt.figure()

    ax = fig.gca()
    ax.plot(snapnum,num_gp,'kx')
    ax.grid()
    ax.set_xlabel("Redshift [snapnum]")
    ax.set_ylabel("Number of CGs")
    ax.set_title("Number of Compact Groups at Different Redshifts")
    fig.savefig("moneyplot.png")
    fig.clf()

    ax = fig.gca()
    ax.plot(snapnum,np.divide(num_mem,num_totGal),'kx')
    ax.grid()
    ax.set_xlabel("Redshift [snapnum]")
    ax.set_ylabel("$\mathrm{Num}_{\mathrm{mem}}/\mathrm{Num}_{\mathrm{galaxies}}$")
    ax.set_title("Fraction of Galaxies at Different Redshifts")
    fig.savefig("galaxyFraction.png")
    fig.clf()

# =============================================================================
#                         Filter Plotting Mode
# =============================================================================

if args.mode == "filter":
    print "MODE: filter\n"
    print "Reading in data files..."
    velfile = np.genfromtxt(args.vel_filter_file, dtype=None, delimiter=',',
                               names=True, comments="#")
    bwfile = np.genfromtxt(args.bandwidth_file, dtype=None, delimiter=',',
                               names=True, comments="#")
    srfile = np.genfromtxt(args.sep_ratio_file, dtype=None, delimiter=',',
                               names=True, comments="#")

    print "Writing 2-D histograms..."
    fig = plt.figure()

    # velocity filter 2d histogram
    ax = fig.gca()
    coll = ax.scatter(velfile['vel_cut'], velfile['snapnum'], c=velfile['num_gp'], lw=0)
    ax.grid()
    ax.set_xlabel("Velocity Filter [$km/s$]")
    ax.set_ylabel("Redshift [snapnum]")
    fig.colorbar(coll, ax=ax).set_label("Number of Compact Groups")
    fig.savefig("hist2d_velCut_z.png")
    fig.clf()

    # bandwidth 2d histogram
    ax = fig.gca()
    coll = ax.scatter(bwfile['bandwidth'], bwfile['snapnum'], c=bwfile['num_gp'], lw=0)
    ax.grid()
    ax.set_xlabel("Bandwidth")
    ax.set_ylabel("Redshift [snapnum]")
    fig.colorbar(coll, ax=ax).set_label("Number of Compact Groups")
    fig.savefig("hist2d_bandwidth_z.png")
    fig.clf()

    # separation ratio 2d histogram
    ax = fig.gca()
    coll = ax.scatter(srfile['sep_ratio'], srfile['snapnum'], c=srfile['num_gp'], lw=0)
    ax.grid()
    ax.set_xlabel("Separation Ratio")
    ax.set_ylabel("Redshift [snapnum]")
    fig.colorbar(coll, ax=ax).set_label("Number of Compact Groups")
    fig.savefig("hist2d_sepRatio_z.png")
    fig.clf()
