"""
FindCompactGroups.py
Search a file from the mini-millenium survey for compact groups
v1.0.0 - Trey Wenger - 20 January 2015
v1.0.1 - Trey Wenger - 04 February 2015
         Added progress indicators
         Re-calculate separation ratio after throwing out bad groups
v1.0.2 - Sophia Xiao - 04 February 2015
         Added galaxy velocity filter
v1.0.3 - TVW 20 March 2015
         added DBSCAN clustering algorithm, nearest neighbors,
         log_file
v1.0.4 - TVW 29 March 2015
         use pandas
"""
vers = "v1.0.3"

import sys
import argparse
import numpy as np
from sklearn.cluster import MeanShift,DBSCAN
from sklearn.neighbors import NearestNeighbors
import time
import pandas

#=====================================================================
# Compact Group Object
#=====================================================================
class CompactGroup:
    """
    Compact Group Object
    """
    def __init__(self, label, members):
        """
        Initialize ComactGroup Object
        """
        self.label = label
        self.members = members
        self.sep_ratio = None
        self.vel_median = None
        self.member_vel = None

    def eliminateDwarfs(self, min_range=3):
        """
        Eliminate all galaxies dimmer than min_range magnitudes in r
        of the brightest galaxy
        """
        # get minimum magnitude
        min_mag = np.min(self.members['mag_r'])
        # get index of all galaxies dimmer than min_range of min mag
        ind = [i for (i,mag) in enumerate(self.members['mag_r'])
               if mag > min_mag + min_range or mag == 99]     # mag=99 indicates galaxies without stars
        # delete them
        self.members = np.delete(self.members, ind, axis=0)

    def calculateMediod(self):
        """
        Calculate the mediod center of the cluster
        """
        x_med = np.median(self.members['x'])
        y_med = np.median(self.members['y'])
        z_med = np.median(self.members['z'])
        self.mediod = dict(x=x_med,y=y_med,z=z_med)

    def calculateRadius(self):
        """
        Calculate radius of cluster
        """
        dists = (self.members['x']-self.mediod["x"])**2.+\
                (self.members['y']-self.mediod["y"])**2.+\
                (self.members['z']-self.mediod["z"])**2.
        self.radius = np.sqrt(np.max(dists))

    def calculateVelocity(self):
        """
        Calculate median galaxy velocity of the group
        """
        velSquared = self.members['velX']**2.+self.members['velY']**2.+self.members['velZ']**2.
        self.member_vel = np.sqrt(velSquared)
        self.vel_median = np.median(self.member_vel)

    def velocityFilter(self, crit_vel=1000.0):
        """
        Eliminate galaxies moving crit_vel faster or slower than median velocity
        """
        ind = [i for (i,vel) in enumerate(self.member_vel)
               if abs(vel - self.vel_median) > crit_vel]
        # delete them
        self.members = np.delete(self.members, ind, axis=0)
        self.member_vel = np.delete(self.member_vel, ind, axis=0)

#=====================================================================
# Calculate Separation Ratio
#=====================================================================
def calcSepRatio(groups,verbose):
    """
    Calculate separation ratio for each group
    """
    # locations of groups
    locs = [[cg.mediod["x"],cg.mediod["y"],cg.mediod["z"]]
            for cg in groups]
    neighbors = NearestNeighbors(n_neighbors=2,algorithm='ball_tree')
    neighbors.fit(locs)
    for c_ind,cg in enumerate(groups):
        if verbose:
            sys.stdout.write("\r[{0:10s}] {1}%".\
                             format('#'*(10*(1+c_ind)/len(groups)),
                                    100*(1+c_ind)/len(groups)))
        closest = neighbors.kneighbors([cg.mediod["x"],
                                        cg.mediod["y"],
                                        cg.mediod["z"]])
        # neighbor[0] is me, neighbor[1] is actual neighbor
        dist = closest[0][0][1]
        cg.sep_ratio = cg.radius/dist

#=====================================================================
# Handle information
#=====================================================================
def logit(message, verbose, log_file):
    if verbose:
        print message
    if log_file is not None:
        with open(log_file,'a') as f:
            f.write(message+'\n')

#=====================================================================
# Find Groups Algorithm
#=====================================================================
def main(filename,bandwidth=0.1,neighborhood=0.5,
         min_members=3,max_sep_ratio=1.0,
         dwarf_range=3.0,
         velocity_filter=1000.0,
         groupfile='groups.txt',
         memberfile='members.txt',
         no_vel_filter=False,use_dbscan = False, verbose=True,
         log_file=None):
    """
    Search filename for compact groups
    """
    # erase log file
    if log_file is not None:
        with open(log_file,'w') as f:
            f.write('FindCompactGroups {0}\n'.format(vers))
    logit('Process started {0}'.format(time.ctime()),False,log_file)
    logit('Using data file {0}'.format(filename),False,log_file)
    logit('Using bandwidth {0}'.format(bandwidth),False,log_file)
    logit('Using neighborhood {0}'.format(neighborhood),False,log_file)
    logit('Using min_members {0}'.format(min_members),False,log_file)
    logit('Using max_sep_ratio {0}'.format(max_sep_ratio),False,log_file)
    logit('Using dwarf_range {0}'.format(dwarf_range),False,log_file)
    logit('Using velocity filter {0}'.format(velocity_filter),False,log_file)
    logit('Using groupfile {0}'.format(groupfile),False,log_file)
    logit('Using memberfile {0}'.format(memberfile),False,log_file)
    logit('no_vel_filter? {0}'.format(no_vel_filter),False,log_file)
    logit('use_dbscan? {0}'.format(use_dbscan),False,log_file)
    logit("Reading in data file...",verbose,log_file)
    # read in data file
    data = pandas.read_csv(filename,header=0,comment='#')
    logit("Found {0} galaxies in data file...".format(len(data)),
          verbose,log_file)
    # set up for clustering
    X = np.array(zip(data['x'],data['y'],data['z']))
    logit("Locating galaxy groups...",verbose,log_file)
    if not use_dbscan:
        logit("Using Mean-shift algorithm",verbose,log_file)
        # set up mean-shift algorithm
        ms = MeanShift(bandwidth=bandwidth,min_bin_freq=min_members,
                       cluster_all=False)
        # perform the fit
        ms.fit(X)
        labels = ms.labels_
    else:
        logit("Using DBSCAN algorithm",verbose,log_file)
        # set up DBSCAN
        db = DBSCAN(eps=neighborhood,min_samples=min_members)
        # perform the fit
        db.fit(X)
        labels = db.labels_
    # get unique labels
    labels_unique = np.unique(labels)
    n_clusters = len(labels_unique)

    # create cluster object for each cluster
    logit("Found {0} potential groups...".format(n_clusters),verbose,
          log_file)
    if n_clusters == 0:
        logit("No good groups!",verbose,log_file)
        logit("Done!",verbose,log_file)
        logit("Process completed {0}".format(time.ctime()),False,log_file)
        return
    logit("Measuring group properties...",verbose,log_file)
    groups = []
    for l_ind,label in enumerate(labels_unique):
        # skip label -1 (ungrouped)
        if label == -1:
            continue
        if verbose:
            sys.stdout.write("\r[{0:10s}] {1}%".\
                             format('#'*(10*(1+l_ind)/n_clusters),
                                    100*(1+l_ind)/n_clusters))
        # get all members of this group
        members_ind = [ind for (ind,this_label) in enumerate(labels)
                       if this_label == label]
        members = data[members_ind]
        # create object
        cg = CompactGroup(label,members)
        # eliminate dwarfs from group
        cg.eliminateDwarfs(min_range=dwarf_range)
        if len(cg.members) == 0: continue
        # calculate median velocity of group and velocity of galaxy
        cg.calculateVelocity()
        # eliminate passing-by galaxies from group
        if not no_vel_filter:
            cg.velocityFilter(crit_vel=velocity_filter)
        if len(cg.members) == 0: continue
        # calculate mediod of group
        cg.calculateMediod()
        # calculate radius of group
        cg.calculateRadius()
        # add it to the list
        groups.append(cg)
    if verbose:
        print

    if len(groups) == 0:
        logit("No good groups!",verbose,log_file)
        logit("Done!",verbose,log_file)
        logit("Process completed {0}".format(time.ctime()),False,log_file)
        return

    # calculate separation ratio for each group, and pull out
    # good groups
    if len(groups) == 1:
        logit("Only one group. Assigning separation ratio of nan",
              verbose,log_file)
        groups[0].sep_ratio = np.nan
        good_groups = groups
    else:
        logit("Calculating separation ratio of groups...",verbose,log_file)
        good_groups = []
        calcSepRatio(groups,verbose)
        if verbose:
            print
        logit("Eliminating groups with",verbose,log_file)
        logit("-- separation ratio > {0}".format(max_sep_ratio),verbose,
              log_file)
        logit("-- less than {0} members".format(min_members),verbose,
              log_file)
        for cg in groups:
            if cg.sep_ratio < max_sep_ratio and len(cg.members) >= min_members:
                good_groups.append(cg)
        logit("Found {0} good groups.".format(len(good_groups)),verbose,
              log_file)

    # re-calculate separation ratio for each group
    if len(good_groups) > 1:
        logit("Re-calculating separation ratio of good groups...",
              verbose,log_file)  
        calcSepRatio(good_groups,verbose)
        if verbose:
            print
    elif len(good_groups) == 0:
        logit("No good groups!",verbose,log_file)
        logit("Done!",verbose,log_file)
        logit("Process completed {0}".format(time.ctime()),False,log_file)
        return
    else:
        logit("Only one good group. Assigning separation ratio of nan",
              verbose,log_file)
        good_groups[0].sep_ratio = np.nan

    # save group statistics
    logit("Saving groups...",verbose,log_file)
    with open(groupfile,'w') as f:
        f.write("group_id,x,y,z,radius,vel_median,sep_ratio,num_members\n")
        for cg in good_groups:
            f.write("{0},{1},{2},{3},{4},{5},{6},{7}\n".\
                    format(cg.label,cg.mediod["x"],cg.mediod["y"],
                           cg.mediod["z"],cg.radius,cg.vel_median,
                           cg.sep_ratio,len(cg.members)))

    # save galaxy statistics
    logit("Saving group members...",verbose,log_file)
    with open(memberfile,'w') as f:
        f.write("group_id,member_id,x,y,z,velX,velY,velZ,velTot,mag_r\n")
        for cg in good_groups:
            vel_ind = 0
            for member in cg.members:
                f.write("{0},{1},{2},{3},{4},{5},{6},{7},{8},{9}\n".\
                        format(cg.label,member['galaxyID'],
                               member['x'],member['y'],
                               member['z'],member['velX'],
                               member['velY'],member['velZ'],cg.member_vel[vel_ind],member['mag_r']))
                vel_ind += 1
    logit("Done!",verbose,log_file)
    logit("Process completed {0}".format(time.ctime()),False,log_file)

#=====================================================================
# Command Line Setup
#=====================================================================
if __name__ == "__main__":
    """
    Set up shell arguments
    """
    # set up argument parser
    parser = argparse.ArgumentParser(
        description="Find Compact Groups in Millenium Simulation",
        prog='FindCompactGroups.py',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--version', action='version',
                        version='%(prog)s '+vers)
    # add required arguments
    required=parser.add_argument_group('required arguments')
    required.add_argument('inputfile',type=str,
                          help="input galaxy file")
    # add optional arguments
    semi_opt = parser.add_argument_group('arguments set to defaults')
    semi_opt.add_argument('--bandwidth',type=float,
                          help='bandwidth of mean shift algorithm',
                          default=0.1)
    semi_opt.add_argument('--neighborhood',type=float,
                          help='neighborhood of DBSCAN algorithm (max distance between two galaxies to be considered part of same group)',
                          default=0.5)
    semi_opt.add_argument('--min_members',type=int,
                          help='minimum number of members in group',
                          default=3)
    semi_opt.add_argument('--dwarf_range',type=float,
                          help='galaxies dimmer than dwarf_range of the brightest galaxy (mag_r) will be excluded',
                          default=3.0)
    semi_opt.add_argument('--velocity_filter',type=float,
                          help='galaxies moving faster than velocity_filter will be excluded',
                          default=1000.0)
    semi_opt.add_argument('--max_sep_ratio',type=float,
                          help='maximum separation ratio to be considered a group',
                          default=1.0)
    semi_opt.add_argument('--groupfile',type=str,
                          help='output file with group statistics',
                          default='groups.txt')
    semi_opt.add_argument('--memberfile',type=str,
                          help='output file with group member statistics',
                          default='members.txt')
    semi_opt.add_argument('--no_vel_filter',action='store_true',
                          help='suppress velocity filter',
                          default=False)
    semi_opt.add_argument('--use_dbscan',action='store_true',
                          help='Use DBSCAN clustering algorithm (default is mean shift)',
                          default=False)
    semi_opt.add_argument('--verbose',action='store_true',
                          help='Print out status along the way',
                          default=False)
    semi_opt.add_argument('--log_file',type=str,
                          help='Print out status to log file of this name',
                          default=None)

    # get command line arguments
    args = vars(parser.parse_args())
    # run with arguments
    main(args['inputfile'],bandwidth=args['bandwidth'],
         neighborhood=args['neighborhood'],
         min_members=args['min_members'],
         dwarf_range=args['dwarf_range'],
         velocity_filter=args['velocity_filter'],
         max_sep_ratio=args['max_sep_ratio'],
         groupfile=args['groupfile'],
         memberfile=args['memberfile'],
         no_vel_filter=args['no_vel_filter'],
         use_dbscan=args['use_dbscan'],
         verbose=args['verbose'],
         log_file=args['log_file'])
