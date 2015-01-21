"""
FindCompactGroups.py
Search a file from the mini-millenium survey for compact groups
v1.0.0 - Trey Wenger - 20 January 2015
"""
vers = "v1.0.0"

import sys
import argparse
import numpy as np
from sklearn.cluster import MeanShift

#=====================================================================
# Compact Group Object
#=====================================================================
class CompactGroup:
    """
    Compact Group Object
    """
    def __init__(self, label,members):
        """
        Initialize ComactGroup Object
        """
        self.label = label
        self.members = members

    def eliminateDwarfs(self, min_range=3):
        """
        Eliminate all galaxies dimmer than min_range magnitudes in r
        of the brightest galaxy.
        """
        # get minimum magnitude
        min_mag = np.min(self.members['mag_r'])
        # get index of all galaxies dimmer than min_range of min mag
        ind = [i for (i,mag) in enumerate(self.members['mag_r'])
               if mag > min_mag + min_range]
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
        dists = np.sqrt((self.members['x']-self.mediod["x"])**2.+
                        (self.members['y']-self.mediod["y"])**2.+
                        (self.members['z']-self.mediod["z"])**2.)
        self.radius = np.max(dists)

    def calculateSeparationRatio(self,groups):
        """
        Calculate the separation ratio between this and other groups
        """
        dists = np.array([np.sqrt((cg.mediod["x"] - self.mediod["x"])**2.+
                                  (cg.mediod["y"] - self.mediod["y"])**2.+
                                  (cg.mediod["z"] - self.mediod["z"])**2.)
                          for cg in groups if cg != self])
        min_sep = np.min(dists)
        self.sep_ratio = self.radius/min_sep

#=====================================================================
# Find Groups Algorithm
#=====================================================================
def main(filename,bandwidth=0.1,min_members=3,max_sep_ratio=1.0,
         dwarf_range=3.0,
         groupfile='groups.txt',
         memberfile='members.txt'):
    """
    Search filename for compact groups
    """
    # read in data file
    data = np.genfromtxt(filename,dtype=None,delimiter=',',
                         names=True,comments="#")
    # set up mean-shift algorithm
    ms = MeanShift(bandwidth=bandwidth,min_bin_freq=min_members,
                   cluster_all=False)
    # perform the fit
    X = np.array(zip(data['x'],data['y'],data['z']))
    ms.fit(X)
    labels = ms.labels_
    # get unique labels
    labels_unique = np.unique(labels)
    n_clusters = len(labels_unique)

    # create cluster object for each cluster
    groups = []
    for label in labels_unique:
        # skip label -1 (ungrouped)
        if label == -1:
            continue
        # get all members of this group
        members_ind = [ind for (ind,this_label) in enumerate(labels)
                       if this_label == label]
        members = data[members_ind]
        # create object
        cg = CompactGroup(label,members)
        # eliminate dwarfs from group
        cg.eliminateDwarfs(min_range=dwarf_range)
        # calculate mediod of group
        cg.calculateMediod()
        # calculate radius of group
        cg.calculateRadius()
        # add it to the list
        groups.append(cg)

    # calculate separation ratio for each group, and pull out
    # good groups
    good_groups = []
    for cg in groups:
        cg.calculateSeparationRatio(groups)
        if cg.sep_ratio < max_sep_ratio and len(cg.members) > min_members:
            good_groups.append(cg)

    print "Found {0} good groups.".format(len(good_groups))

    # save group statistics
    with open(groupfile,'w') as f:
        f.write("group_id,x,y,z,radius,sep_ratio,num_members\n")
        for cg in good_groups:
            f.write("{0},{1},{2},{3},{4},{5},{6}\n".\
                    format(cg.label,cg.mediod["x"],cg.mediod["y"],
                           cg.mediod["z"],cg.radius,cg.sep_ratio,
                           len(cg.members)))

    # save galaxy statistics
    with open(memberfile,'w') as f:
        f.write("group_id,member_id,x,y,z,mag_r\n")
        for cg in good_groups:
            for member in cg.members:
                f.write("{0},{1},{2},{3},{4},{5}\n".\
                        format(cg.label,member['galaxyID'],
                               member['x'],member['y'],
                               member['z'],member['mag_r']))

#=====================================================================
# Command Line Setup
#=====================================================================
if __name__ == "__main__":
    """
    Set up shell arguments
    """
    # set up argument parser
    parser=argparse.ArgumentParser(
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
    semi_opt=parser.add_argument_group('arguments set to defaults')
    semi_opt.add_argument('--bandwidth',type=float,
                          help='bandwidth of mean shift algorithm',
                          default=0.1)
    semi_opt.add_argument('--min_members',type=int,
                          help='minimum number of members in group',
                          default=3)
    semi_opt.add_argument('--dwarf_range',type=float,
                          help='galaxies dimmer than dwarf_range of the brightest galaxy (mag_r) will be excluded',
                          default=3.0)
    semi_opt.add_argument('--max_sep_ratio',type=float,
                          help='maximum separation ratio to be considered a group',
                          default=1.0)
    semi_opt.add_argument('--groupfile',type=str,
                          help='output file with group statistics',
                          default='groups.txt')
    semi_opt.add_argument('--memberfile',type=str,
                          help='output file with group member statistics',
                          default='members.txt')
    # get command line arguments
    args = vars(parser.parse_args())
    # run with arguments
    main(args['inputfile'],bandwidth=args['bandwidth'],
         min_members=args['min_members'],
         dwarf_range=args['dwarf_range'],
         max_sep_ratio=args['max_sep_ratio'],
         groupfile=args['groupfile'],
         memberfile=args['memberfile'])
