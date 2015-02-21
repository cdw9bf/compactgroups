"""
FindCompactGroups.py
Search a file from the mini-millenium survey for compact groups
v1.0.0 - Trey Wenger - 20 January 2015
v1.0.1 - Trey Wenger - 04 February 2015
         Added progress indicators
         Re-calculate separation ratio after throwing out bad groups
v1.0.2 - Sophia Xiao - 04 February 2015
         Added galaxy velocity filter
"""
vers = "v1.0.2"

import sys
import argparse
import numpy as np
from sklearn.cluster import MeanShift
from sklearn.neighbors import NearestNeighbors

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
        of the brightest galaxy.
        """
        # get minimum magnitude
        min_mag = np.min(self.members['mag_r'])
        # get index of all galaxies dimmer than min_range of min mag
        ind = [i for (i,mag) in enumerate(self.members['mag_r'])
<<<<<<< HEAD
               if (mag > min_mag + min_range) or mag == 99]
=======
               if mag > min_mag + min_range]
>>>>>>> FETCH_HEAD
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
def calcSepRatio(groups):
    """
    Calculate separation ratio for each group
    """
    # locations of groups
    locs = [[cg.mediod["x"],cg.mediod["y"],cg.mediod["z"]]
            for cg in groups]
    neighbors = NearestNeighbors(n_neighbors=2,algorithm='ball_tree')
    neighbors.fit(locs)
    for c_ind,cg in enumerate(groups):
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
# Find Groups Algorithm
#=====================================================================
def main(filename,bandwidth=0.1,min_members=3,max_sep_ratio=1.0,
         dwarf_range=3.0,
         velocity_filter=1000.0,
         groupfile='groups.txt',
         memberfile='members.txt'):
    """
    Search filename for compact groups
    """
    print "Reading in data file..."
    # read in data file
    data = np.genfromtxt(filename,dtype=None,delimiter=',',
                         names=True,comments="#")
    print "Found {0} galaxies in data file...".format(len(data))
    # set up mean-shift algorithm
    ms = MeanShift(bandwidth=bandwidth,min_bin_freq=min_members,
                   cluster_all=False)
    # perform the fit
    X = np.array(zip(data['x'],data['y'],data['z']))
    print "Locating galaxy groups..."
    ms.fit(X)
    labels = ms.labels_
    # get unique labels
    labels_unique = np.unique(labels)
    n_clusters = len(labels_unique)

    # create cluster object for each cluster
    print "Found {0} potential groups...".format(n_clusters)
    print "Measuring group properties..."
    groups = []
    for l_ind,label in enumerate(labels_unique):
        # skip label -1 (ungrouped)
        if label == -1:
            continue
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
<<<<<<< HEAD
        if len(cg.members) == 0: continue
=======
>>>>>>> FETCH_HEAD
        # calculate median velocity of group and velocity of galaxy
        cg.calculateVelocity()
        # eliminate passing-by galaxies from group
        cg.velocityFilter(crit_vel=velocity_filter)
        if len(cg.members) == 0: continue
        # calculate mediod of group
        cg.calculateMediod()
        # calculate radius of group
        cg.calculateRadius()
        # add it to the list
        groups.append(cg)
    print

    # calculate separation ratio for each group, and pull out
    # good groups
    print "Calculating separation ratio of groups..."
    good_groups = []
    calcSepRatio(groups)
    print
    print "Eliminating groups with"
    print "-- separation ratio > {0}".format(max_sep_ratio)
    print "-- less than {0} members".format(min_members)
    for cg in groups:
        if cg.sep_ratio < max_sep_ratio and len(cg.members) >= min_members:
            good_groups.append(cg)

    print "Found {0} good groups.".format(len(good_groups))

    # re-calculate separation ratio for each group
    if len(good_groups) > 0:
        print "Re-calculating separation ratio of good groups..."
        calcSepRatio(good_groups)
        print
    else:
        print "No good groups!"
        return

    # save group statistics
    print "Saving groups..."
    with open(groupfile,'w') as f:
        f.write("group_id,x,y,z,radius,vel_median,sep_ratio,num_members\n")
        for cg in good_groups:
            f.write("{0},{1},{2},{3},{4},{5},{6},{7}\n".\
                    format(cg.label,cg.mediod["x"],cg.mediod["y"],
                           cg.mediod["z"],cg.radius,cg.vel_median,
                           cg.sep_ratio,len(cg.members)))

    # save galaxy statistics
    print "Saving group members..."
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
    # get command line arguments
    args = vars(parser.parse_args())
    # run with arguments
    main(args['inputfile'],bandwidth=args['bandwidth'],
         min_members=args['min_members'],
         dwarf_range=args['dwarf_range'],
         velocity_filter=args['velocity_filter'],
         max_sep_ratio=args['max_sep_ratio'],
         groupfile=args['groupfile'],
         memberfile=args['memberfile'])
