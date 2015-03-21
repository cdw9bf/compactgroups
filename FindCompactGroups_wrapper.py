"""
FindCompactGroups_wrapper.py
Use embarrassing parallelization to run FindCompactGroups several
times at once on multiple cores.
v1.0.0 Trey Wenger
"""
vers = "v1.0.0"

import sys
import numpy as np
import multiprocessing as mp
import FindCompactGroups
import time

def worker(snapnum):
    #print snapnum
    """
    The worker function. The argument is the thing that is iterated
    over. Everything else is hard-coded
    """
    filename = "data/datafile_{0}.txt".format(snapnum)
    bandwidth = 1.0
    neighborhood = 1.0
    min_members = 3
    max_sep_ratio = 1.0
    dwarf_range = 3.0
    velocity_filter = 1000.
    groupfile = "groups/groups_{0}.txt".format(snapnum)
    memberfile = "members/members_{0}.txt".format(snapnum)
    no_vel_filter = True
    use_dbscan = True
    verbose = False
    log_file = "logs/log_{0}.txt".format(snapnum)
    try:
        FindCompactGroups.main(filename,bandwidth=bandwidth,
                               neighborhood=neighborhood,
                               min_members=min_members,
                               max_sep_ratio=max_sep_ratio,
                               dwarf_range=dwarf_range,
                               velocity_filter=velocity_filter,
                               groupfile=groupfile,
                               memberfile=memberfile,
                               no_vel_filter=no_vel_filter,
                               use_dbscan=use_dbscan,
                               verbose=verbose,
                               log_file=log_file)
    except:
        print
        print "Problem in snapnum {0}".format(snapnum)
        print
    return 1

def main():
    """
    Main function - change what you need to for multiple calculations
    """
    start_time = time.time()
    print "Starting job at {0}".format(time.ctime())
    # get core count ( = cpus * 2)
    pool_size = mp.cpu_count() * 2
    # for example, let's process all 63 snapnums
    snapnums = np.arange(0,64)
    total = len(snapnums)
    # set up worker pool
    pool = mp.Pool(processes=pool_size)
    # run pool workers
    result = pool.map_async(worker,snapnums)
    while not result.ready():
        remaining = result._number_left
        done = total - remaining
        sys.stdout.write("\r[{0:20s}] {1}% Completed: {2} Remaining: {3}".\
                         format('#'*(20*done/total),100*done/total,
                         done,remaining))
        time.sleep(0.1)
    sys.stdout.write("\r[{0:20s}] {1}% Completed: {2} Remaining: {3}".\
                     format('#'*(20),100,
                     total,0))
    print
    pool.close()
    pool.join()
    end_time = time.time()
    print "Finished job at {0}".format(time.ctime())
    # calculate run-time
    time_diff = end_time - start_time
    hours = int(time_diff)/3600
    time_diff -= hours*3600
    mins = int(time_diff)/60
    time_diff -= mins*60
    secs = time_diff
    print "Runtime: {0}h {1}m {2:.2f}s".format(hours,mins,secs)

if __name__ == "__main__":
    main()
