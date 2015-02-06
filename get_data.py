"""
get_data.py
Download mini-millennium simulation data
v1.0.0 - Trey Wenger - 05 February 2015
"""
vers = "v1.0.0"

import sys
import argparse
import requests

def main(outputfile,redshift_min=None,redshift_max=None,snapnum=None,
         columns='galaxyID,x,y,z,redshift',limit=None):
    """
    Download data from Millennium Simulation
    """
    # data access url
    url = "http://gavo.mpa-garching.mpg.de/Millennium/?action=doQuery&SQL="
    # set up sql query
    query = "select "
    if limit is not None:
        query += "top {0} ".format(limit)
    query += columns
    query += " from millimil..DeLucia2006a "
    if snapnum is not None:
        query += "where snapnum={0}".format(snapnum)
    elif redshift_max is not None and redshift_min is not None:
        query += "where redshift between {0} and {1}".format(redshift_min,
                                                             redshift_max)
    elif redshift_min is not None:
        query += "where redshift > {0}".format(redshift_min)
    elif redshift_max is not None:
        query += "where redshift < {0}".format(redshift_max)
    print "Executing query:"
    print query
    page = requests.get(url+query)
    lines = page.content.split('\n')
    print "Saving file as {0}".format(outputfile)
    with open(outputfile,'w') as f:
        for line in lines:
            if not line.startswith('#'):
                f.write(line+'\n')
    print "Done!"

#=====================================================================
# Command Line Setup
#=====================================================================
if __name__ == "__main__":
    """
    Set up shell arguments
    """
    # set up argument parser
    parser=argparse.ArgumentParser(
        description="Download data from Millennium Simulation",
        prog='get_data.py',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--version', action='version',
                        version='%(prog)s '+vers)
    # add required arguments
    required=parser.add_argument_group('required arguments')
    required.add_argument('outputfile',type=str,
                          help="output data file")
    # add optional arguments
    semi_opt=parser.add_argument_group('arguments set to defaults')
    semi_opt.add_argument('--snapnum',type=int,
                          help='snapnum to retrieve',
                          default=None)
    semi_opt.add_argument('--redshift_min',type=float,
                          help='minimum redshift',
                          default=None)
    semi_opt.add_argument('--redshift_max',type=float,
                          help='maximum redshift',
                          default=None)
    semi_opt.add_argument('--columns',type=str,
                          help='columns to obtain',
                          default='galaxyID,redshift,x,y,z,velX,velY,velZ,mag_r')
    semi_opt.add_argument('--limit',type=int,
                          help='limit number of returned galaxies',
                          default=None)

    # get command line arguments
    args = vars(parser.parse_args())
    # run with arguments
    main(args['outputfile'],redshift_min=args['redshift_min'],
         redshift_max=args['redshift_max'],snapnum=args['snapnum'],
         columns=args['columns'],limit=args['limit'])
