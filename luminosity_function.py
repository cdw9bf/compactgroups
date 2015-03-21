"""
Create luminosity function
Trey Wenger
"""
import numpy as np
import matplotlib as mpl
mpl.rc_file('/home/twenger/matplotlibrc')
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def schechter(mag, numstar, magstar, alpha):
    """
    Return number of galaxies at a magnitude
    """
    num = 0.5*np.log(10.)*numstar
    num *= (10.**(0.4*(magstar - mag)))**(1.+alpha)
    num *= np.exp(-10.**(0.4*(magstar - mag)))
    return num

data = np.genfromtxt('datafile_63.txt',dtype=None,names=True,
                     delimiter=',')
good_data = data['mag_r'] < 99.
# histogram
n,bins,patches = plt.hist(data['mag_r'][good_data],40)
mags = np.array([(bins[i] + bins[i+1])/2. for i in range(len(bins)-1)])
# fit schechter function
