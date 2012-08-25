#!/usr/bin/env python2

# average_ensemble.py
# Author: Jonah Miller (jonah.maxwell.miller@gmail.com)
# Date: August 8, 2012

# This program loads a spacetime ensemble from a file, looks at
# 2-volume as a function of proper time, produces an average, and fits
# to a function

# Modules
#-------------------------------------------------------------------------
import numpy as np # For array control
import scipy as sp
import matplotlib as mpl # For data visualization
import matplotlib.pyplot as plt
import scipy.optimize as opt

import sys

# Written by me for spacetime data
import visualize_spacetime as vs
#-------------------------------------------------------------------------


# Constants
#-------------------------------------------------------------------------
def de_sitter_form_testing(x,s,A,C,N3,num_slices):
    """
    The functional form for de Sitter Space. x is position, s, A, C
    are all fit parameters. N3 is the total number of 3-simplices in
    the spacetime. num_slices is the number of time slices in the spacetime.
    """
    # A normalized cosine. used for clarity.
    n_cos = lambda y: np.cos(y/(s * N3**(1/3.)))

    # The full width/half-max of the spacetime
    fwhm = np.arccos(np.sqrt((np.pi * A/2.)*(s * N3**(1/3.)/N3)))
        
    # The conditional functions
    # for -num_slices/2.0 <= (x - C) < -s * N3**(1/3.) * fwhm 
    # or 
    # s * N3**(1/3.) * fwhm < (x-C) <= num_slices/2.0
    stem = lambda x: float(A)
    # for -s * N3**(1/3.) * fwhm <= (x - C) <= s * N3**(1/3.) * fwhm
    bulk = lambda x: (2/np.pi) * (N3/(s * N3**(1/3.))) * n_cos(x-C)**2

    # List of conditions for the piecewise function
    conds = [(-num_slices/2.0 <= (x-C))&((x-C) < -s * N3**(1/3.) * fwhm),
             (-s * N3**(1/3.) * fwhm <= (x-C))&((x-C) <= s * N3**(1/3.) * fwhm),
             (s * N3**(1/3.) * fwhm < (x-C))&((x-C) <= num_slices/2.0)]

    # List of return functions for the piecewise function
    returnfuncs = [stem,bulk,stem]

    return np.piecewise(x,conds,returnfuncs)
"""
#     
if -num_slices/2.0 <= (x - C) < -s * N3**(1/3.) * fwhm:
return A
elif -s * N3**(1/3.) * fwhm <= (x - C) <= s * N3**(1/3.) * fwhm:
return (2/np.pi) * (N3/(s * N3**(1/3))) * n_cos(x-C)**2
elif s * N3**(1/3.) * fwhm < (x-C) <= num_slices/2.0:
return A

"""
        
#-------------------------------------------------------------------------



# Main Method Functions
#-------------------------------------------------------------------------
def extract_2_volume_ensemble(filename_list):
    """
    Reads in each *.3sx2p1 file in filename_list and produces a list
    of lists of 2-volume information as a function of proper time.

    Returns a vector. The first element contains header data, which is
    assumed to be the same for all elements of the ensemble.
    """
    # Do the first spacetime seperately to extract header data
    first_spacetime = vs.read_3simplices_from_data_file(filename_list[0])
    header_data = first_spacetime[0]
    spacetimes = [first_spacetime[1]]

    # We need to calculate the average of N3-31 and N3-13 for all spacetimes
    N3 = eval(header_data[11])/2

    # Extract the remaining spacetimes
    for f in filename_list[1:]:
        temp = vs.read_3simplices_from_data_file(f)
        spacetimes.append(temp[1])
        N3 += eval(temp[0][11])/2

    # Extract the 2-simplex information 
    sl2simplices = [vs.get_all_sl2simplices([header_data,s]) \
                        for s in spacetimes]

    # Run garbage collection
    del spacetimes

    # Average the N3s instead of summing them
    N3 = N3/len(filename_list)

    # Extract volumes as a function of time
    volumes = [vs.make_v_of_t(s) for s in sl2simplices]

    # Delete the old lists
    del sl2simplices

    return [header_data,volumes,N3]


def statistical_average(volumes):
    """
    Takes an ensemble of 2-volume(proper time) and averages it,
    producing mean and standard deviation at each point. Returns a vector:
    [means, standard_deviations]
    """
    volumes = np.array(volumes).transpose()
    means = []
    stds = []
    for v_list in volumes:
        means.append(np.mean(v_list))
        stds.append(np.std(v_list))

    return [means,stds]


def fit_to_data(mean_volume_data,N3,guesses=[1,0,25]):
    """
    Takes mean volume data and attempts to fit to de Sitter
    spacetime. time_slices is the number of time slices, and N3 is the
    3-volume of the system.
    """

    # Total number of slices
    num_slices = len(mean_volume_data)


    def de_sitter_form(x,s,A,C):
        """
        The functional form for de Sitter Space. x is position, s, A, C
        are all fit parameters. N3 is the total number of 3-simplices in
        the spacetime. num_slices is the number of time slices in the spacetime.
        """
        # A normalized cosine. used for clarity.
        n_cos = lambda y: np.cos(y/(s * N3**(1/3.)))

        # The full width/half-max of the spacetime
        fwhm = np.arccos(np.sqrt((np.pi * A/2.)*(s * N3**(1/3.)/N3)))
        
        # The conditional functions
        # for -num_slices/2.0 <= (x - C) < -s * N3**(1/3.) * fwhm 
        # or 
        # s * N3**(1/3.) * fwhm < (x-C) <= num_slices/2.0
        stem = lambda x: float(A)
        # for -s * N3**(1/3.) * fwhm <= (x - C) <= s * N3**(1/3.) * fwhm
        bulk = lambda x: (2/np.pi) * (N3/(s * N3**(1/3.))) * n_cos(x-C)**2

        # List of conditions for the piecewise function
        """conds = [(-num_slices/2.0 <= (x-C)) & \
                     ((x-C) < -s * N3**(1/3.) * fwhm),
                 (-s * N3**(1/3.) * fwhm <= (x-C)) & \
                     ((x-C) <= s * N3**(1/3.) * fwhm),
                 (s * N3**(1/3.) * fwhm < (x-C)) & \
                     ((x-C) <= num_slices/2.0)]"""

        cond =   (-s * N3**(1/3.) * fwhm <= (x-C)) & \
                     ((x-C) <= s * N3**(1/3.) * fwhm)

        # List of return functions for the piecewise function
        #returnfuncs = [stem,bulk,stem]
        returnfuncs = [bulk,stem]

        #return np.piecewise(x,conds,returnfuncs)
        return np.piecewise(x,cond,returnfuncs)


    # Define the x-axis (time slices)
    xdata = np.array(range(num_slices),dtype=float)
    ydata = np.array(mean_volume_data,dtype=float)
#    err = np.array(std_deviation_data,dtype=float)

    # Fit!
    return opt.curve_fit(de_sitter_form,xdata,ydata,guesses)

def plot_and_fit(filename_list):
    "Read all the elements in the filename list and plot the results."

    # Header data contains relevant information about the insemble,
    # volumes is a list of 2-volumes, for all the spacetimes, and N3 =
    # (N3-31+N3-13)/2.
    header_data,volumes,N3 = extract_2_volume_ensemble(filename_list)

    # Get the mean and standard deviation of the volume list, delete
    # the list of volumes to save memory.
    means,stds = statistical_average(volumes)
    del volumes

    # Fit to data
    # popt = optimized parameters, pcov = covariance
    popt,pcov = fit_to_data(means,N3)

    # Make a curve 
    num_slices = float(len(means)) # The number of spacetime slices
    xdata = np.arange(0.0,num_slices,1.0) # x-axis
    fit = de_sitter_form_testing(xdata,popt[0],popt[1],popt[2],N3,num_slices)

    # Make a plot
    lines = [plt.errorbar(xdata,means,yerr=stds,fmt='ro',label='Monte Carlo'),
             plt.plot(xdata,fit,'b-',label='Fit')]

    # Make a legend
    plt.legend()

    # Set axes
    plt.xlabel('Proper Time')
    plt.ylabel('Spatial Extent')

    # Calculate and set title
    # suptitle is easy
    plt.suptitle('de Sitter Monte Carlo Fit, Fixed Boundaries')

    # subtitle uses header information
    raw_name = "#slices = {}. Volume = {}. k0 = {}. k3 = {}."
    name=raw_name.format(header_data[2],header_data[3],
                         header_data[-3],header_data[-2])
    plt.title(name)

    plt.show()

    # Calculate norm of the residuals: ||f(x,popt) - data||. The
    # smaller the number, the better the fit.
    rnorm = np.linalg.norm(fit - means) / float(np.max(means))

    return [popt,N3,num_slices,rnorm]

def print_popt(popt,N3,num_slices):
    "Prints the optimized parameters nicely."
    output = "N3 = {}. num_slices = {}.\ns = {}\nA = {}\nC = {}."
    print output.format(N3,num_slices,popt[0],popt[1],popt[2])

#-------------------------------------------------------------------------

if __name__ == "__main__":
    popt,N3,num_slices,rnorm = plot_and_fit(sys.argv[1:])
    print "Norm of the residuals: {}.".format(rnorm)
    print_popt(popt,N3,num_slices)

        
