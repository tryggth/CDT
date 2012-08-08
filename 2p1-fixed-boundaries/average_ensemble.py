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

# Written by me for spacetime data
import visualize_spacetime as vs
#-------------------------------------------------------------------------


# Constants
#-------------------------------------------------------------------------
def de_sitter_form(x,s,A,C,N3,num_slices):
    """
    The functional form for de Sitter Space. x is position, s, A, C
    are all fit parameters. N3 is the total number of 3-simplices in
    the spacetime. num_slices is the number of time slices in the spacetime.
    """
    # A normalized cosine. used for clarity.
    n_cos = lambda y: np.cos(y/(s * N3**(1/3)))

    # The full width/half-max of the spacetime
    fwhm = np.arccos(np.sqrt((np.pi * A/2)*(s * N3**(1/3)/N3)))

    # The piecewise function we will return.
    if -num_slices/2.0 <= (x - C) < -s * N3**(1/3.) * fwhm:
        return A
    elif -s * N3**(1/3.) * fwhm <= (x - C) <= s * N3**(1/3.) * fwhm:
        return (2/np.pi) * (N3/(s * N3**(1/3))) * n_cos(x-C)**2
    elif s * N3**(1/3.) * fwhm < (x-C) <= num_slices/2.0:
        return A

    # In case something goes very wrong
    else:
        return A
        
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
    # Extract the spacetimes
    spacetimes = [vs.read_3simplices_from_data_file(f) for f in filename_list]

    # Extract the 2-simplex information 
    sl2simplices = [vs.get_all_sl2simplices(s) for s in spacetimes]

    # Extract volumes as a function of time
    volumes = [make_v_of_t(s) for s in sl2simplices]

    # Extract header data so we can delete old lists
    header_data = spacetimes[0][0]

    # Delete the old lists
    del spacetimes
    del sl2simplices

    return [header_data,volumes]


def statistical_average(volumes):
    """
    Takes an ensemble of 2-volume(proper time) and averages it,
    producing mean and standard deviation at each point. Returns a vector:
    [means, standard_deviations]
    """
    means = []
    stds = []
    for v_list in volumes:
        means.append(np.mean(v_list))
        stds.append(np.std(v_list))

    return [means,stds]

def fit_to_data(mean_volume_data,N3,num_slices):
    """
    Takes mean volume data and attempts to fit to de Sitter
    spacetime. time_slices is the number of time slices, and N3 is the
    3-volume of the system.
    """
    # Cut down on de_sitter_form so that N3 and num_slices are taken
    # into account
    func = lambda x,s,A,C: de_sitter_form(x,s,A,C,N3,num_slices)

    # Define the x-axis (time slices)
    x = range(num_slices)

    # Fit!
    return opt.curve_fit(de_sitter_form,x,mean_volume_data)


