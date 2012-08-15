#!/usr/bin/env python2

"""
Author: Jonah Miller (jonah.miller@colorado.edu)

This Program takes a *.3sx2p1 file and a pair of time slices as input
and generates text files that can be used by the
cdt2p1-fixed-boundaries code to generate a spacetime using the chosen
time slices of a previous simulation as the boundaries.

Example program call:
./make_boundary_conditions.py t_low t_high filename.3sx2p1

t_low is the time slice number in the spacetime for t=0
boundary. Doesn't have to be the low time slice number in the input
data.

t_high is analogous.
"""

# Modules
#-------------------------------------------------------------------------
import sys # For globbing support
import time # For a cute welcome message.
import numpy as np # For array control
import scipy as sp 
import matplotlib as mpl # For data visualization
import matplotlib.pyplot as plt
# Visualize spacetime was written by me to extract 2-simplex data from
# a *.3sx2p1 file. It's useful for this reason. Because we re-use
# code.
import visualize_spacetime as vs 
#-------------------------------------------------------------------------


### Functions
###-------------------------------------------------------------------------
## Functions for input
def get_raw_time_slices(filename,t_low,t_high):
    """
    Given the filename/filepath of a *.3sx2p1 file, extract the
    2-simplices for the high and low simplices we'll use to make our
    boundary conditions as a pair of lists (in a list).
    """

    # Load spacetime
    spacetime = vs.read_3simplices_from_data_file(filename)
    # Assign time slices
    t_slice_low = vs.get_spacelike_2simplices_at_t(spacetime[1],t_low)
    t_slice_high = vs.get_spacelike_2simplices_at_t(spacetime[1],t_high)

    # spacetime[0] contains coupling constant and spacetime
    # size. Useful for file naming conventions.
    return [t_slice_low, t_slice_high,spacetime[0]]


## Functions for output
def generate_filename(simulation_data,tslice, side):
    """
    Makes a file name. Format:
    T<time slices>-V<3-volume>-k<k0>-kkk<k3>-<side>.boundary
    """
    ts = simulation_data[1] # Time slices
    v3 = simulation_data[2] # 3-volume

    # Couping constant position depends on whether data came from HL
    # simulation or 2+1 Einsteinian simulation.
    if len(simulation_data) == 19:
        k0 = simulation_data[-5]
        k3 = simulation_data[-4]
    else:
        k0 = simulation_data[-3]
        k3 = simulation_data[-2]

    filename = 'T'+ts+'-V'+v3+'-k'+k0+'-kkk'+k3+'-tslice'+tslice+'-'+side+'.boundary'

    return filename

def make_lispy_triangle(triangle):
    """
    Takes in a list that represents a triangle. Returns a string that
    lisp interprets as a triangle (a lisp list).
    """
    return '(%i %i %i)' % (triangle[0], triangle[1], triangle[2])

def make_output_string(time_slice):
    """
    Given a time slice and generates a string that the lisp
    interpreter interprets as a list that th 2p1-fixed-boundary code
    can interpret as initial conditions.
    """
    outstring = "(" # The string we'll use as output.  

    # For every triangle in the time slice, make a lispy triangle and
    # add it ot the output string.
    for triangle in time_slice:
        outstring = outstring + make_lispy_triangle(triangle) + ' '

    # Close the output string as a lisp list, removing any trailing
    # whitespace.
    outstring = outstring.rstrip(' ') + ')'

    return outstring

def write_output(out_file,time_slice):
    """
    Given an output file name and a reindexed time_slice, produce a
    file that the lisp interpreter can read as boundary condition
    data.
    """
    
    outstring = make_output_string(time_slice)

    with open(out_file,'w') as f:
        # The endline character (\n) is for linux utilities. Lisp
        # doesn't care.
        f.write(outstring + '\n') 
    

## Functions for actual calculations
def remove_duplicates(L):
    """
    Removes duplicates from the list L. Ordering is lost.
    """
    return list(set(L))

def make_point_list(time_slice):
    """
    Takes a time slice with triangles in it and extracts a sorted list
    of points, with no duplicates.
    """

    # First, take all the points in the time slice, and put them in a
    # list. Since the points are in lists of length 3 (for each
    # triangle), we need two levels of extraction.
    point_list = [] # We'll add points to this
    for triangle in time_slice:
        for point in triangle:
            point_list.append(point)

    # There are duplicate points in our point list, since this is how
    # we connect triangles (several triangles have the same point as a
    # vertex). To do our operations here, we don't want
    # duplicates. Thus, we cull them.
    point_list = remove_duplicates(point_list)

    # Sort the point list, since we want to re-index from lowest to highest.
    point_list.sort()

    # Finally send it out
    return point_list

def reindex_points(time_slice,starting_point = 1):
    """
    Takes a set of triangles (2-simplices) in a time slice and
    re-indexes the points from to N+M, where M is the total number of
    points, and N, is a pre-defined starting point index
    (starting_point). Returns a reindexed time slice in the same
    format as input.
    """
    
    # Get a sorted list of points with no duplicates
    point_list = make_point_list(time_slice)

    # Now, we want to set up a 1:1 correspondance between the points
    # in our time-slice with a new indexing scheme from N to
    # N+M. This is the job for a dictionary.
    reindexing_scheme = {} # initialize the dictionary
    for i in range(len(point_list)): 
        reindexing_scheme[point_list[i]] = i + starting_point

    # Now generating our re-indexed set of simplices is just a matter
    # of going through the original time_slice list and replacing
    # using the dictionary.
    for triangle in time_slice:
        for i in range(len(triangle)):
            triangle[i] = reindexing_scheme[triangle[i]]

    return time_slice

## Main function
def main(filename, t_low, t_high):
    """
    Reads a file, generates the simplices
    """

    # Welcome message:
    print "Welcome to make_boundary_conditions.py!\nThe time is:"
    print time.ctime()
    print "Making your boundary conditions for:\n"+filename+"\n..."

    # Calculations
    raw_data = get_raw_time_slices(filename,t_low,t_high)
    sim_data = raw_data[2]
    raw_t_slices = raw_data[:2]
    reindexed_time_slices = [reindex_points(i) for i in raw_t_slices]

    # Output
    write_output(generate_filename(sim_data,str(t_low),'bottom'),
                 reindexed_time_slices[0]) # Bottom slice
    write_output(generate_filename(sim_data,str(t_high),'top'),
                 reindexed_time_slices[1]) # Top slice

    print "Done! Good luck with your data. :)"

    return              
    
#-------------------------------------------------------------------------

if __name__ == "__main__":
    ts_low = int(sys.argv[1])
    ts_high = int(sys.argv[2])
    for filename in sys.argv[3:]:
        main(filename,ts_low,ts_high)
