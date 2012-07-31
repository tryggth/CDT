#!/usr/bin/env python2

# invert_time.py
# Author: Jonah Miller (jonah.maxwell.miller@gmail.com)
# Date: July 30, 2012

# This program loads a spacetime file and inverts the proper time
# coordinates in each simplex. Useful for certain tests.

# Modules
#-------------------------------------------------------------------------
import sys # For yes/no
import numpy as np # For array control
import scipy as sp
from visualize_spacetime import read_3simplices_from_data_file
#-------------------------------------------------------------------------


# Functions
#-------------------------------------------------------------------------
def extract_t_max(header_data):
    """
    Takes the header data from a spacetime read in by
    read_3simplices_from_data-file
    and extracts the maximum time coordinate.
    """
    openbc = header_data[0] == "OPEN"
    if openbc:
        t_max = int(header_data[2])
    else:
        t_max = int(header_data[1] - 1)
    return t_max

def invert_time_coordinate(t_now,t_max):
    """
    Takes t_now and inverts it so that we count up from t_max.
    """
    t_new = t_max - t_now
    return t_new

def invert_time_for_1_simplex(simplex, t_max):
    """
    Takes in a list defining a simplex, and inverts the time coordinate
    based on the known t_max of the spacetime. This means reversing the 
    the order of the points p1...p4 and the order of neighbors n1...n4.
    It also means renaming the type of simplex.
    """
    
    # The new list of points is the old list of points reversed.
    new_points = simplex[3:7]
    new_points.reverse()
    
    # The new list of neighbors is the old list of neighbors reversed.
    new_neighbors = simplex[7:11]
    new_neighbors.reverse()
    
    # tmlow and tmhigh are inverted and their positions with respect
    # to each other are reversed.
    new_tmlo = invert_time_coordinate(simplex[2],t_max)
    new_tmhi = invert_time_coordinate(simplex[1],t_max)
    
    # The id is identical
    new_sxid = simplex[-1]
    
    # If the simplex is a (2,2), it's type doesn't change. If the
    # simplex is a (1,3) or a (3,1) type, it becomes the other type.
    if simplex[0] == 1:
        new_type = 3
    elif simplex[0] == 2:
        new_type = 2
    elif simplex[0] == 3:
        new_type = 1
    else:
        print "Error! Unkown simplex" \
            + "type for simplex ID {}".format(new_sxid)
    
    # Generate output
    new_simplex = [new_type,new_tmlo,new_tmhi] + new_points + new_neighbors \
        + [new_sxid]
    return new_simplex

def invert_time_for_spacetime(spacetime_data):
    """
    Invert the time coordinate for every 3-simplex in the
    spacetime. The input is the entire spacetime data from
    read_3simplices_from_data_file, not just the simplex information.
    """
    t_max = extract_t_max(spacetime_data[0]) # maximum time index in spacetime
    new_spacetime = [] # An empty list to put our spacetime in
    
    # Fill our spacetime up
    for simplex in spacetime_data[1]:
        new_spacetime.append(invert_time_for_1_simplex(simplex,t_max))
    
    # Return the fancy pants new spacetime with inverted time coordinates:
    return new_spacetime

def print_list_to_file(L,f):
    """
    Prints list L to file f on a single line, space separated. 
    f is an object, not a file name.
    """
    linestring = ""
    for element in L:
        linestring += str(element) + ' '
    linestring = linestring[:-1] + '\n'
    f.write(linestring)

def print_output_file(header_data,simplex_data,filename):
    """
    Makes an output file based on 
    filename.3sx2p1 -> filename.time-inverted.3sx2p1
    and prints simplex data to it.
    """
    new_fname = filename[:-7] + ".time-inverted" + filename[-7:]
    with open(new_fname,'w') as f:
        print_list_to_file(header_data,f)
        for simplex in simplex_data:
            print_list_to_file(simplex,f)

def main(filename):
    """
    Takes a spacetime file, inverts the coordinates, and outputs an
    inverted file.
    """
    spacetime_data = read_3simplices_from_data_file(filename)
    new_spacetime = invert_time_for_spacetime(spacetime_data)
    print_output_file(spacetime_data[0],new_spacetime,filename)
# -------------------------------------------------------------------------


# Main loop
if __name__ == "__main__":
    print "Reversing time coordinates..."
    for i in sys.argv[1:]:
        main(i)
    print "All done! Good luck with your simulations!"
