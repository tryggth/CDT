#!/usr/bin/env python2

# get_subspacetime.py
# Author: Jonah Miller (jonah.maxwell.miller@gmail.com)
# Time-stamp: <2013-10-12 19:55:28 (jonah)>

# This program loads a spacetime file and generates a new file based
# on the old one. The new file is for a spacetime with initial
# boundary and final boundaries which are selected slices in the
# initial spacetime. A program call looks like this:
#
# pyton2 get_subspacetime.py 5 10 my_spacetime.3sx2p1
#
# You can also call it on as many spacetimes as you like if you want
# to generate the same subspacetime for each one:
#
# python2 get_subspacetime.py 10 15 my_firs_spacetime.3sx2p1 my_second_spacetime.3sx2p1
#
# And it supports globbing
#
# python2 get_subspacetime.py 10 15 *.3sx2p1
#
# The new spacetime file is the same as the old one, but with the
# initial and final slice numbers appended to the end. Like so:
#
# my_spacetime.subspacetime-5-10.3sx2p1
#
# IMPORTANT. The first line of the file is unchanged. It exists solely
# as a record of the original simulation. DONT'T RELY ON IT! IT IS NOT
# ACCURATE!

# Modules
# ----------------------------------------------------------------------
import sys # For globbing support
from copy import copy
# ----------------------------------------------------------------------


# Functions
# ----------------------------------------------------------------------
def new_file_name(old_file_name,initial_slice,final_slice):
    """
    This function takes old_file_name and makes it conform to the
    format for an output file name. initial slice is the initial slice
    of the subspacetime. final_slice is the final slice of the
    subspacetime.
    """
    subspacetime_designation='subspacetime-{}-{}.'.format(initial_slice,
                                                          final_slice)
    output=old_file_name[:-6]+subspacetime_designation+old_file_name[-6:]
    return output

def make_output_string(simplex,initial_slice):
    """
    takes a simplex list of the form

    [type, tmlo, tmhi, point1, point2, point3, point4, neighbor1,
     neighbor2, neighbor3, neighbor4, simplexID]

    and outputs it as a space-separated list that's a string.
    """
    # The simplex may have the incorrect the incorrect tmhi and tmlo
    # so wee need to fix that.
    output_list = copy(simplex)
    output_list[1] = int(output_list[1])-initial_slice
    output_list[2] = int(output_list[2])-initial_slice

    return reduce(lambda x,y: "{} {}".format(x,y),output_list) + "\n"

def in_bulk(simplex,initial_slice,final_slice):
    """
    takes a simplex list of the form

    [type, tmlo, tmhi, point1, point2, point3, point4, neighbor1,
     neighbor2, neighbor3, neighbor4, simplexID]

    And checks whether or not initial_slice < tmlo and final_slice > tmhi.

    This excludes simplices on the boundary of the new subspacetime.    
    """
    tmlo = simplex[1]
    tmhi = simplex[2]
    return tmlo > initial_slice and tmhi < final_slice

def on_boundary(simplex, initial_slice,final_slice):
    """
    takes a simplex list of the form

    [type, tmlo, tmhi, point1, point2, point3, point4, neighbor1,
     neighbor2, neighbor3, neighbor4, simplexID]

    And checks whether or not the simplex is on the boundary of the
    subspacetime.
    """
    tmlo = simplex[1]
    tmhi = simplex[2]
    return tmlo == initial_slice or tmhi == final_slice

def make_new_first_line(first_line,initial_slice, final_slice):
    """
    Takes the string that's the first line of the file that has
    parameter information and resets the number of time slices.
    """
    line_list = first_line.split(' ')
    line_list[2] = str(final_slice - initial_slice)
    return reduce(lambda x,y: "{} {}".format(x,y),line_list)    

def extract_subspacetime(old_file_name,initial_slice,final_slice):
    """
    Takes old_spacetime and extracts the subspacetime bounded by
    initial_slice and final_slice. Outputs a file with the subspacetime.
    """

    # Open the old and new spacetimes in an error safe way.
    outfile_path=new_file_name(old_file_name,initial_slice,final_slice)
    with open(old_file_name,'r') as infile:
        with open(outfile_path,'w') as outfile:

            # Copy the first line if the infile to the outfile. This
            # is the list of simulation parameters.
            first_line = make_new_first_line(infile.readline(),
                                            initial_slice,
                                            final_slice)
            outfile.write(first_line)

            # Now iterate through the file and read each remaining
            # line. Each line represents a 3-simplex. It's a space
            # separated list containing the following simplex
            # information in this order:
            #
            # type tmlo tmhi point1 point2 point3 point4 neighbor1 neighbor2
            # neighbor3 neighbor4 simplexID
            #
            # We copy the line from the old file to the new file if
            # and only if tmlo is greater than initial_slice and tmhi
            # is less than final_slice.
            #
            # Unfortunately, we need to copy the spacetime
            # carefully. If a simplex has a "neighbor" that's outside
            # the boundary time slices, we need to remove this
            # neighbor. Neighbors "zero" are considered void pointers
            # and are ignored. So we need to ensure that neighbors are
            # treated this way.

            # First we read in the file and store each line
            lines = infile.readlines()
            # Then we make a hash table that maps a simplex id to the
            # apppropriate simplex and the string in associated with
            # it in the file
            simplices = {}
            for line in lines:
                simplex = [int(i) for i in line.split()]
                simplices[simplex[-1]] = simplex

            # Now, we go back through each line and, if the simplex is
            # in the bulk, print it to the new outfile. If the simplex
            # is on the boundary, ensure all its neighbors are on the
            # boundary too. Then print it to the file.
            for simplex_id in simplices.keys():
                simplex = simplices[simplex_id]
                if in_bulk(simplex,initial_slice,final_slice):
                    outfile.write(make_output_string(simplex,initial_slice))
                if on_boundary(simplex,initial_slice,final_slice):
                    for i in range(7,11):
                        if simplex[i] != 0:
                            neighbor = simplices[simplex[i]]
                            if not (in_bulk(neighbor,initial_slice,
                                            final_slice)\
                                        or on_boundary(neighbor,
                                                       initial_slice,
                                                       final_slice)):
                                simplex[i] = 0
                    outfile.write(make_output_string(simplex,initial_slice))
    return outfile_path

def main(system_arguments):
    """
    Takes the system arguments (i.e., sys.argv), which include the
    command line arguments, and extracts the spacetimes for the input
    files.
    """

    print "Welcome to get_subspacetime.py!"

    # script_name=system_arguments[0]
    initial_slice=int(system_arguments[1])
    final_slice=int(system_arguments[2])
    input_spacetime_filenames=system_arguments[3:]

    print "You want spacetimes sandwiched between slices {} and {}.".format(initial_slice,final_slice)
    print "Extracting subspacetimes..."
    print "The new file names are:"

    for spacetime_file in input_spacetime_filenames:
        new_file_name=extract_subspacetime(spacetime_file,
                                           initial_slice,
                                           final_slice)
        print new_file_name

    print "All finished. Have fun!"

# ----------------------------------------------------------------------

if __name__ == "__main__":
    main(sys.argv)
