#!/usr/bin/env python2

"""
make_thermalization_scripts_initial_final_boundary_pairs.py
Author: Jonah Miller (jonah.maxwell.miller@gmail.com)

This program takes the fullpaths fo two directories as input, each
directory containing an equal number of *.boundary2p files. If one
directory has more files than the other, only makes pairs until it
runs out of files in the smaller directory. The program randomly
matches the files into ordered pairs. It then generates thermalization
scripts that use the *boundary2p1 files as initial and final
boundaries. The first argument is the initial boundary and the second
argument is the final boundary. File names are a concatenation of the
two boundary files.

Example call:
-------------------------------------
python2 make_thermalization_scripts_initial_final_boundary_pairs.py $(pwd)/initial_boundaries $(pwd)/final boundaries

"""

# Import modules
#----------------------------------------------------------------------
from random import shuffle # for randomizing the order ofcontents of a
                           # directory
import sys,os # For the operating system/file system commands
import make_thermalization_scripts as script_maker # Uses my previous codebase
#----------------------------------------------------------------------

# Functions we need
#----------------------------------------------------------------------
def make_boundary_filename(boundary_file_1,boundary_file_2):
    "Utility to help interface with "
    boundary_file_name = "IBN_"+\
        script_maker.make_boundary_file_out(boundary_file_1) \
         + "_FBN_" + boundary_file_2
    return boundary_file_name

def main(command_line_arguments):
    "Takes command linke arguments and makes the scripts you need."

    # Ensure the correct number of arguments given.
    if len(command_line_arguments) < 2:
        print "You have the wrong number of command line arguments!"
        print command_line_arguments
        raise ValueError("Invalid command line arguments.")
    
    # Local variables
    initial_boundary_directory = command_line_arguments[0]
    final_boundary_directory = command_line_arguments[1]
    initial_boundaries = os.listdir(initial_boundary_directory)
    final_boundaries = os.listdir(final_boundary_directory)

    print "Making scripts!"
    print "Initial boundary directory: {}".format(initial_boundary_directory)
    print "Final boundary directory: {}".format(final_boundary_directory)

    # Randomize the order of the directories
    shuffle(initial_boundaries)
    shuffle(final_boundaries)

    # Make new scripts with ordered-pair filenames
    for i in range(min(len(initial_boundaries),len(final_boundaries))):
        boundary_filename = make_boundary_filename(initial_boundaries[i],
                                                   final_boundaries[i])
        script_maker.make_output(boundary_filename,initial_boundaries[i],
                                 final_boundaries[i],i)
    print "All done! Happy hacking!"

#----------------------------------------------------------------------

# Main loop
#----------------------------------------------------------------------
if __name__ == "__main__":
    main(sys.argv[1:])
                                 
    
