"""
output.py

Time-stamp: Time-stamp: <2012-09-27 21:54:35 (jonah)>

Author: Jonah Miller (jonah.maxwell.miller@gmail.com)

This file is part of the sphere_generator program which generates
spheres of (as close as possible to) uniform curvature for a given
surface area by monte carlo simulation.

This file contains the functions for printing a sphere that's been
generated to a file.
"""

# Dependancies
#------------------------------------------------------------------------
import numpy as np
import scipy as sp
import random
import datetime
# Class data structures we need
import simplex_ancestors as sa
import simplex_descendants as sd
import utilities as ut
import state_tracking as st
#------------------------------------------------------------------------


### CONSTANTS
###----------------------------------------------------------------------
output_suffix = ".boundary2p1" # The file ending for files generated.
tracking_suffix = ".boundaryprg2p1" # The file ending for progress files
output_prefix = "S2_" # The file prefix for files generated.
###----------------------------------------------------------------------



### TEMPLATE FUNCTIONS
###----------------------------------------------------------------------
def make_sphere_data_list():
    """
    Generates a list of sets, each set containing the vertices of a
    triangle. Useful for output. Can be reindexed by
    reindex_sphere_data_list.
    """
    return [t.get_vertex_ids() for t in sd.triangle.instances.values()]

def reindex_sphere_data_list(sphere_data_list):
    """
    Takes an input list of the form [set([v1,v2,v3]),...] and ensures
    that the list of vertices has no hooes. i.e., if the maximum
    vertex id is n, then the set of vertices is
    set(range(1,n+1)). Returns a reindex list. Used for output so that
    the boundary file contains a list we can give to
    2p1-fixed-boundaries as a boundary.
    """
    # Generate a sorted list of all the vertex ids. The reindexed id
    # of each vertex is just the id of the vertex in the list.
    vertex_ids = list(ut.set_union(sphere_data_list))
    vertex_ids.sort()
    # It's troublesome to find the index of the id in the list this
    # way. A better data structure is a dictionary:
    reindexing_scheme = {v:vertex_ids.index(v)+1 for v in vertex_ids}
    # We can now reconstruct our reindexed list by mapping v to
    # reindexing_scheme[v]
    sphere_data_list = [set([reindexing_scheme[v] for v in s]) \
                            for s in sphere_data_list]
    return sphere_data_list
    
def make_triangle_point_string(triangle_points):
    """
    Takes a collection of 3 vertices, (v1,v2,v3), for instance, and
    makes it into a string that LISP can handle.
    """
    assert len(triangle_points) == 3
    triangle_points = list(triangle_points)
    return "({} {} {})".format(triangle_points[0],triangle_points[1],
                               triangle_points[2])

def make_output_string():
    """
    Generates a string containing information about the current sphere
    in a format that 2p1-fixed-boundary simulations can accept.
    """
    sphere_data_list = reindex_sphere_data_list(make_sphere_data_list())
    outstring = "("
    for t in sphere_data_list:
        outstring += make_triangle_point_string(t)+" "
    outstring = outstring[:-1] + ")"
    return outstring

def write_sphere_to_file(filename,overwrite=True):
    """
    Writes the current sphere to a file named filename. No progress is
    preserved. The output file is exactly what's required as input for
    2p1-fixed-boundary simulations. Overwrites the current file. If
    overwrite is false, it appends, instead.
    """
    if overwrite:
        writestyle = 'w'
    else:
        writestyle = 'a'
    with open(filename,writestyle) as f:
        f.write(make_output_string())
        f.write('\n')

def write_progress_file(metropolis_state,final_sweep,save_every_n_sweeps,
                        filename):
    """
    Makes the progress file for the program. metropolis_state is an
    instance of a descendant class of the metropolis class.
    """
    with open(filename,'w') as f:
        f.write(metropolis_state.make_progress_file_output(final_sweep,save_every_n_sweeps))

def make_time_string():
    """
    Makes a string out of the current time
    """
    now = datetime.datetime.today()
    nowstr = str(now).split(' ')[0] + '_' + str(now).split(' ')[1]
    return nowstr

###----------------------------------------------------------------------


### DATA GATHERING FUNCTIONS
###----------------------------------------------------------------------
def gather_data_to_1_file(metropolis_state,final_sweep,
                          initial_sweep=0,
                          save_every_n_sweeps=10):
    """
    Runs a monte carlo simulation using the metropolis state (a class
    instance of a descendant class of metropolis). 

    Saves one sphere file and one progress file.
    """
    # Reset the current sweep so that we start from the specified sweep
    metropolis_state.reset_current_sweep(initial_sweep)

    # Make the filename
    filename = metropolis_state.make_file_name_v1(final_sweep)
    # Immediately make a progress file and a current state
    write_sphere_to_file(filename+output_suffix)
    write_progress_file(metropolis_state,final_sweep,save_every_n_sweeps,
                        filename+tracking_suffix)

    # Rewrite each file every save_every_n_sweeps
    while metropolis_state.get_current_sweep() < final_sweep:
        metropolis_state.sweep(save_every_n_sweeps)
        write_sphere_to_file(filename+output_suffix)
        write_progress_file(metropolis_state,final_sweep,save_every_n_sweeps,
                            filename+tracking_suffix)

    # Print a happy message
    print "All done! :)"

def gather_data_to_n_files(metropolis_state, final_sweep,
                          initial_sweep=0,
                          save_every_n_sweeps=10):
    """
    Runs a monte carlo simulation using the metropolis state (a class
    instance of a descendant class of metropolis). 

    Saves as many sphere files as necessary. There's no progress file.
    """
    # Reset the current sweep so that we start from the specified sweep
    metropolis_state.reset_current_sweep(initial_sweep)

    # Start time!
    start_time = make_time_string()

    # Immediately save a file
    fname = metropolis_state.make_file_name_v2(final_sweep,start_time)\
        + output_suffix
    write_sphere_to_file(fname)
    while metropolis_state.get_current_sweep() < final_sweep:
        metropolis_state.sweep(save_every_n_sweeps)
        fname = metropolis_state.make_file_name_v2(final_sweep,start_time)\
            + output_suffix
        write_sphere_to_file(fname)

    # Print a happy message
    print "All done! :)"

###----------------------------------------------------------------------
