"""
initialization.py

Author: Jonah Miller (jonah.maxwell.miller@gmail.com)

This file is part of the sphere_generator program which generates
spheres of (as close as possible to) uniform curvature for a given
surface area by monte carlo simulation.

This module initializes the sphere at the start of runtime.
"""


### Dependencies
#-------------------------------------------------------------------------
import numpy as np
import scipy as sp
import random
# Class data structures we need
from simplex_ancestors import *   
from simplex_descendants import *
import state_manipulation as sm
import utilities as ut
import error_checking
import moves

#-------------------------------------------------------------------------


### CONSTANTS
###----------------------------------------------------------------------
tetrahedron_point_data = [set([4,3,2]),
                          set([4,1,3]),
                          set([1,4,2]),
                          set([2,1,3])]
###----------------------------------------------------------------------


### Functions required during initialization
### -------------------------------------------------------------------------
def connect_all_triangles():
    """
    Likely only useful at initialization. This function ensures taht
    every vertex knows what triangles it is part of and every triangle
    knows what neighbors it has.
    """
    for point in vertex.instances.values():
        point.find_and_set_triangles()
        point.find_and_set_edges()
        point.connect_surrounding_triangles()

def build_sphere_from_data(sphere_data):
    """
    sphere_data should be a list of sets, each with 3 points in
    it. Resets the sphere and then builds it from this data.
    """
    # First we need to check that sphere_data is what we want.
    for s in sphere_data:
        if not (type(s) == set and len(s) == 3):
            raise TypeError("Sphere data needs to be a list of sets, "
                            +"each of length 3!")
        for i in s:
            if not (type(i) == int and i > 0):
                raise TypeError("The vertices must all be "
                                +"positive integers!")

    # The vertices in existence will be the union of the triangle
    # point lists:
    vertex_ids = ut.set_union(sphere_data)
    # The last used vertex id is the maximum of the vertex ids
    last_used_vertex_id = max(vertex_ids)

    # make sure there are no vertices are missing (print sphere to
    # file shouldn't let there be any)
    assert set(range(1,last_used_vertex_id+1)) - vertex_ids == set([])

    # Once we're sure the set is what we want, we need to make sure no
    # sphere is currently initialized.
    sm.delete_all_geometries()
    # Make the points we need
    sm.make_n_vertices(last_used_vertex_id)
    # Ensure the points we made are the points in our list
    assert set(vertex.instances.keys()) == vertex_ids

    # Make the triangles we need out of the points.
    triangle_ids = [sm.build_triangle_and_edges(t) for t in sphere_data]

    # Connect all the triangles and tell points what triangles contain them
    connect_all_triangles()

    return triangle_ids
    
            

def build_first_tetrahedron():
    """
    The base of initialization. Builds a tetrahedron to operate moves on.

    The points of a triangle, grouped by triangle, are:
    ((4 3 2) (4 1 3) (1 4 2) (2 1 3))
    """
    return build_sphere_from_data(tetrahedron_point_data)

def increase_volume_once():
    """
    Uses random volume increasing moves to increase the surface area once.
    """
    # Try a random volume increasing move on a random simplex
    mdata = moves.try_random(moves.list_of_volume_increasing_functions)

    # If the movedata exists, apply the move
    if mdata:
        moves.apply_move(mdata)

def increase_volume_to_target(target_volume):
    """
    Uses random volume increasing moves to increase the surface area
    to the target surface area.
    """
    while triangle.count_instances() < target_volume:
        increase_volume_once()

def initialize_sphere(target_volume):
    """
    Initialize a sphere with surface area equal to target_volume.
    """
    # Ensure that there's no sphere currently in place.
    sm.delete_all_geometries()
    # Build the first tetrahedron
    build_first_tetrahedron()
    # Increase set the surface area randomly
    increase_volume_to_target(target_volume)

### -------------------------------------------------------------------------

    
# Functions for loading a sphere from a file
###-------------------------------------------------------------------------
def parse_sphere_file(file_contents):
    """
    Takes a string, file_contents, as input and parses it into point
    lists that can be made into triangle objects.

    INPUT FORMAT:
    "((v1 v2 v3) (v1 v4 v3) ... )"

    RETURN FORMAT:
    [set([v1, v2, v3]), set([v1, v4, v3]), ... ]
    """
    # We want to make sure we're operating on what we think we're
    # operating on.
    assert type(file_contents) == str
    # The first two characters of the string will be "((" we don't
    # want them.  The last three characters of the string will be
    # "))\n". We don't want them either.
    file_contents = file_contents.lstrip("(").rstrip('\n').rstrip(")")
    # The string is now space-separated numbers partitioned by ") (".
    file_contents = file_contents.split(") (")
    # We now have a list, where each element is 3 space-separated numbers.
    file_contents = [s.split(' ') for s in file_contents]
    # We now have a list, where each element is a list of string
    # numbers. We need to convert them into real numbers (integers
    # only. If they're not integers, raise an error.)
    file_contents = [set([int(s) for s in l]) for l in file_contents]
    # We're done! Return it!
    return file_contents

def get_data_from_file(filename):
    """
    Gets the data for a sphere (a list of triangle points) from the
    file named filename. filename is a string.
    """
    with open(filename,'r') as f:
        file_contents = f.read()
    return parse_sphere_file(file_contents)

def load_sphere_from_file(filename):
    """
    Takes sphere data from a file and builds a brand new sphere using
    the sphere data.
    """
    return build_sphere_from_data(get_data_from_file(filename))

###-------------------------------------------------------------------------
