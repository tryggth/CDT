#!/usr/bin/env python2

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
# Class data structures we need
from simplex_ancestors import *   
from simplex_descendants import *
import state_manipulation as sm 
import utilities as ut
import error_checking
#-------------------------------------------------------------------------


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
        point.connect_surrounding_triangles()

def build_first_tetrahedron():
    """
    The base of initialization. Builds a tetrahedron to operate moves on.

    The points of a triangle, grouped by triangle, are:
    ((4 3 2) (4 1 3) (1 4 2) (2 1 3))
    """
    # Ensure no sphere is currently initialized
    sm.delete_all_geometries()

    # The ids of the points we will use
    point_names = [1,2,3,4]

    # The number of points in a tetrahedron
    n_points = len(point_names)
    
    # The triangles we'll build, defined by their points. Each element
    # is a list of points for a triangle.
    triangles_by_points = [[4,3,2], 
                           [4,1,3],
                           [1,4,2],
                           [2,1,3]]

    # Make the points we need
    sm.make_n_vertices(n_points)
    
    # Double check the points we have are what we think they are
    assert list(vertex.instances.keys()) == point_names

    # Make the triangles we need out of those points.
    triangle_ids = [sm.build_triangle_and_edges(t) \
                        for t in triangles_by_points]

    # Connect all the triangles and tell points what triangles contain them
    connect_all_triangles()

    return triangle_ids
### -------------------------------------------------------------------------

