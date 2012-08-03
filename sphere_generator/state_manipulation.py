#!/usr/bin/env python2

"""
simplex_descendants.py

Author: Jonah Miller (jonah.maxwell.miller@gmail.com)

This file is part of the sphere_generator program which generates
spheres of (as close as possible to) uniform curvature for a given
surface area by monte carlo simulation.

This module contains classes, functions, and methods which give or
manipulate the state of the simulation. The primary class is the
sphere class, which is meant to be initialized and hold the general
state of the simulation. Name spaces are cool.
"""



### Dependencies
#-------------------------------------------------------------------------
import numpy as np
import scipy as sp
# Class data strucbures we need
from simplex_ancestors import *   
from simplex_descendants import *
#-------------------------------------------------------------------------


####---------------------------------------------------------------------####
#                               Classes                                     #
####---------------------------------------------------------------------####
# The sphere base class, which gives the state of the simulation
#---------------------------------------------------------------------------
class sphere:
    """
    Contains and calculates the state information for the sphere at a
    given time-step. Meant to be initialized.

    It's worth noting that printing the sphere instances is useful.

    Methods:
    euler_characteristic() = Calculates the euler characteristic of 
                             the sphere.
    curvature_total(normalized=False) = Calculates the total Gauss curvature 
                                        of the sphere. If normalized=True, 
                                        also divide by the number of vertices
                                        to get an average.
    curvature_std() = Calculates standard deviation of the Gauss curvature 
                      over the sphere. 
    """
    
    def __init__(self):
        # perhaps add sphere initialization routine call here
        pass 
    
    def euler_characteristic(self):
        "Calculates the euler characteristic of the sphere."
        v = vertex.count_instances()   # vertices
        e = edge.count_instances()     # edges
        f = triangle.count_instances() # faces
        return v - e + f
    
    def curvature_total(self,normalized=False):
        """
        Calculates the total Gauss curvature of the sphere. If
        normalized=True, also divide by the number of vertices to get
        an average.
        """
        # Sum up the curvature over all points
        integrated_curvature = 0
        for point in vertex.instances.values():
            integrated_curvature += point.curvature()
        # Maybe normalize
        if normalized:
            integrated_curvature /= float(vertex.count_instances())
        return integrated_curvature

    def curvature_std(self):
        """
        Calculates standard deviation of the Gauss curvature over the
        sphere. 
        """
        # Make a list of all local curvatures
        local_curvatures = [point.curvature() for point in vertex.instances.values()]
        # Get the standard deviation
        return np.std(local_curvatures)
          
    def __str__(self):
        "The state of the system at a given time."
        outstring = """Sphere Current state:
---------------------------------TOPOLOGY----------------------------------
Number of Vertices:   {}
Number of Edges:      {}
Number of Triangles:  {}
Euler Characteristic: {}
---------------------------------CURVATURE---------------------------------
Total: {}
Mean:  {}
Std:   {}
""".format(vertex.count_instances(),
           edge.count_instances(),
           triangle.count_instances(),
           self.euler_characteristic(),
           self.curvature_total(),
           self.curvature_total(True),
           self.curvature_std())
        return outstring
#---------------------------------------------------------------------------


### Functions that generate or remove objects en bulk
###-------------------------------------------------------------------------
def make_n_vertices(n):
    """
    Makes n new vertices and adds them to the hash table. Returns the
    vertex IDS generated.
    """
    new_vertices = []
    for i in range(n):
        temp = vertex()
        vertex.add(temp)
        new_vertices.append(temp.id)
    return new_vertices

def delete_all_geometries():
    "Delete all instances of all geometric objects."
    vertex.delete_all()
    edge.delete_all()
    triangle.delete_all()
###-------------------------------------------------------------------------


### Functions that manipulate objects at a higher level. These are
### basically syntactic sugar.
#-------------------------------------------------------------------------
def bisect_edge(edge_id):
    """
    Uses the edge.bisect() function to bisect an edge and update the
    hash tables properly. Takes an edge ID.
    """
    # Bisect the edge and generate two pairs of points, each
    # representing one of the two new edges.
    new_edges,newpoint = edge.instances[edge_id].bisect() 

    # Delete the old edge
    edge.delete(edge_id)

    # Add the new point to the hash table
    vertex.instances[newpoint] = vertex([],newpoint)
    
    # Add the new edges to the edges hash table. We want to return the
    # new idge IDs (because why not), so a temporary list is generated
    # for this purpose.
    new_ids = []
    for endpoint_pair in new_edges:
        e = edge(endpoint_pair)
        edge.instances[e.id] = e
        new_ids.append(e.id)
    return new_ids
# -------------------------------------------------------------------------



### Functions that make simplices and sub-simplices together. The
### workhorses of state manipulation.
### -------------------------------------------------------------------------
def build_triangle(point_list):
    """
    Builds a triangle by creating the triangle object (if it doesn't
    already exist) and all subsimplices---vertices and edges---if they
    don't already exist. Checks for redundancies.
    """
### -------------------------------------------------------------------------
