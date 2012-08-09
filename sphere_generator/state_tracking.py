#!/usr/bin/env python2

"""
state_tracking.py

Author: Jonah Miller (jonah.maxwell.miller@gmail.com)

This file is part of the sphere_generator program which generates
spheres of (as close as possible to) uniform curvature for a given
surface area by monte carlo simulation.

This file contains classes and functions for measuring properties of
the generated sphere.
"""


### Dependencies
#-------------------------------------------------------------------------
import numpy as np
import scipy as sp
# Class data structures we need
import simplex_ancestors as sa
import simplex_descendants as sd
import state_manipulation as sm
import utilities as ut
import error_checking
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
        v = sd.vertex.count_instances()   # vertices
        e = sd.edge.count_instances()     # edges
        f = sd.triangle.count_instances() # faces
        return v - e + f
    
    def curvature_total(self,normalized=False):
        """
        Calculates the total Gauss curvature of the sphere. If
        normalized=True, also divide by the number of vertices to get
        an average.
        """
        # Sum up the curvature over all points
        integrated_curvature = 0
        for point in sd.vertex.instances.values():
            integrated_curvature += point.curvature()
        # Maybe normalize
        if normalized:
            integrated_curvature /= float(sd.vertex.count_instances())
        return integrated_curvature

    def curvature_std(self):
        """
        Calculates standard deviation of the Gauss curvature over the
        sphere. 
        """
        # Make a list of all local curvatures
        local_curvatures = [point.curvature() for \
                                point in sd.vertex.instances.values()]
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
           sd.edge.count_instances(),
           sd.triangle.count_instances(),
           self.euler_characteristic(),
           self.curvature_total(),
           self.curvature_total(True),
           self.curvature_std())
        return outstring
#---------------------------------------------------------------------------

