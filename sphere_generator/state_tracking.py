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
""".format(sd.vertex.count_instances(),
           sd.edge.count_instances(),
           sd.triangle.count_instances(),
           self.euler_characteristic(),
           self.curvature_total(),
           self.curvature_total(True),
           self.curvature_std())
        return outstring
#---------------------------------------------------------------------------


# The imaginary vertex class, which contains information used for move
# attempts
# ---------------------------------------------------------------------------
class imaginary_vertex:
    """
    A class used by try-move functions to calculate the "action" which
    determines whether or not to accept a move. Tells whether a vertex
    is being added or removed from the spacetime, how many triangles
    it has, and what the curvature around it would be if it were real.
    """

    def __init__(self,num_triangles,volume_increasing):
        """
        Initializes an imaginary vertex. Tell it the number of
        triangles that contain it, and whether its a volume increasing
        move or a volume decreasing move. For volume increasing, pass
        True. Otherwise, pass False.  
        """

        self.num_triangles = int(num_triangles)
        self.volume_increasing = bool(volume_increasing)

    
    def __str__(self):
        "What string typecasting for an imaginary vertex reveals."
        return str(self.curvature())

    def __len__(self):
        "Length reveals the number of attached triangles."
        return int(self.num_triangles)

    def curvature(self):
        "Curvature is directly proportional to deficit angle."
        return 2 * (2 * np.pi - sd.triangle.angle * self.num_triangles)
# ---------------------------------------------------------------------------



# The move_data class contains many instances of imaginary vertices
# and knows how to calculate expected curvature mean and standard
# deviation."
# ---------------------------------------------------------------------------
class move_data:
    """
    The move_data class contains many instances of imaginary vertices
    and knows how to calculate expected curvature mean and standard
    deviation.  
    """

    def __init__(self,i_vertex_list):
        """
        Set the imaginary vertices. Take a list of imaginary vertices
        as input.
        """
        self.imaginary_vertices = i_vertex_list

    def __len__(self):
        "Just returns the number of imaginary vertices in the move data."
        return len(self.imaginary_vertices)

    def __str__(self):
        "Prints info on the imaginary vertices and calculation results."
        # Initialize string
        outstring = 'Volume Increasing Vertices:\n'
        
        # Vertices for volume increasing and volume decreasing
        v_increasing = [v for v in self.imaginary_vertices \
                            if v.volume_increasing]
        v_decreasing = [v for v in self.imaginary_vertices \
                            if not v.volume_increasing]
        for v in v_increasing:
            outstring += "Curvature: {}.\n".format(v.curvature())
        outstring += "Volume Decreasing Vertices:\n"
        for v in v_decreasing:
            outstring += "Curvature: {}.\n".format(v.curvature())
        
        outstring += "\nExpected Mean: {}\n".format(self.predicted_mean_curvature())
        outstring += "Expeced std dev: {}\n".format(self.predicted_curvature_std_dev())
    
        return outstring

    def predicted_mean_curvature(self):
        """
        Predict the mean curvature by calculating total curvature and
        subtracting off curvature from volume decreasing moves and
        adding curvature from volume increasing simplices. Divide by
        the "expected" number of vertices.
        """
        curvature_total = 0 # Initialize curvature total
        # Initialize total number of vertices
        total_vertices = sd.vertex.count_instances()

        # Get total from real vertices
        for vertex in sd.vertex.instances.values():
            curvature_total += vertex.curvature()

        # Get total from imaginary vertices
        for imaginary_vertex in self.imaginary_vertices:
            if imaginary_vertex.volume_increasing:
                total_vertices += 1 # More vertices
                curvature_total += imaginary_vertex.curvature()
            else:
                total_vertices -= 1 # Fewer vertices
                curvature_total -= imaginary_vertex.curvature()
            # If we have fewer than 0 vertices, something went VERY wrong.
            assert total_vertices >= 0 

        # Mean is expected total curvature divided by expected total vertices
        return curvature_total/float(total_vertices)

    def predicted_curvature_std_dev(self):
        """
        Predict the standard deviation of curvature. This is more
        complicated than predicted_mean_curvature, but the idea is the
        same.
        """

        # Need the mean to calculate the standard deviation
        sample_mean = self.predicted_mean_curvature()

        # For convenience, we define the "deviation" of a single point
        # from the mean.
        dev = lambda x: (x-sample_mean)**2

        # First get the total "deviation" from real vertices
        total_dev = sum([dev(vertex.curvature()) \
                             for vertex in sd.vertex.instances.values()])
        # Initialize total vertices (will be changed)
        total_vertices = sd.vertex.count_instances()

        # Now we need to account for imaginary vertices
        for imaginary_vertex in self.imaginary_vertices:
            if imaginary_vertex.volume_increasing:
                total_vertices += 1 # More vertices
                total_dev += dev(imaginary_vertex.curvature())
            else:
                total_vertices -= 1 # Fewer vertices
                total_dev -= dev(imaginary_vertex.curvature())
            # If we have fewer than 0 vertices something went VERY WRONG.
            assert total_vertices >= 0
            # If the standard deviation is less than zero, the move is
            # not topologically acceptable.
            if total_dev < 0:
                return False

        # Now, calculate the standard deviation as the square root of
        # the average of the deviations:
        ave_dev = total_dev / total_vertices
        return np.sqrt(ave_dev)
# ---------------------------------------------------------------------------
