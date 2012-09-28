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
          
    def surface_area(self):
        """
        Calculates the surface area of the sphere. Basically syntactic sugar.
        """
        return sd.triangle.count_instances()

    def get_vertices(self):
        "Returns a string each vertex in the sphere. Output could be long."
        outstring = ''
        for v in sd.vertex.instances.values():
            outstring += str(v) + '\n'
        return outstring

    def get_edges(self):
        "Returns a string each edge in the sphere. Output could be long."
        outstring = ''
        for e in sd.edge.instances.values():
            outstring += str(e) + '\n'
        return outstring

    def get_triangles(self):
        "Returns a string each triangle in the sphere. Output could be long."
        outstring = ''
        for t in sd.triangle.instances.values():
            outstring += str(t) + '\n'
        return outstring
        
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

    def __init__(self,i_vertex_list,cmpx,move_type,change_in_surface_area=0):
        """
        Set the imaginary vertices. Take a list of imaginary vertices
        as input.
        """
        self.imaginary_vertices = i_vertex_list
        self.complex = cmpx
        self.move_type = move_type
        self.change_in_surface_area = change_in_surface_area

    def __len__(self):
        "Just returns the number of imaginary vertices in the move data."
        return len(self.imaginary_vertices)

    def __str__(self):
        "Prints info on the imaginary vertices and calculation results."
        outstring = 'Move type: {}\n'.format(self.move_type)

        # Initialize string
        outstring += 'Volume Increasing Vertices:\n'
        
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
        outstring += "Expected std dev: {}\n".format(self.predicted_curvature_std_dev())
    
        return outstring

    def get_move_type(self):
        return self.move_type

    def get_complex(self):
        return self.complex

    def predicted_mean_curvature(self):
        """
        Predict the mean curvature by calculating total curvature and
        subtracting off curvature from volume decreasing moves and
        adding curvature from volume increasing simplices. Divide by
        the "expected" number of vertices.
        """
        # Initialize total number of vertices
        total_vertices = sd.vertex.count_instances()

        # Get total from real vertices
        curvature_total = sum([v.curvature() \
                                   for v in sd.vertex.instances.values()])

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
        # A finite acceptable negative number to account for
        # discritization error.
        acceptable_min = -0.01

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

        # Try to account for rounding error
        total_dev = ut.round_to_zero(total_dev)

        # If the standard deviation is less than zero, the move is
        # not topologically acceptable.
        if total_dev < 0:
            return False


        # Now, calculate the standard deviation as the square root of
        # the average of the deviations:
        ave_dev = total_dev / total_vertices
        return np.sqrt(np.abs(ave_dev))

    def predicted_surface_area(self):
        "Predicts the number of triangles on the sphere."
        # The triangle number as it is now
        triangle_number = sd.triangle.count_instances()

        # Now we need to account for imaginary vertices
        triangle_number += self.change_in_surface_area

        # If we have fewer than 4 triangles, something went wrong
        assert triangle_number >= 4

        return triangle_number
        
# ---------------------------------------------------------------------------


# Functions that look at local properties:
# ---------------------------------------------------------------------------
def find_opposite_vertices(triangle_set):
    """ 
    A function that looks for a pair of vertices in a set of triangles
    that don't share a triangle. Returns all such vertex pairs. Takes
    ids or instances as input. But only takes collections.
    """
    # Parse input
    triangle_objects = set([sd.triangle.parse_input(t) for t in triangle_set])
    # Extract vertices
    vertices = ut.set_union([set(t.get_vertices()) for t in triangle_objects])
    # Finds pairs of vertices that don't share a triangle
    opposite_vertices = [set([v1,v2]) for v1 in vertices \
                             for v2 in vertices \
                             if not v1.shares_a_triangle_with(v2)]
    # We need to remove duplicates from the list of opposite vertices
    for vertex_pair in opposite_vertices:
        if vertex_pair not in filtered_vertices:
            filtered_vertices.append(vertex_pair)

    return filtered_vertices

def shared_triangles_in_complex(vertex1,vertex2,triangle_complex):
    """
    Meant for internal use with
    find_opposite_vertices_in_complex. Given a pair of vertices,
    finds if they are both vertices of the same triangle in a
    complex. The complex must be a collection of objects or ids.
    """
    # Parse input
    triangles = set([sd.triangle.parse_input(t) for t in triangle_complex])
    v1 = sd.vertex.parse_input(vertex1)
    v2 = sd.vertex.parse_input(vertex2)
    
    # Shared triangles for the vertex pair that intersect with the complex
    shared_triangles = set(v1.triangles_shared_with(v2)) & triangles
    
    return shared_triangles

def shares_triangles_in_complex(vertex1,vertex2,triangle_complex):
    "Like shared_triangles_in_complex, but returns a boolean."
    return bool(shared_triangles_in_complex(vertex1,vertex2,triangle_complex))

def find_opposite_vertices_in_complex(triangle_complex):
    """
    Similar to find_opposite_vertices, but looks vertices that don't
    share a triangle contained in the triangle complex.
    """
    # Parse input
    triangle_objects = set([sd.triangle.parse_input(t) \
                                for t in triangle_complex])
    # Extract vertices
    vertices = ut.set_union([set(t.get_vertices()) for t in triangle_objects])
    # Vertices that don't share triangles are:
    shared = [set([v1,v2]) for v1 in vertices for v2 in vertices\
                  if not shares_triangles_in_complex(v1,v2,triangle_objects)]
    # Filter out duplicates
    filtered = []
    for vertex_pair in shared:
        if vertex_pair not in filtered:
            filtered.append(vertex_pair)

    return filtered

# ---------------------------------------------------------------------------
