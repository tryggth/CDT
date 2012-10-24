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

    
# Vertex Count class
# ---------------------------------------------------------------------------
class vertex_count():
    """
    The vertex count class holds information about vertices. You call
    the constructor on a list of vertex objects (or IDS), and you can
    then ask it how many vertices are of order n (how many triangles
    are attached to them.
    """
    def __init__(self, vertex_list=[], imaginary_vertex_list = []):
        """
        Constructs a vertex_count instance. Calculates the number of
        vertices in the list vertex_list of order n.
        """
        # The container for vertex order information. The key is the
        # order, the value is the number of vertices.
        self.orders = {}
        # Now, we calculate the vertex orders. If the vertex list is
        # empty, nothing to do.
        if vertex_list:
            #If we have a list of imaginary vertices, take those into
            # account.
            if imaginary_vertex_list:
                self.calculate_from_imaginary_vertices(vertex_list,
                                                       imaginary_vertex_list)
            else: # Otherwise, just use the real vertices.
                self.calculate_vertex_orders(vertex_list)

    def reset_orders(self):
        """
        Empties the orders field. Useful if we want to recycle this object.
        """
        self.orders = {}

    def calculate_vertex_orders(self, vertex_list):
        """
        Goes through the list of vertex objects or ids and, if a
        vertex is of a given order, add it to the count of vertices of
        that order.
        """
        # We want objects for our vertices, so
        vertex_object_list = [sd.vertex.parse_input(v) for v in vertex_list]
        self.reset_orders()
        for v in vertex_object_list:
            if not len(v) in self.orders.keys():
                self.orders[len(v)] = 1
            else:
                self.orders[len(v)] += 1
        self.enforce_vertex_orders(vertex_object_list)

    def check_vertex_orders(self,vertex_object_list):
        """
        Check that the total number of vertices counted in orders is
        the same as the number of vertices in
        vertex_object_list. Returns a boolean
        """
        return sum(self.orders.values()) == len(vertex_object_list)

    def enforce_vertex_orders(self,vertex_object_list):
        """
        Checks that the total number of vertices counted in orders is
        the same as the total number of vertex_object_list. Raises an
        error if this is not the case.
        """
        if not self.check_vertex_orders(vertex_object_list):
            print "The number of triangles counted in vertex_object_list " +\
                "do not match the number in self.orders."
            print "Object list "+vertex_object_list
            print "Orders "+self.orders
            raise ValueError("self.order not the correct length.")

    def get_vertex_count_for_order(self,order):
        """
        Returns the number of vertices of the specified order.
        """
        if order in self.orders.keys():
            return self.orders[order]
        else:
            return 0

    def get_total_vertex_count(self):
        """
        Returns the the total number of vertices we keep track of.
        """
        return sum(self.orders.values())

    # The following method is useful for working with imaginary
    # vertices. Not in the constructor to avoid buggy behavior.
    def calculate_from_imaginary_vertices(self,real_vertex_list,
                                          imaginary_vertex_list):
        """
        Calculates the vertex count using both real and imaginary
        vertices. Useful for working with move data.
        """
        self.calculate_vertex_orders(real_vertex_list)
        for v in imaginary_vertex_list:
            num_triangles = len(v)
            if v.volume_increasing:
                if num_triangles in self.orders.keys() \
                        and self.orders[num_triangles] > 0:
                    self.orders[num_triangles] -= 1
            else:
                if num_triangles in self.orders.keys():
                    self.orders[num_triangles] += 1
                else:
                    self.orders[num_triangles] = 1
    
    # Data conventions methods
    def __len__(self):
        """
        The length of the vertex count object.
        """
        return len(self.orders)

    def __repr__(self):
        """
        The debugging form of a vertex_count object.
        """
        return str(self.orders)

    def __str__(self):
        """
        Return a useful output for computation.
        """
        output = "#Order\tn-vertices\n"
        for i in self.orders.keys():
            output+="{}\t{}\n".format(i, self.orders[i])
        return output

    def __eq__(self,other):
        """
        Tests whether or not two sets of vertex counts are the
        same. Returns true if they give the same vertex counts for
        each order.
        """
        return self.orders == other.orders

    def __ne__(self,other):
        """
        Not equal condition.
        """
        return not self.__eq__(other)

#---------------------------------------------------------------------------

# Subclasses of vertex_count
#---------------------------------------------------------------------------
class vertex_count_selection_optimal(vertex_count):
    """
    This subclass of vertex_count has a fitness function for finding
    out how close to optimal a sphere is. It also has a boolean
    function for deciding when a sphere is close enough to optimal.

    Here, optimal means that there are N-12 vertices of order 6 and 12
    vertices of order 5, where N is the total number of vertices.

    "close enough" means, for some damping integers chosen at the
    start, call them D1 and D2, let V5 be the vertices of order 5 and
    V6 be the vertices of order 6. Then,

    N-12-D1 <= V6 <= N-12+D1 and 12-D2 <= V5 <= 12+D2
    """
    def __init__(self, order_5_damping, order_6_damping,
                 vertex_list = [],
                 imaginary_vertex_list = []):
        """
        Same as super(self).__init__, except sets the damping integers
        for order 5 and order 6 damping. Enforces that they're
        positive integers. If they aren't typecasts them.
        """
        vertex_count.__init__(self,vertex_list,imaginary_vertex_list)
        self.order_5_damping = abs(int(order_5_damping))
        self.order_6_damping = abs(int(order_6_damping))


    def optimal_count_order_6(self):
        """
        We want N-12 vertices of order 6, where N is the total number
        of vertices.
        """
        return self.get_total_vertex_count() - 12

    def optimal_count_order_5(self):
        """
        We want 12 vertices of order 5. This is only a function in
        case we want to change these criteria later.
        """
        return 12

    def fitness_order_6(self):
        """
        How close we are to the optimal count for order 6.
        """
        return abs(self.optimal_count_order_6() \
                       - self.get_vertex_count_for_order(6))

    def fitness_order_5(self):
        """
        How close we are to the optimal count for order 5.
        """
        return abs(self.optimal_count_order_5() \
                       - self.get_vertex_count_for_order(5))

    def fitness_function(self):
        """
        A fitness function similar to those used for the metropolis
        algorithm. Possibly useful if we decide to select for
        this. Gives a quick idea how close we are at a clance.
        """
        return np.exp(-(self.fitness_order_5()+self.fitness_order_6()))

    def is_close_enough(self):
        """
        Tests whether or not a given sphere is close enough to
        microscopically optimal. Let V5 be the total number of
        vertices of order 5, V6 the total number of vertices of order
        6, and N the total number of vertices. Then we test whether or
        not

        N-12-D1 <= V6 <= N-12+D1 and 12-D2 <= V5 <= 12+D2

        Returns a boolean.
        """
        return self.fitness_order_5() <= self.order_5_damping\
            and self.fitness_order_6() <= self.order_6_damping
    
#---------------------------------------------------------------------------



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
