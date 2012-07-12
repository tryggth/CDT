#!/usr/bin/env python2

"""
simplex.py 

Author: Jonah Miller (jonah.miller@colorado.edu) 

This file is part of the sphere_generator program which generates
spheres of (as close as possible to) uniform curvature for a given
surface area by monte carlo simulation.

Data structures and classes are contained in this file. We need a
class to contain information about triangle vertices, triangle edges,
and triangles. Each one class instance will be stored in a dictionary
with a value to access it given by an ID number, which we cycle
through and recycle. There is a shared "static" attribute for each
class to keep track of this.
"""



### Dependencies
#-------------------------------------------------------------------------
import numpy as np
import scipy as sp
from copy import * # For copy and deepcopy
#-------------------------------------------------------------------------



### Global Constants Relating to Classes. On the top for convenience only.
#-------------------------------------------------------------------------
edge_length = 1 # For triangles. Triangles are equilateral.

# area = (1/2) * base * height
# base = a
# height^2 = hypotenuse^2 - (base/2)^2 = sqrt(3)(edge_length)/2
triangle_area = np.sqrt(3) * edge_length**2 / 4
# Triangle Dihedral Angle. For calculating curvature.
triangle_angle = np.pi/3 
#-------------------------------------------------------------------------



### Dictionaries that contain collections of class objects. The keys
### are ids. The values are class instances. These hold all
### vertices,edges,points that exist on the surface of the sphere at
### any given time.
#-------------------------------------------------------------------------
vertices = {}  # For points
edges = {}     # For edges
triangles = {} # For triangles

# So that certain functions (delete, for instance) can access these
# hash tables using strings as input, we have this dictionary. It uses
# the same keys as the object_map dictionary (below).
hash_tables_reference = {'vertex':vertices,
                         'edge':edges,
                         'triangle':triangles}
#-------------------------------------------------------------------------



### Functions that act on all classes
###-------------------------------------------------------------------------
def increment_id(class_name):
    """Increment the staticvariable last_used_id for a class by 1. 
    Then return the new id for use."""
    # Accesses the class by using a string.
    object_map[class_name].last_used_id += 1
    return object_map[class_name].last_used_id

# Take an unused ID name and add it to the list of recycled IDs.
def reclaim_id(class_name,object_id):
    """"Take an input class ID and add it to the list of recycled IDs for
    that class."""
    object_map[class_name].recycled_ids.append(object_id)
    # Not strictly necessary to have a return value, but why not?
    return object_map[class_name].recycled_ids 

# Take an element of recycled ids for a given class, pop it out, and
# return it.
def recycle_id(class_name):
    """Looks for the minimum element of your recycled IDs of a given object, 
    returns it, and removes it from the list."""
    # Find the minimum element
    minval = min(object_map[class_name].recycled_ids)  
    # Remove the minimum element from the list.
    object_map[class_name].recycled_ids.remove(minval) 
    return minval # Return the returned element for reuse.

# Generate an ID name (new or recycled, depending on availability)
# for use.
def make_id(class_name):
    """Generate an ID name for use with the given class. First looks in 
    recycled ids. If none are available, makes a brand new ID."""
    if len(object_map[class_name].recycled_ids) > 0:
        new_id = recycle_id(class_name)
    else:
        new_id = increment_id(class_name)
    return new_id

# Delete an object from the correct hash table and reclaim its id.
def delete(class_name, object_id):
    """
    Delete an object from the correct hash table and reclaim its id.
    """
    # Delete from the hash table. Finds the correct hash
    # table/dictionary containing that object type in
    # hash_tables_reference, and then deletes the object id from the
    # correct hash table.
    del hash_tables_reference[class_name][object_id]
    # Reclaims the id:
    reclaim_id(class_name,object_id)
###-------------------------------------------------------------------------




# Classes
#-------------------------------------------------------------------------
# The points class
class vertex:
    """
    A class for keeping track of points.
    
    Attributes: 
    -- last_used_id = the last used point ID. global/static. 
                      Shared by all class instances
    -- recycled_ids = A list of used IDs, which, after deletion, become 
                      available. Global/static. Shared by all class instances.
    -- triangles    = A list of triangles that contain this vertex. Used to 
                      calculate curvature. Defined by the init function
    -- curvature() = a function that calculates the local curvature
                     of the manifold at a vertex.
    -- edges() = A function that calculates the edges that contain the
                 vertex. Uses triangle properties. (NOT YET IMPLIMENTED)

    Functions that act on the class:
    -- increment_id(class_name) = A function that increments last_used_id
    -- reclaim_id(class_name,id) = Takes the id and adds it to recycled_ids 
                                 list.
    -- recycle_id(class_name) = A function that returns the least element of 
                               recycled_ids and deletes it from the list.
    -- make_id(class_name) = A function that generates a new id for a vertex.
                            First looks in recycled IDs. If there are non
                            available, uses increment ID.

    To initialize: v = vertex()
    """               
    
    # Static/global variables
    last_used_id = 0  # Last used point ID.
    recycled_ids = [] # List of unused point IDS <= last_used_id
    
    # Initialization function
    def __init__(self, triangle_list = []):
        """
        Initialize a vertex. You can give the vertex a list of
        triangles you know are connected to it to add them to the list
        of triangles. However, by default, this list is empty.
        """
        # Local triangle list
        self.triangles = copy(triangle_list)
    
    # Other functions
    def curvature(self):
        """Curvature is directly proportional to the deficit angle
        2*pi - area_of_triangle * sum(triangles attached to vertex)."""
        return 2 * (2 * np.pi - triangle_angle * len(self.triangles))

# The edge class
class edge:
    """
    A class for keeping track of triangle edges
    
    Attributes: 
    -- last_used_id = the last used point ID. global/static. 
                      Shared by all class instances
    -- recycled_ids = A list of used IDs, which, after deletion, become 
                      available. Global/static. Shared by all class instances.
    -- vertices     = A list of point IDs for edge endpoints.
    -- bisect       = Bisects the edge and returns two pairs of points, each 
                      defining a new edge.

    Functions that act on the class:
    -- increment_id(class_name) = A function that increments last_used_id
    -- reclaim_id(class_name,id) = Takes the id and adds it to recycled_ids 
                                 list.
    -- recycle_id(class_name) = A function that returns the least element of 
                               recycled_ids and deletes it from the list.
    -- make_id(class_name) = A function that generates a new id for a vertex.
                            First looks in recycled IDs. If there are non
                            available, uses increment ID.
    -- points = a function that returns the points the edge contains
    
    Example of call: e = edge(vertex_pair)
    where vertex_pair = [vertex1, vertex2]
    """

    # Static/global variables
    last_used_id = 0 # Last used edge ID.
    recycled_ids = [] # List of unused edge IDS <= last_used_id

    # Functions
    def bisect(self):
        """
        If the edge length of the triangle is 
        """
        if len(self.vertices) != 2:
            print "Your edge doesn't have the right number of points!"
            return
        newpoint = make_id('vertex')
        return [[self.vertices[0],newpoint],[newpoint,self.vertices[1]]]


    # Initialization function
    def __init__(self,vertex_pair):
        """
        Initialize an edge. Takes a vertex pair (list) as input.  If
        you attempt to initialize an edge with the wrong number of
        endpoints, an edge with no endpoints is generated and an error
        message is given.
        """
        # Points contained in edge (IDS)
        if len(vertex_pair) != 2:
            print "Your edge doesn't have the right number of points!"
            print "Generating edge with no endpoints."
            self.vertices = []
        else:
            self.vertices = copy(vertex_pair)


# The triangle class
class triangle:
    """
    A class for keeping track of triangle edges
    
    Attributes: 
    -- last_used_id = the last used point ID. global/static. 
                      Shared by all class instances
    -- recycled_ids = A list of used IDs, which, after deletion, become 
                      available. Global/static. Shared by all class instances.
    -- points       = A list of points contained by the triangle

    Functions that act on the class:
    -- increment_id(class_name) = A function that increments last_used_id
    -- reclaim_id(class_name,id) = Takes the id and adds it to recycled_ids 
                                 list.
    -- recycle_id(class_name) = A function that returns the least element of 
                               recycled_ids and deletes it from the list.
    -- make_id(class_name) = A function that generates a new id for a vertex.
                            First looks in recycled IDs. If there are non
                            available, uses increment ID.
    -- edges() = A list of edges contained by the triangle

    
    Example of call: t = triangle(point_list)
    where point_list = [p1,p2,p3]
    """

    # Static/global variables
    last_used_id = 0 # Last used triangle ID.
    recycled_ids = [] # List of unused vertex IDS <= last_used_id

    # Initialization function
    def __init__(self):
        "Initialize a triangle"
        self.points = [] # Vertexes of the triangle (IDS)
        self.edges = [] # Edges of the triangle (IDS)


# A map from strings to class names. Used for all functions in this section.
object_map = {'vertex':vertex,
              'edge':edge,
              'triangle':triangle}
#-------------------------------------------------------------------------



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
    new_edges = edges[edge_id].bisect() 

    # Delete the old edge
    delete('edge',edge_id)
    
    # Add the new edges to the edges hash table. We want to return the
    # new idge IDs (because why not), so a temporary list is generated
    # for this purpose.
    new_ids = []
    for endpoint_pair in new_edges:
        new_id = make_id('edge')
        edges[new_id] = edge(endpoint_pair)
        new_ids.append(new_id)
    return new_ids

# -------------------------------------------------------------------------

