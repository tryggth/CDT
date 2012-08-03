#!/usr/bin/env python2

"""
simplex_descendants.py

Author: Jonah Miller (jonah.maxwell.miller@gmail.com)

This file is part of the sphere_generator program which generates
spheres of (as close as possible to) uniform curvature for a given
surface area by monte carlo simulation.

Data structures and classes are contained in this file. We need a
class to contain information about triangle vertices, triangle edges,
and triangles. Each one class instance will be stored in a dictionary
with a value to access it given by an ID number, which we cycle
through and recycle. There is a shared "static" attribute for each
class to keep track of this.

The functions that manipulate dictionaries are part of the goemetry
class, which is in the simplex_ancestors.py
module.
"""



### Dependencies
#-------------------------------------------------------------------------
import numpy as np
import scipy as sp
from simplex_ancestors import * # Required for parent classes
#-------------------------------------------------------------------------



####---------------------------------------------------------------------####
#                               Classes                                     #
####---------------------------------------------------------------------####
# The points class
#-------------------------------------------------------------------------
class vertex(geometry):
    """
    A class for keeping track of points.
    
    Attributes: 
    -- last_used_id = the last used point ID. global/static. 
                      Shared by all class instances
    -- recycled_ids = A list of used IDs, which, after deletion, become 
                      available. Global/static. Shared by all class instances.
    -- triangles    = A list of triangles that contain this vertex. Used to 
                      calculate curvature. Defined by the init function
    -- id = The identifying number of the simplex. Useful. 
    -- instances = The dictionary containing all instances of the class.

    Methods:
    -- find_duplicates(triangle_list) = A function that searches for other 
                                        vertices that have the same triangle
                                        list as triangle_list.
    -- find_triangles() = A function that calculates triangles that contain 
                          this vertex.
    -- set_triangles(triangle_list) = A way to set self.triangles
    -- find_and_set_triangles() = A function that uses find triangles and 
                                  set triangles to do exactly that.
    -- curvature() = a function that calculates the local curvature
                     of the manifold at a vertex.
    -- connect_surrounding_triangles() = Ensures surrounding triangles 
                                         see each other as neighbors.
    -- find_edges() = A function that calculates the edges that contain the
                     vertex.
    -- increment_id() = A function that increments last_used_id
    -- reclaim_id(id) = Takes the id and adds it to recycled_ids 
                                 list.
    -- recycle_id() = A function that returns the least element of 
                     recycled_ids and deletes it from the list.
    -- make_id() = A function that generates a new id for a vertex.
                   First looks in recycled IDs. If there are non
                   available, uses increment ID.


    To initialize: v = vertex()
    """               
    
    # Static/global variables
    last_used_id = 0  # Last used point ID.
    recycled_ids = set([]) # List of unused point IDS <= last_used_id
    instances = {} # Contains all instances of the class

    # Duplicate checking
    @classmethod
    def find_duplicates(self,triangle_list):
        """
        Searches through self.instances for vertices with the same
        triangles as in triangle_list. Requires set equality. Returns
        the duplicates.
        """
        duplicates = []
        for v in self.instances.values():
            if v.triangles == set(triangle_list):
                duplicates.append(v.id)
        return duplicates


    # Initialization function
    def __init__(self, triangle_list = [],new_id=0):
        """
        Initialize a vertex. You can give the vertex a list of
        triangles you know are connected to it to add them to the list
        of triangles. However, by default, this list is empty.
        """
        # ID number
        if new_id == 0:
            self.id = self.make_id()
        else:
            self.id = new_id
        # Local triangle list
        self.triangles = set(triangle_list)

    def __str__(self):
        "String conversion reveals id and the number of triangles."
        t_string = str([x for x in self.triangles])
        idstring = str(self.id)
        return "Vertex object.\nID: {}\nTriangles: {}".format(idstring, 
                                                                t_string)

    def __len__(self):
        "Length reveals the number of attached triangles."
        return len(self.triangles)

    
    # Other functions
    def find_edges(self):
        "Calculates what edges contain this vertex."
        containing_edges = set([])
        for e in edge.instances.values():
            if self.id in e.vertices:
                containing_edges.add(e.id)
        return containing_edges

    def find_triangles(self):
        """
        Calculates what triangles contain this vertex. Can use to set
        triangles.
        """
        containing_triangles = set([])
        for t in triangle.instances.values():
            if self.id in t.vertices:
                containing_triangles.add(t.id)
        return containing_triangles

    def set_triangles(self,triangle_list):
        "Sets self.triangles"
        self.triangles = set(triangle_list)

    def find_and_set_triangles(self):
        """
        Finds attached triangles and sets self.triangles to include them.
        """
        new_triangles = self.find_triangles()
        self.set_triangles(new_triangles)
        return new_triangles

    def curvature(self):
        """Curvature is directly proportional to the deficit angle
        2*pi - area_of_triangle * sum(triangles attached to vertex)."""
        return 2 * (2 * np.pi - triangle.angle * len(self.triangles))

    def connect_surrounding_triangles(self):
        """
        For each triangle around the point, check their
        intersections to find which triangles share an edge.
        Only works properly once simplex_countainers.py is imported.
        """
        for id1 in self.triangles:
            for id2 in self.triangles:
                triangle.instances[id1].connect_to_triangle(triangle.instances[id2])
#-------------------------------------------------------------------------


# The edge class
#-------------------------------------------------------------------------
class edge(geometry):
    """
    A class for keeping track of triangle edges
    
    Attributes: 
    -- last_used_id = the last used point ID. global/static. 
                      Shared by all class instances
    -- recycled_ids = A list of used IDs, which, after deletion, become 
                      available. Global/static. Shared by all class instances.
    -- vertices     = A list of point IDs for edge endpoints.
    -- points = a function that returns the points the edge contains
    -- id = The identifying number of the simplex. Useful. 
    -- length = The length of an edge in the physical geometry.
    -- instances = The dictionary containing all instances of the class.

    Methods:
    -- bisect()  = Bisects the edge and returns two pairs of points, each 
                   defining a new edge.
    -- increment_id() = A function that increments last_used_id
    -- reclaim_id(id) = Takes the id and adds it to recycled_ids 
                        list.
    -- recycle_id() = A function that returns the least element of 
                      recycled_ids and deletes it from the list.
    -- make_id() = A function that generates a new id for a vertex.
                   First looks in recycled IDs. If there are non
                   available, uses increment ID.
    -- check_topology() = A function to make sure the edge has 
                          either 2 or 0 endpoints.
    -- set_vertices(vertices) = A function to set self.vertices. Takes a 
                                list of length 2 as input.

    Example of call: e = edge(vertex_pair)
    where vertex_pair = [vertex1, vertex2]
    """

    # Static/global variables
    last_used_id = 0 # Last used edge ID.
    recycled_ids = set([]) # List of unused edge IDS <= last_used_id
    instances = {} # Contains all instances of the class

    # Constant value length
    length = 1

    # Functions
    def bisect(self):
        """
        If the edge length of the triangle is equilateral, cuts it in
        two and generates two edges with a point between them. Returns
        the edges and the new point ID generated.
        """
        if len(self.vertices) != 2:
            print "Your edge doesn't have the right number of points!"
            return
        newpoint = vertex.make_id()
        new_edges = [set([newpoint,point]) for point in self.vertices]
        return [new_edges,newpoint]

    def check_topology(self):
        """
        Check to ensure that the number of endpoints is correct.
        """
        assert len(self.vertices) == 2 or len(self.vertices) == 0
        print "Topology is okay."

    def set_vertices(self,new_vertices):
        "Takes an input list or set and resets self's vertices/endpoints."
        if len(new_vertices) == 2 or len(new_vertices) == 0:
            self.vertices = set(new_vertices)
            return self.vertices
        else:
            print "Error! Wrong number of vertices! Nothing changed."
            return self.vertices

    # Duplicate checking
    @classmethod
    def find_duplicates(self,vertex_list):
        """
        Searches through self.instances for edges with the same
        vertices as in triangle_list. Requires set equality. Returns
        the duplicates.
        """
        duplicates = []
        for e in self.instances.values():
            if e.vertices == set(vertex_list):
                duplicates.append(e.id)
        return duplicates


    # Initialization function
    def __init__(self,vertex_pair=[]):
        """
        Initialize an edge. Takes a vertex pair (list) as input.  If
        you attempt to initialize an edge with the wrong number of
        endpoints, an edge with no endpoints is generated and an error
        message is given.
        """
        # Edge id
        self.id = self.make_id()

        # Points contained in edge (IDS)
        if len(vertex_pair) != 2 and len(vertex_pair) != 0:
            print "Your edge doesn't have the right number of points!"
            print "Generating edge with no endpoints."
            self.vertices = set([])
        else:
            self.vertices = set(vertex_pair)

    def __str__(self):
        "String conversion reveals id and the number of vertices"
        vertex_string = str([x for x in self.vertices])
        idstring = str(self.id)
        return "Edge object.\nID: {}\nVertices: {}".format(idstring,
                                                           vertex_string)

    def __len__(self):
        "Length reveals the number of vertices."
        return len(self.vertices)
#-------------------------------------------------------------------------


# The triangle class
#-------------------------------------------------------------------------
class triangle(geometry):
    """
    A class for keeping track of triangle edges
    
    Attributes: 
    -- last_used_id = the last used point ID. global/static. 
                      Shared by all class instances
    -- recycled_ids = A list of used IDs, which, after deletion, become 
                      available. Global/static. Shared by all class instances.
    -- points       = A list of points contained by the triangle
    -- id = The identifying number of the simplex. Useful. 
    -- area = The physical area of an equilateral triangle with the 
              giben edge length.
    -- angle = The angle of each angle of an equilateral triangle. 
               60 degrees.
    -- instances = The dictionary containing all instances of the class.

    Methods:
    -- increment_id() = A function that increments last_used_id
    -- reclaim_id(id) = Takes the id and adds it to recycled_ids 
                                 list.
    -- recycle_id() = A function that returns the least element of 
                               recycled_ids and deletes it from the list.
    -- make_id() = A function that generates a new id for a vertex.
                            First looks in recycled IDs. If there are non
                            available, uses increment ID.
    -- edges() = A list of edges contained by the triangle
    -- check_topology() = Checks to ensure the number of neighbors, 
                          edges, and vertices is okay.
    -- connect_to_triangle(other_triangle) = Sees if another triangle is 
                                             a neighbor. If it is, set each 
                                             triangle in the other 
                                             triangle's neighbors list.
    
    Example of call: t = triangle(point_list)
    where point_list = [p1,p2,p3]
    """

    # Static/global variables
    last_used_id = 0 # Last used triangle ID.
    recycled_ids = set([]) # List of unused vertex IDS <= last_used_id
    instances = {} # Contains all instances of the class

    # Geometric quantities
    area = np.sqrt(3) * edge.length**2 / 4 # surface area
    angle = np.pi/3 # Angle between edges

    # Functions
    def check_topology(self):
        """
        Ensure that the triangle has the correct numbers of
        vertices
        edges
        neighbors.
        """
        assert len(self.vertices) == 0 or len(self.vertices) == 3
        assert len(self.edges) == 0 or len(self.edges) == 3
        assert 0 <= len(self.neighbors) <= 3
        print "Topology is okay."

    def connect_to_triangle(self,other_triangle):
        """
        Looks at another triangle and checks to see if it shares an edge
        (or equivalently two points) with self.
        """
        # Do the triangles have an intersection?
        intersection = self.vertices & other_triangle.vertices
        # If the triangles share exactly 2 points, acknowledge their
        # neighborliness.
        if len(list(intersection)) == 2:
            self.neighbors.add(other_triangle.id)
            other_triangle.neighbors.add(self.id)
        self.check_topology()
        other_triangle.check_topology()
        
    # Duplicate checking
    @classmethod
    def find_duplicates(self,vertex_list=False,
                        edge_list=False,neighbor_list=False):
        """
        Searches through self.instances for triangles with the same
        vertices or the same edges or the same neighbors as the given lists.
        """
        duplicates = set([])
        for t in self.instances.values():
            if vertex_list:
                if t.vertices == set(vertex_list):
                    duplicates.add(t.id)
            if edge_list:
                if t.edges == set(edge_list):
                    duplicates.add(t.id)
            if neighbor_list:
                if t.neighbors == set(neighbor_list):
                    duplicates.add(t.id)
        return duplicates

    # Initialization function
    def __init__(self,point_list=[],edge_list=[],neighbor_list=[]):
        "Initialize a triangle"
        # Id
        self.id = self.make_id()
        # Vertexes of the triangle (IDS)
        if len(point_list) != 0 and len(point_list) != 3:
            print "Your triangle has the wrong number of points!"
            print "Creating an empty point list."
            self.vertices = set([])
        else:
            self.vertices = set(point_list) 
        # Edges of the triangle (IDS)
        if len(edge_list) != 0 and len(edge_list) != 3:
            print "Your triangle has the wrong number of edges!"
            print "Creating an empty list."
            self.edges = set([])   
        else:
            self.edges = set(edge_list)
        # The triangles this triangle is connected to.
        if len(neighbor_list) != 0 and len(neighbor_list) != 3:
            print "Your triangle has the wrong number of neighbors!"
            print "Creating an empty list."
            self.neighbors = set([])
        self.neighbors = set(neighbor_list) 

    def __str__(self):
        "String conversion reveals id, neighbors, and number of vertices"
        vertex_string = "Vertices: {}".format(str([x for x in self.vertices]))
        idstring = "ID: {}".format(str(self.id))
        edgestring = "Edges: {}".format(str(self.edges))
        neighbor_string = "Neighbors: {}".format(str(self.neighbors))
        outstring =  "Triangle Object.\n" + idstring + '\n' + \
            vertex_string + '\n' + edgestring + '\n' + neighbor_string
        return outstring

    def __len__(self):
        "Length reveals the number of vertices."
        return len(self.vertices)

    def __add__(self,other_triangle):
        "Equivalent to self.connect_to_triangle(other_triangle)."
        self.connect_to_triangle(other_triangle)
#-------------------------------------------------------------------------

