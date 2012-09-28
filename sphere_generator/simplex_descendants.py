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
import utilities as ut
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

    @classmethod
    def is_redundant_id(self,vertex_id):
        """
        Checks to see if a vertex_id is currently in use or not. If
        the vertex is currently in use, return True. If the vertex is
        currently in recycled IDs, return the string
        'recycled'. Otherwise, returns False.
        """
        v = vertex_id # For compactness
        if v <= self.last_used_id and v not in self.recycled_ids:
            return True
        elif v <= self.last_used_id and v in self.recycled_ids:
            return 'recycled'
        else:
            return False
            
    @classmethod
    def isinstance(self,object_instance):
        "Tests if an object_instnace is an instance of the vertex class."
        return isinstance(object_instance,vertex)

    # Initialization function
    def __init__(self, triangle_list = [],new_id=0,edge_list=[]):
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
        # Local edge list
        self.edges = set(edge_list)

    def __str__(self):
        "String conversion reveals id and the number of triangles."
        t_string = str([x for x in self.triangles])
        e_string = str([x for x in self.edges])
        idstring = str(self.id)
        outstring = "Vertex object.\nID: {}\nTriangles: {}\nEdges: {}"
        outstring = outstring.format(idstring, t_string, e_string)
        return outstring

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

    def get_edges(self):
        "Returns the objects, not the ids."
        return [edge.instances[e] for e in self.edges]

    def get_edge_ids(self):
        "Returns the ids, not the objects."
        return self.edges

    def set_edges(self,edge_list):
        "Sets self.edges"
        self.edges = set(edge_list)

    def find_and_set_edges(self):
        """
        Finds attached edges and sets self.edges to include them.
        """
        new_edges = self.find_edges()
        self.set_edges(new_edges)
        return new_edges

    def get_triangles(self):
        "Returns the objects, not the ids."
        return [triangle.instances[t] for t in self.triangles]

    def get_triangle_ids(self):
        "Returns the ids, not the objects."
        return self.triangles

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

    def triangles_shared_with(self,other_instance_or_id):
        """
        Tests to see whether self is a vertex of one of the same
        triangles as other. Other can be a vertex instance or a vertex
        ID. If triangles exist, returns them. Otherwise, returns false.
        """
        # Parse input 
        other = self.parse_input(other_instance_or_id)
        # The triangles shared with self:
        shared_triangles = ut.set_intersection([set(self.get_triangles()),
                                                set(other.get_triangles())])
        return shared_triangles
        
    def shares_a_triangle_with(self,other_instance_or_id):
        """
        Tests to see whether self is a vertex of one of the same
        triangles as other. Returns a boolean value.
        """
        return bool(self.triangles_shared_with(other_instance_or_id))
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
    @classmethod
    def isinstance(self,object_instance):
        "Tests if an object_instnace is an instance of the edge class."
        return isinstance(object_instance,edge)

    @classmethod
    def exists(self,endpoint_pair):
        """
        Tests to see if the edge defined by a vertex pair exists. If
        it does, return the edge ID.
        """
        if len(endpoint_pair) != 2:
            raise ValueError("Must be a pair of vertices or vertex IDs!")
        endpoints = set([vertex.parse_input(i) for i in endpoint_pair])
        possible_edges = list(endpoints)[0].get_edges()
        same_edges = [e.get_id() for e in possible_edges \
                          if endpoints == set(e.get_vertices())]
        if same_edges:
            return same_edges
        else:
            return False
        

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

    def check_topology(self,return_value=False):
        """
        Check to ensure that the number of endpoints is correct.
        """
        assert len(self.vertices) == 2 or len(self.vertices) == 0
        if return_value:
            print "Topology is okay."

    def set_vertices(self,new_vertices):
        "Takes an input list or set and resets self's vertices/endpoints."
        if len(new_vertices) == 2 or len(new_vertices) == 0:
            self.vertices = set(new_vertices)
            return self.vertices
        else:
            print "Error! Wrong number of vertices! Nothing changed."
            return self.vertices

    def get_vertices(self):
        "Returns the instances, not the ids."
        return [vertex.instances[v] for v in self.vertices]

    def get_vertex_ids(self):
        "Returns the ids."
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
    @classmethod
    def isinstance(self,object_instance):
        "Tests if an object_instnace is an instance of the vertex class."
        return isinstance(object_instance,triangle)

    def check_neighbor_edge_correlation(self):
        """
        Ensures there is a bijection between neighbors and edges.
        """
        edges = set([])
        for n in self.get_neighbors():
            shared_edges = ut.set_intersection([n.edges,self.edges])
            if len(shared_edges) == 1:
                edges.add(list(shared_edges)[0])
        if edges != self.edges:
            print "\nErroneous Triangle!\n"
            print self
            raise ValueError("Each neighbor does not share exactly 1 edge!")
                        

    def check_topology(self,return_value=False):
        """
        Ensure that the triangle has the correct numbers of
        vertices
        edges
        neighbors.
        """
        assert len(self.vertices) == 0 or len(self.vertices) == 3
        assert len(self.edges) == 0 or len(self.edges) == 3
        if not 0 <= len(self.neighbors) <= 3:
            print "\nErroneous Triangle!\n"
            print self
            raise ValueError("Number of neighbors is wrong!")
        if return_value:
            print "Topology is okay."

    def check_topology_v2(self,return_value=False):
        """
        Like check_topology, but only allows 3 for each value. Useful
        after initialization.
        """
        if len(self.vertices) != 3:
            print self
            raise ValueError("Number of vertices is wrong!")
        if len(self.edges) != 3:
            print self
            raise ValueError("Number edges is wrong!")
        if len(self.neighbors) != 3:
            print self
            raise ValueError("Number vertices is wrong!")
        if return_value:
            print "Topology is okay."

    def check_edge_validity(self):
        """
        Checks to make sure that the number of endpoints of edges is
        the same as the number of vertices.
        """
        edges = self.get_edges()
        endpoints = [e.get_vertex_ids() for e in edges]
        endpoint_union = ut.set_union([e.get_vertex_ids() for e in edges])
        assert endpoint_union == self.get_vertex_ids()

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
        
    def remove_neighbor(self,other_triangle):
       """
       other_triangle is in self.neighbors, remove it from the
       neighbor list. Accepts ids only.
       """
       if other_triangle in self.neighbors:
           self.neighbors.remove(other_triangle)

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

    def get_vertices(self):
        "Returns the instances, not the ids."
        return [vertex.instances[v] for v in self.vertices]

    def get_vertex_ids(self):
        "Returns the ids."
        return self.vertices

    def get_edges(self):
        "Returns the instances, not the ids."
        return [edge.instances[e] for e in self.edges]

    def get_edge_ids(self):
        "Returns the ids."
        return self.edges

    def get_neighbors(self):
        "Returns the instances, not the ids."
        return [triangle.instances[t] for t in self.neighbors]

    def get_neighbor_ids(self):
        "Returns the ids."
        return self.neighbors

    def make_outstring(self):
        """
        Makes a string of the form '(v1, v2, v3)' for vertices. Useful
        for output.
        """
        assert len(self.vertices) == 3
        v = list(self.vertices)
        return "({} {} {})".format(v[0],v[1],v[2])

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

