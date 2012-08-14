#!/usr/bin/env python2

"""
moves.py

Author: Jonah Miller (jonah.maxwell.miller@gmail.com)

This file is part of the sphere_generator program which generates
spheres of (as close as possible to) uniform curvature for a given
surface area by monte carlo simulation.

This file contains move functions for monte carlo moves.

Each move has 3 functions associated with it: 
     --- A complex generator which takes a vertex as input and calculates 
         all possible combinations of triangles around that vertex which 
         would be topologically acceptable for the move to operate on.
     --- A try function which takes a vertex as input and 
         calls the complex generator function and randomly selects 
         an acceptable complex. The try function then uses the methods 
         defined in state_tracking to calculate how the move will effect 
         global values such as the mean and standard deviation of the 
         sphere's curvature. This information can then be sent to the 
         Monte Carlo Metropolis algorithm to calculate whether or not to 
         accept the move. The try function returns the expected change in
         state and the complex that will be modified.
     --- An apply function actually applies a move to the spacetime. 
         It accepts a subcomplex as input. It has side effects.
"""


### Dependencies
#-------------------------------------------------------------------------
import numpy as np
import scipy as sp
import random
# Class data structures we need
import simplex_ancestors as sa
import simplex_descendants as sd
import state_manipulation as sm
import utilities as ut
import error_checking
import initialization
import state_tracking as st
#-------------------------------------------------------------------------


####---------------------------------------------------------------------####
#                               Classes                                     #
####---------------------------------------------------------------------####
# The complex class contains information about a set of simplices onto
# which it is topologically acceptable to apply a move to.
#---------------------------------------------------------------------------
class generalized_complex:
    """
    The complex class contains information about a set of simplices onto
    which it is topologically acceptable to apply a move to.

    Accepts triangle instances, lists of triangle ids, or lists of
    lists of triangle ids. This is the primary advantage of the
    complex class over just a list.

    The generalized_complex class is the parent class for complex
    child classes, which are for specialized purposes (for special
    moves, for instance). It should rarely if ever be called.
    """
    
    def __init__(self,triangle_list):
        triangle_id_list = []
        for t in triangle_list:
            if isinstance(t,sd.triangle):
                triangle_id_list.append(t.id)
            elif type(t) == list or type(t) == set or type(t) == tuple:
                triangle_id_list += list(t)
            elif type(t) == int:
                triangle_id_list.append(t)
            else:
                raise TypeError("Input to moves.complex must " + 
                                "be a collection, a triangle, or an ID.")
            self.triangles = set(triangle_id_list)
  
    def __str__(self):
        return str(self.triangles)

    def __len__(self):
        return len(self.triangles)

    def get_triangles(self):
        return self.triangles

class complex(generalized_complex):
    """
    The complex class contains information about a set of simplices onto
    which it is topologically acceptable to apply a move to.

    Accepts triangle instances, lists of triangle ids, or lists of
    lists of triangle ids. This is the primary advantage of the
    complex class over just a list.
    """
    
class complex42(generalized_complex):
    """
    The complex42 class is the same as the complex class, but includes
    a method to set the pair of vertices which have >=4 triangles
    attached to them and which we treat differently than the other 2
    boundary vertices. For use with the complex_4_to_2, try_4_to_2 and
    move_4_to_2 functions.
    """
    def set_vertices(self,vertices_with_4_triangles):
        """
        The 4->2 move requires a pair of vertices with 4
        triangles. These are they.
        """
        self.vertices_with_4_triangles = vertices_with_4_triangles

    def get_vertices(self):
        """
        The 4->2 move requires a pair of vertices with 4
        triangles. These are they.
        """
        return self.vertices_with_4_triangles
#---------------------------------------------------------------------------

   

####---------------------------------------------------------------------####
#                               Functions                                   #
####---------------------------------------------------------------------####

# Utility functions
#---------------------------------------------------------------------------
def extract_triangle(simplex_id_or_simplex):
    "Given a simplex id or simplex, returns the simplex."
    if type(simplex_id_or_simplex) == int:
        triangle = sd.triangle.instances[simplex_id_or_simplex]
    elif isinstance(simplex_id_or_simplex,sd.triangle):
        triangle = simplex_id_or_simplex
    else:
        raise ValueError("Move functions can only accept ints or triangles.")
    return triangle

def check_area_decreasing_validity(decrease):
    """
    If the area decrease of a move, decrease, would cause the sphere
    to cease to be topologically accpetable, return False. Otherwise,
    return True.  
    """
    min_acceptable_area = 4
    if sd.triangle.count_instances() - decrease < min_acceptable_area:
        return False
    else:
        return True

def set_neighbors_and_triangles(changed_vertices,changed_triangle_ids):
    """
    For each changed vertex in a move, finds and sets its triangles.
    For each triangle added in a move, connects it to its neihgbors.
    Finally, runs some error checking. Accepts collections only.
    """
    # Set triangles associated with a vertex
    for v in changed_vertices:
        v.find_and_set_triangles()
        
    # Set triangle neighbors
    for v in changed_vertices:
        v.connect_surrounding_triangles()

    # Finally run some error checking:
    for t_id in changed_triangle_ids:
        sd.triangle.instances[t_id].check_topology()
        sd.triangle.instances[t_id].check_edge_validity()
#---------------------------------------------------------------------------


# Functions for the 1->3 move
#---------------------------------------------------------------------------
"""
Move behavior: Add a vertex in the center of the triangle and connect it to 
               each of the other 3 vertices.
      0                 0
     / \               /|\
    /   \     -->     / | \
   /     \           /  0  \
  /       \         /  / \  \
 0---------0       0---------0

Volume increase: 2
"""
def complex_1_to_3(simplex_id_or_simplex):
    """
    Takes a simplex or simplex id as input and calculates what, if
    any, complices are topologically acceptable to operate on.

    The 1->3 move is by far the easiest of these, since it is
    trivially topologically acceptable for one and only one
    simplex. However, for completeness, it's included.
    """
    # Extract info
    triangle = extract_triangle(simplex_id_or_simplex)

    return complex([triangle.id])

def try_1_to_3(simplex_id_or_simplex):
    """
    Tries a 1->3 move and returns the move data that the metropolis
    algorithm will use to decide whether or not to accept a move.
    """
    # The complex
    cmpx = complex_1_to_3(simplex_id_or_simplex)

    # The triangle id. There should only be one.
    if len(cmpx.get_triangles()) == 1:
        triangle_id = list(cmpx.get_triangles())[0]
    else:
        raise ValueError("There should be only one "+
                         "triangle for the 1->3 complex.")

    # The triangle
    triangle = sd.triangle.instances[triangle_id]

    # The original vertices
    original_vertices = [sd.vertex.instances[p] for p in triangle.vertices]

    # We're adding 1 vertex, but each of 3 other vertices is connected
    # to 1 additional triangle each.

    # "replace" the 3 original vertices
    removed_vertices = [st.imaginary_vertex(len(v),False) \
                            for v in original_vertices]
    added_vertices = [st.imaginary_vertex(len(v)+1,True) \
                          for v in original_vertices]
    # Add the center vertex
    added_vertices.append(st.imaginary_vertex(3,True))

    return st.move_data(removed_vertices + added_vertices,cmpx,move_1_to_3)

def move_1_to_3(cmpx):
    "Applies a 1->3 move to the input complex, cmpx."
    # Extract the single triangle required for this move.
    if len(cmpx.get_triangles()) == 1:
        triangle_id = list(cmpx.get_triangles())[0]
    else:
        raise ValueError("There should be only one "+
                         "triangle for the 1->3 complex.")
    triangle = sd.triangle.instances[triangle_id]

    # Vertices
    original_vertices = triangle.get_vertices()
    
    # Edges 
    original_edges = triangle.get_edges()

    # Endpoints of edges
    endpoints = [e.get_vertex_ids() for e in original_edges]

    # Make the point that will trisect the triangle
    v = sd.vertex()
    sd.vertex.add(v)

    # Generate the vertices for the three new triangles to be created
    vertex_list = [points | set([v.get_id()]) for points in endpoints]
    
    # Make the new triangles
    new_triangles = [sm.build_triangle_and_edges(tri) for tri in vertex_list]

    # Delete the old triangle
    sm.remove_triangle(triangle)

    # Set neighbors, triangles for each vertex, and do some error checking
    vertex_ids = ut.set_union(vertex_list)
    vertices = [sd.vertex.instances[i] for i in vertex_ids]
    set_neighbors_and_triangles(vertices,new_triangles)

    return True
#---------------------------------------------------------------------------


# Functions for the 3->1 move
#---------------------------------------------------------------------------
"""
Move behavior: Undo the 1->3 move.

      0                 0
     /|\               / \
p   / | \     -->     /   \
   /  0  \           /     \
  /  / \  \         /       \
 0---------0       0---------0

Volume decrease: 2
"""
def complex_3_to_1(simplex_id_or_simplex):
    """
    Takes a simplex or simplex id as input and calculates what, if
    any, complices are topologically acceptable to operate on.

    The 3->1 move requires a vertex that is attached to 3
    simplices. There is a slight danger that volume decreasing moves
    can make a system topologically unacceptable. If this is the case,
    then the function returns False.

    If there is more than one acceptable complex, returns one at
    random.
    """
    volume_decrease = 2
    if not check_area_decreasing_validity(volume_decrease):
        return False

    # Extract triangle
    triangle = extract_triangle(simplex_id_or_simplex)

    # Extract vertices
    vertices = triangle.get_vertices()

    # If a vertex has only three triangles, consider it
    ids = lambda x: x.get_triangle_ids() # For convenience
    possibilities = [ids(v) for v in vertices if len(ids(v)) == 3]
    
    # Each boundary point on a possible complex must be attached to at
    # least 4 triangles
    # Function to extract vertex class objects from a triangle id
    obj = lambda i: sd.triangle.instances[i].get_vertices()
    # function to extract the boundary points of a complex possibility
    boundaries = lambda c: ut.set_union([set(obj(t)) for t in c]) \
        - ut.set_intersection([set(obj(t)) for t in c])
    # Function that tells us if every boundary vertex of a possibility
    # is connected to greater than 3 simplex
    acceptable = lambda p: len([b for b in boundaries(p) if len(b) > 3]) == 3

    # Only accept a possibility if each boundary vertex has 4 or more
    # triangles.
    complices = [complex(p) for p in possibilities if acceptable(p)]

    if len(complices) == 1:
        return complices[0]
    elif len(complices) > 1:
        return random.choice(complices)
    else:
        return False

def try_3_to_1(simplex_id_or_simplex):
    """
    Tries a 3->1 move and returns the move data that the metropolis
    algorithm will use to decide whether or not to accept the move. If
    the move is simply not topologically acceptable, returns false.
    """
    # The complex
    cmpx = complex_3_to_1(simplex_id_or_simplex)
    
    # If there are no topologically acceptable complices, stop right
    # now and return False. Otherwise, extract the triangles.
    if cmpx:
        triangles = cmpx.get_triangles()
    else:
        return False

    # Now that we have the triangle, extract the list of points of
    # each triangle.
    old_vertices = [sd.triangle.instances[t].get_vertices() \
                        for t in triangles]
    # The central vertex is the intersection of all the vertices in
    # the old triangles
    central_vertex = ut.set_intersection([set(t) for t in old_vertices])
    # The boundary vertices are the union of the vertices of the old
    # triangles minus the intersection:
    boundary_vertices = ut.set_union([set(t) for t in old_vertices]) \
        - central_vertex

    # We're removing the central vertex, but each boundary vertex will
    # be attached to 1 fewer triangle.
    
    # "Replace" the 3 original vertices
    removed_vertices = [st.imaginary_vertex(len(v),False) \
                            for v in boundary_vertices]
    added_vertices = [st.imaginary_vertex(len(v)-1,True) \
                          for v in boundary_vertices]
    # "Remove" the center vertex
    removed_vertices.append(st.imaginary_vertex(3,False))

    return st.move_data(removed_vertices + added_vertices,cmpx,move_3_to_1)

def move_3_to_1(cmpx):
    "Applies a 3->1 move to the input complex, cmpx."

    # Extract the complex, which contains three triangles.
    if len(cmpx.get_triangles()) == 3:
        original_triangles = [sd.triangle.instances[i] \
                                  for i in cmpx.get_triangles()]
    else:
        raise ValueError("There should be exactly 3 triangles for the " +
                         "3->1 complex.")

    # Extract the boundary vertices
    original_vertices = [set(t.get_vertices()) for t in original_triangles]
    central_vertex = list(ut.set_intersection(original_vertices))[0]
    boundary_vertices = ut.set_union(original_vertices) \
        - set([central_vertex])
    boundary_vertex_ids = [v.get_id() for v in boundary_vertices]

    # Extract the edges to be removed (there are 3)
    # Lambda functions for clarity
    get_edges = lambda i: set(sd.triangle.instances[i].get_edges())
    intersected_element = lambda L: ut.only_element(ut.set_intersection(L))
    # Nested list comprehensions: take the intersection of the set of
    # edges associated with pairs of the triangles we care
    # about. These edges are flagged for deleteion.
    shared_edges = [intersected_element([get_edges(i) for i in c]) \
                        for c in ut.k_combinations(cmpx.get_triangles())]
    # There should only be 3 shared edges. If there are more, we messed up.
    assert len(shared_edges) == 3
    
    # Make the new triangle
    new_triangle = sm.build_triangle_and_edges(boundary_vertex_ids)

    # Clean up
    for t in original_triangles: # Delete the old triangles:
        sm.remove_triangle(t)
    for e in shared_edges: # Delete the old edges.
        sm.remove_edge(e)
    sd.vertex.delete(central_vertex) # Delete the old vertex

    # Set triangles, neihgbors, and check topology
    set_neighbors_and_triangles(boundary_vertices,[new_triangle])

    return True
#---------------------------------------------------------------------------


# Functions for the 2->4 move
#---------------------------------------------------------------------------
"""
Move behavior: Bisect the edge the triangles share and add a vertex there.
               Connect the new vertex to each of the other 4 vertices.
      0                   0
     / \                 /|\
    /   \               / | \
   /     \             /  |  \
  /       \           /   |   \
 0---------0   -->   0----0----0 
  \       /           \   |   /
   \     /             \  |  /
    \   /               \ | /
     \ /                 \|/
      0                   0 

Volume increase: 2
"""
def complex_2_to_4(simplex_id_or_simplex):
    """
    Takes a simplex or simplex id as input and calculates what, if
    any, complices are topologically acceptable to operate on.

    The 2->4 move requires only that two triangles are neighbors. Thus
    its is almost always (if not always) topologically acceptable to
    use. If there is more than on acceptable complex (there will
    almost always be 3), return one at random.
    """
    # Extract data for the first triangle
    # Extract triangle
    t1 = extract_triangle(simplex_id_or_simplex)
    # Extract neighbors
    t1_neighbors = t1.get_neighbors()
    # t1 and one of its neighbors will always form the complex we need.
    complices = [complex([t1,neighbor]) for neighbor in t1_neighbors]
    # Since this is a volume increasing move, there are no other
    # topological concerns, so we just pick a complex at random.
    return random.choice(complices)

def extract_from_complex_2_to_4(cmpx):
    """
    A sub-function for try_2_to_4. Extracts vertex information from
    the complex passed to it.
    """
    # Extract information from the complex:
    triangles = [sd.triangle.instances[i] for i in cmpx.get_triangles()]
    # Extract vertices for each triangle
    old_vertices = [set(triangle.get_vertices()) for triangle in triangles]
    # Extract edges for each triangle
    old_edges = [set(triangle.get_edges()) for triangle in triangles]

    # Find the edge shared by the two triangles. There should be only
    # one of these, but we should be careful and check anyways.
    shared_edges = ut.set_intersection(old_edges) 
    if len(shared_edges) == 1:
        shared_edge = list(shared_edges)[0]
    else:
        raise ValueError("The triangles should only share one edge!")

    # The shared vertices are just the endpoints of the edge
    shared_vertices = set(shared_edge.get_vertices())

    # We can now extract the vertices not on the shared edge using set
    # difference.
    unshared_vertices = ut.set_union(old_vertices) - shared_vertices
    # There should be only 2 unshared vertices.
    assert len(unshared_vertices) == 2

    # We have everything we need now, so return it:
    return [unshared_vertices,shared_edge,triangles]

def try_2_to_4(simplex_id_or_simplex):
    """
    Tries a 2->4 move and returns the move data that the metropolis
    algorithm will use to determine whether or not to accept the
    move. If the move is simply not topologically acceptable, returns
    false.
    """
    # Increases/decreases per imaginary vertex
    unshared_increase = 1
    central_increase = 4

    # The complex
    cmpx = complex_2_to_4(simplex_id_or_simplex)
    # If there are no topologically acceptable complices, stop right
    # now and return False. Otherwise, extract the triangles.
    if cmpx:
        unshared_vertices = extract_from_complex_2_to_4(cmpx)[0]
    else:
        return False
    
    # The shared vertices remain unchanged and we add a vertex with 4
    # triangles in the center. The unshared vertices each get one
    # additional triangle attached to them.
    
    # First, "replaces" the original unshared vertices
    removed_vertices = [st.imaginary_vertex(len(v),False) \
                            for v in unshared_vertices]
    added_vertices = [st.imaginary_vertex(len(v)+unshared_increase,True) \
                          for v in unshared_vertices]

    # Add the center vertex
    added_vertices.append(st.imaginary_vertex(central_increase,True))

    return st.move_data(removed_vertices + added_vertices,cmpx,move_2_to_4)

def extract_from_complex_move_2_to_4(cmpx):
    "Extracts useful triangle,edge,and vertex information for move_2_to_4."
    # Extract the triangles required for this move.
    # Error checking
    triangles = set([sd.triangle.instances[i] for i in triangle_ids])
    # Get edges
    edges = [set(t.get_edges()) for t in triangles]
    # Get vertices
    vertices = [set(t.get_vertices()) for t in triangles]
    # Find the shared ed

def move_2_to_4(cmpx):
    "Applies a 2->4 move to the input complex, cmpx."
    
    # We need to extract vertex and edge information. Check to make
    # sure the complex is correct too.
    if len(cmpx.get_triangles()) == 2:
        extracted_complex = extract_from_complex_2_to_4(cmpx)
    else:
        raise ValueError("There should be 2 triangles for the 2->4 complex.")
    unshared_vertices = extracted_complex[0] # Vertices opposite the shared edge
    shared_edge = extracted_complex[1] # The edge shared between the
                                       # two triangles
    shared_vertices = set(shared_edge.get_vertices())
    triangles = extracted_complex[2] # The two triangle objects used.

    # With all this information at our disposal, we're ready to start
    # changing things:
    
    # Bisect the edge, get a set of 2 new edges and the new
    # vertex. The old shared edge now no longer exists.
    new_edges,new_vertex = sm.bisect_edge(shared_edge)
    
    # We need to make 4 triangles: 2 connected to each unshared
    # vertex, one per new edge.
    new_triangles = set([])
    for vertex in unshared_vertices:
        vertex_id = vertex.get_id()
        for edge in new_edges:
            endpoints = set(edge.get_vertex_ids())
            assert len(endpoints) == 2 # If this is not true,
                                       # something went wrong.
            triangle_vertices = endpoints | set([vertex_id])
            new_triangles.add(sm.build_triangle_and_edges(triangle_vertices))

    # With the new triangles constructed, the old triangles must be destroyed.
    for t in triangles:
        sm.remove_triangle(t)

    # Set the triangles for each vertex that was changed:
    changed_vertices = ut.set_union([shared_vertices,
                                     unshared_vertices,
                                     set([new_vertex])])
    set_neighbors_and_triangles(changed_vertices,new_triangles)
    
    return True
# ---------------------------------------------------------------------------


# Functions for the 4->2 move
#---------------------------------------------------------------------------
"""
Move behavior: Undo a 4->2 move. Take the central vertex of a complex of 4
               triangles and remove it and one edge to form 2 triangles.
       0                 0
      /|\               / \
     / | \             /   \
    /  |  \           /     \
   /   |   \         /       \
  0----0----0  -->  0---------0
   \   |   /         \       /
    \  |  /           \     /
     \ | /             \   /
      \|/               \ /
       0                 0
 Volume decrease: 2
"""
def complex_4_to_2(simplex_id_or_simplex):
    """
    Takes a simplex or simplex id as input and calculates what, if
    any, complices are topologically acceptable to operate on.

    The 4->2 move requires that 4 triangles are connected to a single
    vertex. If there is more than one acceptable complex, then the
    function returns one at random.
    """

    #      0
    #     /|\       
    #    / | \       For a complex to be acceptable, at least one pair of
    #   /  |  \      vertices opposite each other in the complex (marked 1
    #  1---0---1     on the diagram to the left), must have at least 4 
    #   \  |  /      triangles attached to each vertex.
    #    \ | /       
    #     \|/       
    #      0        

    # The accepted minimum number of triangles attached to a vertex in
    # the chosen opposite pair of vertices:
    acceptable_minimum = 4

    # Local syntactic sugar:
    vertex_pairs = lambda c: st.find_opposite_vertices_in_complex(c)
    # Number of triangles attached to a given vertex
    num_ts = lambda vertex: len(vertex.get_triangles())
    # A vertes pair is okay if each vertex in it is attached to 4 triangles:
    pair_ok = lambda pair: bool([v for v in pair \
                                     if num_ts(v) >= acceptable_minimum])

    # Ensure that the volume decrease is accepable for the sphere at all.
    volume_decrease = 2
    if not check_area_decreasing_validity(volume_decrease):
        return False

    # Extract data for the first triangle
    t1 = extract_triangle(simplex_id_or_simplex) # Extract triangle
    t1_vertices = set(t1.get_vertices()) # Vertices surrounding the
                                         # first triangle.
    
    # Vertices to center the complex have 4 triangles associated with them
    possible_vertices = [v for v in t1_vertices if len(v) == 4]
    # Thus possible complices are built around the possible vertices
    possibilities = [set(v.get_triangles()) for v in possible_vertices]
    
    # For a possible complex to be acceptable, it must have at least
    # one pair of vertices opposite each other for which each vertex
    # is attached to at least 4 triangles. 
    possibilities = [(cmpx,pair) for cmpx in possibilities \
                         for pair in vertex_pairs(cmpx)
                         if pair_ok(pair)]

    # If there are any possibilities left, return one at
    # random. Otherwise, return false.
    if possibilities:
        chosen_complex = random.choice(possibilities)
        output_complex = complex42(chosen_complex[0])
        output_complex.set_vertices(pair)
        return output_complex
    else: 
        return False
# ---------------------------------------------------------------------------
