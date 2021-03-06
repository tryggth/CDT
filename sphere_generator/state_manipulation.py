"""
state_manipulation.py

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
# Class data structures we need
from simplex_ancestors import *   
from simplex_descendants import *
import utilities as ut
import error_checking
#-------------------------------------------------------------------------


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
def bisect_edge_using_id(edge_id):
    """
    Uses the edge.bisect() function to bisect an edge and update the
    hash tables properly. Takes an edge id.
    """
    # Bisect the edge and generate two pairs of points, each
    # representing one of the two new edges.
    new_edges,newpoint = edge.instances[edge_id].bisect() 

    # Delete the old edge
    edge.delete(edge_id)

    # Add the new point to the hash table
    new_vertex = vertex([],newpoint)
    vertex.add(new_vertex)
    
    # Add the new edges to the edges hash table. We want to return the
    # new idge IDs (because why not), so a temporary list is generated
    # for this purpose.
    new_edge_instances = set([])
    for endpoint_pair in new_edges:
        e = edge(endpoint_pair)
        edge.add(e)
        new_edge_instances.add(e)
    return [new_edge_instances,new_vertex]

def bisect_edge(edge_id_or_instance):
    """
    Uses the edge.bisect() function to bisect an edge and update the
    hash tables properly. Takes an edge id or instance.
    """
    if type(edge_id_or_instance) == int \
            and edge_id_or_instance in edge.instances.keys():
        edge_id = edge_id_or_instance
    elif isinstance(edge_id_or_instance,edge):
        edge_id = edge_id_or_instance.get_id()
    else:
        raise TypeError("bisect_edge accepts only edge ids "+
                        "or edge instances.")
    return bisect_edge_using_id(edge_id)
# -------------------------------------------------------------------------



### Functions that make simplices and sub-simplices together. The
### workhorses of state manipulation.
### -------------------------------------------------------------------------
def build_sub_edges_of_triangle(point_list_or_triangle_object):
    """
    Takes either a triangle class instance or a set of three
    points as input and makes the edges that a triangle defined this
    way would contain. Checks for redundancies.
    """
    if type(point_list_or_triangle_object) == list:
        vertices = set(point_list_or_triangle_object)
    elif type(point_list_or_triangle_object) == tuple:
        vertices = set(point_list_or_triangle_object)
    elif type(point_list_or_triangle_object) == set:
        vertices = point_list_or_triangle_object
    elif isinstance(point_list_or_triangle_object,triangle):
        vertices = set(point_list_or_triangle_object.vertices)
    else: # This is a way to give more verbose error messages
        print "ERROR: You have passed the function"
        print "'build_sub_edges_of_triangle'"
        print "an object that is not a list, a set, a tuple," + \
            "or a triangle class instance."
        print "The object was: {}".format(point_list_or_triangle_object)
        assert isinstance(point_list_or_triangle_object,triangle) \
            or type(point_list_or_triangle_object) == list \
            or type(point_list_or_triangle_object) == tuple \
            or type(point_list_or_triangle_object) == set
    # A triangle should only have 3 vertices
    assert error_checking.check_length(vertices,3,'triangle',
                                       'build_sub_edges_of_triangle')

    # Now calculate the points the subsimplices contain
    ordered_pairs = ut.k_combinations(vertices,2)
    
    # The ids of our subsimplices. This will be a return value.
    edge_ids = set([]) 
    
    # For each edge, check to see if it exists. If it does not, make it.
    for pair in ordered_pairs:
        duplicates = edge.find_duplicates(pair)
        # If there are no duplicates, then we need to make this edge
        if len(duplicates) == 0: 
            # And tell our triangle it exists
            new_edge = edge(pair)
            edge_ids.add(new_edge.id)
            # And add it to the edge instances hash table
            edge.add(new_edge)
        # If there is exactly one duplicate, then we simply add its id
        # to the triangle's edge ID list
        elif len(duplicates) == 1:
            edge_ids.add(duplicates[0])
        else:  # If there is more than one duplicate, something went
               # very wrong.
            error_checking.too_many_duplicates('edge',pair,duplicates,1,
                                               'build_sub_edges_of_triangle')
    
    # There should be exactly 3 edge ids.
    assert error_checking.check_length(edge_ids,3,'edge_ids',
                                       'build_sub_edges_of_triangle')
    return edge_ids


def build_triangle_and_edges(point_list):
    """
    Builds a triangle by creating the triangle object (if it doesn't
    already exist) and all subsimplices---just edges---if they don't
    already exist. Checks for redundancies. If the triangle already
    exists, returns its id.
    """
    # local constants
    triangle_length = 3 # The total number of vertices in a triangle.

    # Make sure we actually have a collection of 3 triangle vertices.
    assert len(point_list) == triangle_length

    # Generate a list of point ids.
    point_list = [vertex.parse_input(t).get_id() for t in point_list]
    
    # Typecast to eliminate duplicates
    vertices = set(point_list)

    # Build the edges contained by the triangle and return their ids.
    # build_sub_edges_of_triangle runs some error checking too.
    edge_ids = build_sub_edges_of_triangle(vertices)

    # Check to see if the triangle we want to build already exists
    duplicates = triangle.find_duplicates(vertices,edge_ids)

    # If there is one duplicate, work with it. If there is more than
    # one, raise an error. If there are no duplicates, make that
    # triangle!
    if len(duplicates) == 0:
        t = triangle(vertices,edge_ids) # Make the triangle
        triangle_id = t.id # Ensure we have the id for the triangle
        triangle.add(t) # Add the triangle to its hash table
    elif len(duplicates) == 1:
        triangle_id = duplicates[0] # The id for our triangle
    else: # If the length of duplicates is not 0 or 1, something went
          # very wrong.
        assert error_checking.too_many_duplicates('triangle',vertices,
                                                  duplicates,1,
                                           'build_triangle_and_edges')
        
    return triangle_id

def remove_triangle(triangle_id_or_instance):
    """
    Delets a triangle from the corresponding hash table. Moreover,
    removes it from its' neighbors list of neighbors.
    """
    # Interpret input
    if type(triangle_id_or_instance) == int:
        triangle_id = triangle_id_or_instance
        triangle_instance = triangle.instances[triangle_id]
    elif isinstance(triangle_id_or_instance,triangle):
        triangle_instance = triangle_id_or_instance
        triangle_id = triangle_instance.get_id()
    else:
        raise ValueError("We need a triangle instance or an ID here!")

    # Delete the triangle's id from neighbors
    for n in triangle_instance.get_neighbors():
        n.remove_neighbor(triangle_id)

    # Delete the triangle itself
    triangle.delete(triangle_id)

def remove_edge(edge_id_or_instance,triangle_list=False):
    """
    Deletes an edge from the corresponding hash table. If a triangle
    list was given, removes the id from each triangle's edge id
    list.
    """

    # Interpret input
    if type(edge_id_or_instance) == int:
        edge_id = edge_id_or_instance
        edge_instance = edge.instances[edge_id]
    elif isinstance(edge_id_or_instance,edge):
        edge_instance = edge_id_or_instance
        edge_id = edge_instance.get_id()
    else:
        raise ValueError("We need an edge instance or an ID here!")

    # If triangle list given, delete the edge id from triangle.edges
    if triangle_list:
        if type(triangle_list) == int:
            triangle.instances[triangle_list].edges.remove(edge_id)
        elif type(triangle_list) == list:
            if isinstance(triangle_list[0],triangle):
                for t in triangle_list:
                    t.edges.remove(edge_id)
            else:
                for i in triangle_list:
                    triangle.instances[i].edges.remove(edge_id)
        elif type(triangle_list) == tuple or type(triangle_list) == set:
            if isinstance(list(triangle_list)[0],triangle):
                for t in triangle_list:
                    t.edges.remove(edge_id)
            else:
                for i in triangle_list:
                    triangle.instances[i].edges.remove(edge_id)
        else:
            raise TypeError("You did not pass a list of triangle ids.")
        
    # In any case, delete the edge.
    edge.delete(edge_instance)
### -------------------------------------------------------------------------

