"""
simplex_ancestors.py

Author: Jonah Miller (jonah.maxwell.miller@gmail.com)

This file is part of the sphere_generator program which generates
spheres of (as close as possible to) uniform curvature for a given
surface area by monte carlo simulation.

This module contains the the geometry parent class. geometry contains
functions that manipulate the dictionaries that contain triangles,
edges, and vertices.
"""



### Dependencies
#-------------------------------------------------------------------------
import numpy as np
import scipy as sp
#-------------------------------------------------------------------------



####---------------------------------------------------------------------####
#                               Classes                                     #
####---------------------------------------------------------------------####
# The Geometric object base class
#---------------------------------------------------------------------------#
class geometry:
    """
    The ancestor class for the geometric objects in
    simplex_descendants.py. Contains the methods that act on all of
    them to control id recycling and construction.

    Do not initialize this class. It is not useful to do so.
    """
    
    @classmethod
    def increment_id(self):
        """
        Increment the static variable last_used_id by 1. Then return the
        new id for use.
        """
        self.last_used_id += 1
        return self.last_used_id
    
    @classmethod
    def reclaim_id(self,object_id):
        """
        Take an input ID and add it to the list of recycled IDs.
        Return the list of recycled ids. Because why not?
        """
        self.recycled_ids.add(object_id)
        return self.recycled_ids

    @classmethod
    def recycle_id(self):
        """
        Looks for the minimum element of your recycled IDs. Returns it and
        removes it from the set.
        """
        # Find the minimum element of recycled Ids.
        minval = min(self.recycled_ids)
        # Remove the minimum element from the set
        self.recycled_ids.remove(minval)
        # Return the returned element for reuse
        return minval
    
    @classmethod
    def make_id(self):
        """
        Generate an ID for use. First look in recycled IDs. If none
        are available, make a brand new ID.
        """
        if len(self.recycled_ids) > 0:
            new_id = self.recycle_id()
        else:
            new_id = self.increment_id()
        return new_id
    
    @classmethod
    def add(self,object_instance):
        """
        Adds an object of the same type as self to the dictionary
        containing those objects.
        """
        self.instances[object_instance.id] = object_instance

    @classmethod
    def delete(self,argument):
        """
        Deletes an object of the same type as self from the dictionary
        containing those objects and reclaims its id.

        Accepts an id, a list of ids, an object, or a list of objects.
        """
        if type(argument) == list: # i.e., if we were passed a list
            if len(argument) == 0: # If list is empty, do nothing
                pass
            if type(argument[0]) == int: # i.e., if we're working with ids
                for object_id in argument:
                    del self.instances[object_id] # remove the id from
                                                  # the hash table
                    self.reclaim_id(object_id) # Reclaim the id
                pass
            else: # I.e., if we're working with objects
                for object_instance in argument:
                    # Remove the id from the hash table
                    del self.instances[object_instance.id]
                    self.reclaim_id(object_instance.id) # Reclaim the id
                pass
        elif type(argument) == int: # I.e., we were passed a single object id
            del self.instances[argument]
            self.reclaim_id(argument)
            pass
        else: # i.e., we were passed a single object instance
            del self.instances[argument.id]
            self.reclaim_id(argument.id)
            pass
        
    @classmethod
    def delete_all(self):
        """
        Deletes all instances of the same type as self.
        """
        to_delete = list(self.instances.keys())
        if len(to_delete) > 0: # Only delete stuff if there's stuff to
                               # delete.
            self.delete(to_delete)

    @classmethod
    def list_ids(self):
        "Returns a list of all object ids of type self."
        return list(self.instances.keys())

    @classmethod
    def count_instances(self):
        "Counts the number of instances of object of type self."
        return len(list(self.instances))

    def get_id(self):
        "Gets the object id."
        return self.id

    @classmethod
    def parse_input(self,id_or_instance):
        """
        Enables a function to take an id or an instance as
        input. Always returns a class instance rather than an
        id. Takes in either an id or an instance as input and returns
        the corresponding class instance.
        """
        if type(id_or_instance) == int \
                and id_or_instance in self.instances.keys():
            other = self.instances[id_or_instance]
        elif self.isinstance(id_or_instance):
            other = id_or_instance
        else:
            raise TypeError("I need a geometry subclass instance or ID.")
        return other
#-------------------------------------------------------------------------

