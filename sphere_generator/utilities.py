#!/usr/bin/env python2

"""
utilities.py

Author: Jonah Miller (jonah.maxwell.miller@gmail.com)

This file is part of the sphere_generator program which generates
spheres of (as close as possible to) uniform curvature for a given
surface area by monte carlo simulation.

This module contains functions not necessarily linked to the sphere
generation program that are simply useful little tricks.
"""


# Dependancies
#-------------------------------------------------------------------------
import itertools # For combinations
import math # For factorials
#-------------------------------------------------------------------------


# Functions

# Combinatorics stuff
#-------------------------------------------------------------------------
def swap(a,b):
    "sets a = b, and b = a."
    temp = a
    a = b
    b = temp
    del temp
    return [a,b]

def binomial_coefficient(n,k):
    """
    Returns n choose k.
    """
    if k > n: # Make sure n is greater than k.
        print "n < k. Reversing your inputs."
        n,k = swap(n,k)
    return math.factorial(n)/(math.factorial(n-k)*math.factorial(k))

def k_combinations(big_set,k=2):
    """
    Finds and returns the k-combinations of big_set---i.e., the
    subsets of big_set with cardinality k. For convenience reasons,
    the default value of k is 2.
    """
    S = set(big_set) # In case big_set was given to us as a
                     # list. Eliminates redundant elements.
    if k > len(S): # If k > len(S), combinations don't make sense.
        print "k > len(S). Setting k = len(S)."
        k = len(S)

    combinations = [] # The k_combinations to return
    count = 0 # Test to make sure that the number of k_combinations is
              # equal to len(S) choose k.
    for c in itertools.combinations(S,k): # Get the combinations we need.
        count += 1
        combinations.append(set(c))

    # Just to prevent unexpected errors
    assert count == binomial_coefficient(len(S),k) 

    return combinations
#-------------------------------------------------------------------------


