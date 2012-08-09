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
# Class data structures we need
import simplex_ancestors as sa
import simplex_descendants as sd
import state_manipulation as sm
import utilities as ut
import error_checking
import state_tracking as st
#-------------------------------------------------------------------------


# 
