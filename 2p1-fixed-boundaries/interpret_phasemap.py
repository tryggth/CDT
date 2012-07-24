#!/usr/bin/env python2

"""
interpret_phasemap.py
Author: Jonah Miller
"""

# Modules
import numpy as np
import scipy as sp
import matplotlib as mpl
import matplotlib.pyplot as plt
import sys


# Global constants
A_RANGE_MIN = 0.9 # Volume must be >=90% of target volume to be examined.
A_RANGE_MAX = 1.1 # Volume myust be <=110% of target volume to be examined.


# Functions
def set_parameters(cli_arguments):
    """
    Given a list of command-line arguments, (sys.argv[1:]), set and
    return the parameters we care about: target volume, acceptable
    variational range.
    """
    target_volume = sys.argv[1]
    file_names = sys.argv[2:]
    return [target_volume,acceptable_range,file_names]

def acceptable_range(target_volume):
    """
    Defines the acceptable range in spacetime tetrahedra based on
    target volume.
    """
    return target_volume * np.array([A_RANGE_MIN,A_RANGE_MAX])
