#!/usr/bin/env python2
"""
k3_of_k0.py
Author: Jonah Miller

This program takes *.critical_surface_fit data and constructs k3 as a
function of k0 using it.

call example:
./k3_of_k0.py phasemaps/T028-V01024-BItetra.txt-BFtetra.txt-k0.5-4.5-kkk0.6-1.6.phasemap.critical_surface_fit 0.5 1.0 1.5 2.0 2.5
"""

import sys

def main(k0,filename):
    # Read in the data
    with open(filename,'r') as f:
        data = f.read().split('\n')
    # If there are any empty lines, remove them
    data = [line for line in data if len(line) > 0]
    # Prune the data.
    data = [line for line in data if line[0] != '#']
    # The first line is now an equation for k3 as a function of k0
    # with some coupling constants. The second line is the names of
    # the coupling constants. The third is those coupling constants.
    eqn = data[0]
    constant_names = data[1].split(' ') # A list of coupling
                                        # constant names. They're strings
                                        # right now.
    constants = data[2].split(' ')      # A list of coupling
                                        # constants. They're strings
                                        # right now.
    # Put the constants into the equation
    for i in range(len(constant_names)):
        eqn = eqn.replace(constant_names[i],constants[i])
    eqn = eqn.replace('k0',k0) # Put in k0.
    k3 = eval(eqn)
    print k3
    return k3

if __name__ == "__main__":
    filename = sys.argv[1]
    for k0 in sys.argv[2:]:
        main(k0,filename)
