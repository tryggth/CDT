# make_spacetime_movie.sage.py
# Author: Jonah Miller (jonah.maxwell.miller@gmail.com)
# Date: July 30, 2012

# This is a script meant for Sage mathematics to generate a movie
# using *.mov2p1 files.

# Import statements
#-----------------------------------------------------------------------------
import numpy as np
import sys
#-----------------------------------------------------------------------------


# Functions
#-----------------------------------------------------------------------------
def load_data(filename):
    "Reads *.mov2p1 data from a file."
    with open(filename,'r') as f:
        data = f.read().rstrip().split('\n')
    data = [line.rstrip().split(' ') for line in data if len(line)]
    data = [[eval(simplex_count) for simplex_count in line] for line in data]
    return data
        
def make_movie(filename):
    moviedata = load_data(filename)
    frames = [list_plot(frame,plotjoined=True,
                        axes_labels=['Proper Time','Spatial Extent']) \
                  for frame in moviedata]
    anim = animate(frames)
    return anim

    
