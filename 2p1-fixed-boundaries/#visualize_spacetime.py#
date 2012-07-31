#!/usr/bin/env python2

# visualize_spacetime.py
# Author: Jonah Miller (jonah.miller@colorado.edu)
# Date: June 29, 2012

# This program loads a spacetime or a spacetime and looks at
# observables such as volume as a function of proper time.

# TODO: Integrate looking at an entire ensemble


# Modules
#-------------------------------------------------------------------------
import sys # For yes/no
import numpy as np # For array control
import scipy as sp
import matplotlib as mpl # For data visualization
import matplotlib.pyplot as plt
import matplotlib.animation as animation
#-------------------------------------------------------------------------

# Yes/no function. From a recipe:
# http://code.activestate.com/recipes/577058/
def query_yes_no(question, default="yes"):
    """Ask a yes/no question via raw_input() and return their answer.

    "question" is a string that is presented to the user.
    "default" is the presumed answer if the user just hits <Enter>.
        It must be "yes" (the default), "no" or None (meaning
        an answer is required of the user).

    The "answer" return value is one of "yes" or "no".
    """
    valid = {"yes":True,   "y":True,  "ye":True,
             "no":False,     "n":False}
    if default == None:
        prompt = " [y/n] "
    elif default == "yes":
        prompt = " [Y/n] "
    elif default == "no":
        prompt = " [y/N] "
    else:
        raise ValueError("invalid default answer: '%s'" % default)
    
    while True:
        sys.stdout.write(question + prompt)
        choice = raw_input().lower()
        if default is not None and choice == '':
            return valid[default]
        elif choice in valid:
            return valid[choice]
        else:
            sys.stdout.write("Please respond with 'yes' or 'no' "\
                             "(or 'y' or 'n').\n")


# Reads in a *.3sx2p1 file and produces a vector of two lists. The
# first list contains information from the first line of the data
# file. It tells us:

# [BCTYPE, STOPOLOGY, NUM-T, N-INIT, *LAST-USED-POINT*, *LAST-USED-3SXID*, 
#  N0, N1-SL, N1-TL, N2-SL, N2-TL, N3-TL-31, N3-TL-22,
#  N1-SL-BOUNDARY, N3-31-BOUNDARY, N3-22-BOUNDARY,
#  *eps* *k0* *k3* *alpha*]

# The second element of the vector is a list of lists of simplex information.
# index contains simplexes, and the second dimension contains simplex
# information for that simplex. The information is:

# [type, tmlo, tmhi, point1, point2, point3, point4, neighbor1, neighbor2, 
#  neighbor3, neighbor4, simplexID]
def read_3simplices_from_data_file(filepath):
    "Reads in a *.3sx2p1 file."
    SIMPLEX_LIST_LENGTH = 12
    with open(filepath) as f:
        data_parsed_by_line = f.read().split('\n')
    header_data = data_parsed_by_line[0].split(' ')
    simplex_data_parsed_by_line = data_parsed_by_line[1:]
    simplex_data = []
    for i in simplex_data_parsed_by_line:
        simplex_data.append(i.split(' '))
    for i in range(len(simplex_data)):
        if len(simplex_data[i]) == SIMPLEX_LIST_LENGTH:
            for j in range(len(simplex_data[i])):
                simplex_data[i][j] = int(simplex_data[i][j])
        else: 
            simplex_data.pop(i)
    # JM: This line causes problems later. Commenting it out.
    #simplex_data = np.array(simplex_data)
    return [header_data, simplex_data]

# Reads in data from a movie file. Each element of the array is a
# single moment in "time" after making a move. The element is a list
# containing the 2-volume of each proper time slice in the
# simulation. This lets you make a plot of spatial volume as a
# function of time slice and simulation/thermalization time.
def read_movie_data(inputfile):
    "Reads in data from a movie file."
    return np.loadtxt(inputfile)

# Extracts the points of a 2 simplex from a 3 simplex datum (from a
# line in *.3sx2p1). Returns a list.
def extract_2simplex(simplex_datum):
    """Extracts the points of a 2 simplex from a 3 simplex datum 
    (from a line in *.3sx2p1). Returns a list."""
    
    # Defines simplex types for readiblity
    type13 = 1
    type31 = 3

    # The first element of a simplex datum list, is the type.
    simplex_type = simplex_datum[0]

    # If the simplex is of type (3,1), look at the first 3 points. If
    # it is of type (1,3), look at the last 3 points. In the unlikely
    # case the function was not passed a simplex datum, then return
    # "FAILURE" and print an error message.
    if simplex_type == type13:
        spacelike_subsimplex = simplex_datum[4:7]
    elif simplex_type == type31:
        spacelike_subsimplex = simplex_datum[3:6]
    else:
        print "This doesn't seem right..."
        return "FAILURE"
    return spacelike_subsimplex



# Given a spacetime output from read_3simplices_from_data_file (just
# the second element!) and a time t, find all space-like 2-simplices
# at time t. The approach here is that each 2simplex is is the face of
# (3,1) or a (1,3) simplex. For the first time-slice, we look at
# (3,1)-simplices with tmlo at time=0. For each other time slice, we
# look at (1,3)-simplices with tmlo at t=(time-1). This function uses
# extract_2simplex
def get_spacelike_2simplices_at_t(spacetime,time):
    """Reads in the spacelike 2-simplices on a given timeslice, 
    given a spacetime array."""
    
    # Defines simplex types for readiblity
    type13 = 1
    type31 = 3

    # Generate a list of 2-simplices at the correct proper time. If
    # the proper time is 0, then look at (3,1)-simplices with
    # tlo=0. Otherwise look at (1,3)-simplices in the tlow = (t-1).
    if time == 0:
        # Generate list of (3,1)-simplices with at the proper time:
        list3simplices = [i for i in spacetime if i[0]==type31 and i[1]==time]
        # Generate list of 2-simplex faces of the (3,1) simplices at time 0.
        list2simplices = [extract_2simplex(i) for i in list3simplices]
    else:
        # Generate list of (1,3)-simplices with faces at the proper time:
        list3simplices = [i for i in spacetime if i[0]==type13 and i[2]==time]
        # Generate list of 2-simplex faces of the (1,3)-simplices at time.
        list2simplices = [extract_2simplex(i) for i in list3simplices]

    # Return the list of 2 simplices:
    return list2simplices


# Given a spacetime output read from 3simplices_from_data_file (the
# entire output, including the header line in element zero, not just
# the second element), extracts all the 2 simplices and generates a
# vector indexing the simplices by time-slice. Assumes no empty
# time-slices.
def get_all_sl2simplices(spacetime):
    """Reads in the entire spacetime from 
    a data file and generates all spacelike 2-simplices."""
    
    # The list of 2-simplices. Each 2-simplex is set as a list of 3
    # points. The output of this function is a data structure of
    # nested lists. The top level is indexed by time slice and
    # contains lists of simplices. Thus, we have something like:
    # [[simplex_1_at_time_0,[p4,p5,p6],...],...,[[simplex_1_at_time_n],...]]
    simplices_by_time_slice = []

    # First we need to figure out how many time slices there are. This
    # is different for periodic and fixed boundary conditions
    
    # Name boundary condition types, for clarity.
    open_bc = 'OPEN'
    periodic_bc = 'PERIODIC'

    # Extract boundary condition type from spacetime
    boundary_condition_type = spacetime[0][0]

    # Figure out the number of time slices based on the boundary
    # condition type. It's possible these could be different.
    if boundary_condition_type == open_bc:
        n_time_slices = int(spacetime[0][2]) + 1
    elif boundary_condition_type == periodic_bc:
        n_time_slices = int(spacetime[0][2])
    else: # In case something goes horribly wrong. 
        print "ERROR: Unknown boundary condition type!"
        return "FAILURE" 
    
    # For each time slice, extract a list of the 2-simplices there
    for i in range(n_time_slices):
        s_at_t = get_spacelike_2simplices_at_t(spacetime[1],i)
        simplices_by_time_slice.append(s_at_t)

    # Return the output
    return simplices_by_time_slice


# Given a set of 2-simplices indexed by proper time (i.e., the output
# of get all 2-simplices), generate 2-volume (i.e., area) indexed by
# proper time. You can use this to compare to movie data if you want.
def make_v_of_t(triangle_set):
    "Makes a list of spatial 2-volume as a function of proper time."
    return [len(i) for i in triangle_set]
    

# Plots 2-volume as a function of proper time. Takes the output of
# make_v_of_t.
def plot_v_of_t(volume_list,name,iteration):
    """Plots 2-volume as a function of proper time. Takes the output of
    make_v_of_t. 
    name = name of simulation
    iteration = number of spacetime in ensemble. Might be sweep# instead."""
    
    # Defines the plot
    vplot = plt.plot(volume_list, 'bo', volume_list, 'r-')

    # plot title is made of name+iteration
    plot_title = name+' '+str(iteration)
    
    # Labels and Titles
    plt.title(plot_title)
    plt.xlabel('Proper Time')
    plt.ylabel('2-Volume Per Time Slice')

    # Ensure the y range is appropriate 
    plt.ylim([np.min(volume_list)-.5,np.max(volume_list)+.5])

    # Turn on minor ticks
    plt.minorticks_on()

    # Show the plot
    plt.show()
    return 

# Animates the spacetime volume as a function of proper time. The
# animation parameter is the "sweep."

# JM: On my implimentation of matplotlib, this doesn't work right. I'm
# possibly missing dependencies.
def animate_spacetime_evolution(moviedata,filename=''):
    """Animates the spacetime volume as a function of proper time.
    Animation parameter is the sqeep."""

    # The range for the x-axis
    plotdomain = len(moviedata[-1])

    # The range for the y-axis
    plotrange = np.max(moviedata)

    # The figure to use
    fig = plt.figure()
   
    # The axis
    ax = fig.add_subplot(111)
    ax.set_xlim(0,plotdomain-1)   # Set the x-axis
    ax.set_ylim(0,plotrange)      # Set the y-axis

    # Titles and labels
    ax.set_xlabel('Proper Time')
    ax.set_ylabel('2-Volume')

    # The list that will contain each frame
    frames = []

    # Make the list of frames. The for loop is a clumsy C-like loop
    # because order matters.
    for i in len(moviedata):
        temp_plot = plt.plot(moviedata[i])
        frames.append((temp_plot,))
    
    # Define the animation
    ani = animation.ArtistAnimation(fig,frames)

    # Show the animation
    plt.show()

    # Asks if you want to save it. If you say yes. Save it.
    do_save = query_yes_no("Do you want to save the animation? ","no")
    if do_save:
        ani.save(filename)
        print "Okay. The filename is "+str(filename)+".mp4"
        return filename+".mp4"
    return

# This one visualizes volume as a function of proper time.
# TODO: Make more visualize spacetime functions.
def visualize_spacetime_v1(filepath):
    "Visualizes a spacetime using the above commands."
    
    # Read in the spacetime
    spacetime = read_3simplices_from_data_file(filepath)
    
    # Generate a clever name for the plot:
    bc_type = spacetime[0][0]
    if bc_type == 'OPEN':
        ts = int(spacetime[0][2])+1
    elif bc_type == 'PERIODIC':
        ts = int(spacetime[0][2])
    else:
        print "Unrecogniced boundary conditions."
        print "Setting number of time slices to "+spacetime[0][2]
        ts = int(spacetime[0][2])
    
    N22 = int(spacetime[0][12]) # Number of (2,2)-simplices
    N31 = int(spacetime[0][11]) # Number of (1,3)- and(3,1)-simplices
    
    v3 = N22 + N31 # Total number of (3,1)-simplices
    
    # Phase space
    k0 = float(spacetime[0][-3])
    k3 = float(spacetime[0][-2])

    top = spacetime[0][1]

    # Magic numbers extracted from file name format. May not work forever.
    start_date = filepath[-26:-7]
    plottitle = bc_type+': '+str(ts)+' time slices.'+' 3-volume='+str(v3)+'. k0='+str(k0)+' k3='+str(k3)

    # Make the plot
    triangle_list = get_all_sl2simplices(spacetime)
    volume_list = make_v_of_t(triangle_list)
    plot_v_of_t(volume_list,plottitle,1)

# If file is run from command line. Execute. Otherwise do
# not. Supports globbing.
if __name__ == "__main__":
    for i in sys.argv[1:]:
        visualize_spacetime_v1(i)
