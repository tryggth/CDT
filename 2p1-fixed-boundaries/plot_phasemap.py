#!/usr/bin/env python2
"""
plot_phasemap.py
Author: Jonah Miller

Makes a 3D plot of 3-volume of a CDT simulation in 2+1
dimensions. Also places the critical surface on the plot for
comparison.

Takes two command line arguments: a docstring and a filename.
    If the command string contains certain letters, the function does
    more. Each letter adds a functionality.

    contains   function does
    'f'        complete 3D plot of 3-volume(k0,k3)
    's'        3D plot of 3-volume(k0,k3), but with values 0<=3-volume<=2.
    'c'        Plots the critical surface and fits k3 as a function of k0.
                     Outputs this information to a file.
 
The filename is *.phasemap and is some phasemap.

Example calls:
./plot_phasemap.py s T028_V030850.phasemap 
./plot_phasemap.py scf T028_V030850.phasemap 
"""

# Import modules
#--------------------------------------------------------------------------
import numpy as np                      # Arrays
from scipy.optimize import curve_fit    # Data fitting
import matplotlib.pyplot as plt         # Basic plot utilities
from mpl_toolkits.mplot3d import Axes3D # 3D utilities
from copy import copy                   # To copy arrays
import sys
#--------------------------------------------------------------------------


# Function definitions
#--------------------------------------------------------------------------
def resize_volume(data,target_volume):
    """" 
    Rescales the volume average and standard deviation data so
    volumes are v-avg/target-volume and v-std/target-volume.
    """
    data[2] = data[2]/float(target_volume)
    data[3] = data[3]/float(target_volume)
    data[4] = data[4]/float(target_volume)
    return data
    
def extract_target_volume(filename):
    "By looking at the file name, extract the target volume used."
    # The index the volume string starts at
    v_string_beginning = filename.find('-V0')+3
    # The indeexes the volume string ends at.  
    v_string_end = filename.find('-BI') 
    # The string containing the target volume
    v_string = filename[v_string_beginning:v_string_end]
    return int(v_string)  

def make_plot(data,plotname):
    """
    Makes a 3D plot of the data. Uses a scatter plot. Marks the
    critical surface.
    """

    # Initialize plot
    fig = plt.figure()
    ax = fig.add_subplot(111,projection='3d')
    ax.scatter(data[0],data[1],data[2])

    # Labels
    ax.set_xlabel('k0')
    ax.set_ylabel('k3')
    ax.set_zlabel('Average 3-Volume')
    ax.set_title(plotname)

    # Show the plot
    plt.show()

def extract_sane_data(normalized_data,minval,maxval):
    """
    Takes a normalized data set (so that target-volume = 1) and finds
    points around 1 (between minval and maxval).
    """
    normalized_data = normalized_data.transpose()
    new_data =  [i for i in normalized_data if minval<=i[2]<=maxval]
    return np.array(new_data).transpose()

def fit_polynomial(x,a,b,c):
    "A polynomial for polynomial least-squares fitting."
    return a * x**2 + b * x + c

def fit_to_critical_surface(data):
    "A polynomial fit to the shape of the critical surface as k3(k0)."
    return curve_fit(fit_polynomial,data[0],data[1])

def print_fit_file(filename,optimized_constants,k0min,k0max,k3min,k3max):
    """
    After fitting to the critical surface, this function makes a file
    containing information about the fit parameters so that the line
    can be reproduced.
    """
    outfile_name = filename+'.critical_surface_fit'
    with open(outfile_name,'w') as f:
        f.write('# Original file')
        f.write(filename+'\n')
        f.write('# k0 = {} - {}. k3 = {} - {}.\n'.format(k0min,
                                                         k0max,k3min,k3max))
        f.write('#k3:\n')
        f.write('a * k0 * k0 + b * k0 + c\n')
        f.write('a b c\n')
        f.write('{} {} {}\n'.format(optimized_constants[0],
                                    optimized_constants[1],
                                    optimized_constants[2]))
        

def complete_plot(filename):
    "Given a file, does the data analysis and generates a plot."
    # The phasemap data
    # k0 k3 volume-average volume-standard-deviation, total-change-in-volume
    data = np.loadtxt(filename).transpose() 
    # The target volume
    target_volume = extract_target_volume(filename)
    # Renormalize the data so that the everything is proportional to
    # target volume.
    data = resize_volume(data,target_volume)
    # Plot
    plotname = filename[filename.index('T0'):]
    make_plot(data,plotname)

def sane_plot(filename):
    """
    Given a file, does the data analysis, and generates a more
    readable plot.
    """
    
    sanevals = [0,2] # Sane [min,max] values to look at 3-volume for.
    # The phasemap data
    # k0 k3 volume-average volume-standard-deviation, total-change-in-volume
    data = np.loadtxt(filename).transpose() 
    # The target volume
    target_volume = extract_target_volume(filename)
    # Renormalize the data so that the everything is proportional to
    # target volume.
    data = resize_volume(data,target_volume)
    data = extract_sane_data(data,sanevals[0],sanevals[1])
    # Plot
    plotname = filename[filename.index('T0'):]
    make_plot(data,plotname)

def plot_critical_surface(filename):
    "Plots the critical surface of a file only."
    tolerance = 0.1 # How close to the target volume we're willing to accept.
    # The phasemap data
    # k0 k3 volume-average volume-standard-deviation, total-change-in-volume
    data = np.loadtxt(filename).transpose() 
    # Minimum and maximum values of k0 and k3.
    k0min = min(data[0])
    k0max = max(data[0])
    k3min = min(data[1])
    k3max = max(data[1])
    # The target volume
    target_volume = extract_target_volume(filename)
    # Renormalize the data so that the everything is proportional to
    # target volume.
    data = resize_volume(data,target_volume)
    data = extract_sane_data(data,1-tolerance,1+tolerance)
    # Curve fitting
    popt, pcov = fit_to_critical_surface(data)
    k0_fitted = np.linspace(min(data[0]),max(data[0]),250)
    k3_fitted = fit_polynomial(k0_fitted,popt[0],popt[1],popt[2])
    fiteqn = "k3 = {}$k0^2$ +{}$k0$ + {}".format(popt[0],popt[1],popt[2])
    # Plot
    plotname = "Critical Surface.\nV: {}. Tolerance: {}".format(target_volume,tolerance) + '\n' + fiteqn
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plots = [ax.plot(data[0],data[1],'bo'),ax.plot(k0_fitted,k3_fitted,'r-')]
    ax.set_xlabel('k0')
    ax.set_ylabel('k3')
    ax.set_title(plotname)
    plt.show()
    # Output fit information to file
    print_fit_file(filename,popt,k0min,k0max,k3min,k3max)

def main(command_string,filename):
    """
    For use with command line input.

    If the command string contains certain letters, the function does
    more. Each letter adds a functionality.

    contains   function does
    'f'        complete 3D plot of 3-volume(k0,k3)
    's'        3D plot of 3-volume(k0,k3), but with values 0<=3-volume<=2.
    'c'        Plots the critical surface and fits k3 as a function of k0.
                     Outputs this information to a file.
    """
    if command_string.count('f') > 0:
        complete_plot(filename)
    if command_string.count('s') > 0:
        sane_plot(filename)
    if command_string.count('c') > 0:
        plot_critical_surface(filename)
#--------------------------------------------------------------------------


# Main command line input
if __name__ == "__main__":
    command_string = sys.argv[1]
    for filename in sys.argv[2:]:
        main(command_string,filename)
