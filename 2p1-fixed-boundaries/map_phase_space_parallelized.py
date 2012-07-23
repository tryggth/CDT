#!/usr/bin/env python2
# Change the above line depending on your system

# tune.py
# Authors: Christian Anderson (original author)
#          Jonah Miller (jonah.miller@colorado.edu) updater 

# Given an alpha, and a k0, the system tunes k3 to ensure a critical
# surface. Outputs to a tuning data file.

# Import statements
#-----------------------------------------------------------------------------
import multiprocessing
import subprocess
import time
import numpy
import sys
import os
#-----------------------------------------------------------------------------

# Global variables
#-----------------------------------------------------------------------------
# The number of cores the computer has available for use.
num_cores = multiprocessing.cpu_count()

# The number of time slices used in the simulation
t_slices = 28
# Target 3 volume we want:
target_volume = 1024
# Initiali and final boundaries
initial_boundary = '"tetra.txt"'
final_boundary = '"tetra.txt"'

# The script used to initialize spaceitime
initialization_script = """
(load "cdt2p1.lisp")
(load "tuning_utilities.lisp")

(initialize-t-slices-with-v-volume :num-time-slices {}
				   :target-volume {}
				   :spatial-topology "S2"
				   :boundary-conditions "OPEN"
                                   :initial-spatial-geometry {}
                                   :final-spatial-geometry {})
""".format(t_slices,target_volume,initial_boundary,final_boundary)

# Number of thermalization sweeps used
thermalization_sweeps = 100

# The number of sweeps we use for data taking
data_sweeps = 5

# The remainder of the script... used to tell lisp how to take data.
data_taking_script = """
(set-k0-k3-alpha {} {} -1)

(run-thermalization-sweeps {})

(take-data-and-print-to-file {} {})
"""

# File name before the number
default_file_name = "tuning"

# How close we're willing to go
acceptable_range = target_volume * numpy.array([0.9,1.1])

# Range for k3
k3min = 0.75
k3max = 0.76
k3range = numpy.arange(k3min,k3max,0.001)

# Range for k0
k0min = 1.0
k0max = 1.1
k0range = numpy.arange(k0min,k0max,0.1)

# The implimentation of lisp in use on a given computer. You need to
# tell it the script too.
lispcommand = ["nice","sbcl","--dynamic-space-size","1024","--script"]
#-----------------------------------------------------------------------------

# Functions
#-----------------------------------------------------------------------------
def make_outfile_name():
    fn= 'T0{}-V0{}-BI{}-BF{}-k{}-{}-kkk{}-{}.phasemap'.format(t_slices,target_volume,initial_boundary[1:-1],final_boundary[1:-1],k0min,k0max,k3min,k3max)
    return fn

def make_script_name(k0, k3):
    "Generate a single file name for the startup script."
    return default_file_name + '-k' + str(k0) + '-kkk' + str(k3) + '.script'

def make_script_names(k0range,k3range):
    "Generate a large number of names for the startup scripts."
    filenames = []
    for k0 in k0range:
        for k3 in k3range:
            filenames.append(make_script_name(k0,k3))
    return filenames

def make_scripts(k0range,k3range,outfilename):
    "Generate the scripts we'll call with sbcl."
    scriptnames = []
    ofname = '"{}"'.format(outfilename)
    for k0 in k0range:
        for k3 in k3range:
            scriptname = make_script_name(k0,k3)
            script = initialization_script + data_taking_script.format(k0,k3,thermalization_sweeps,ofname,data_sweeps)
            with open(scriptname,'w') as f:
                f.write(script)
            scriptnames.append(scriptname)
    return scriptnames

def run_scripts(scriptlist):
    """
    Runs the scripts in scriptlist 
    When they all finish, returns true.
    """
    # Starts all the scripts
    procs = [] # the list of subprocesses we run
    running = 1 # 0 when the subprocesses are all done.
    for script in scriptlist:
        sbcl_call = lispcommand + [script]
        procs.append(subprocess.Popen(sbcl_call))
    while running:
        for proc in procs:
            proc.poll()
            if proc.returncode == 0:
                running = proc.returncode * running
            if proc.returncode == 1:
                print "ERROR. Program number {} failed.".format(procs.index(proc))
                return 1
    return 0
    
def make_slices(big_scriptlist):
    """
    Partitions big_scriptlist into smaller lists of script names with
    length less than or equal to the number of cores of the computer.
    """
    list_of_scriptlists = [] # will be our output
    incrementlist = range(0,len(big_scriptlist),num_cores) # how we increment
    for i in incrementlist:
        list_of_scriptlists.append(big_scriptlist[i:i+num_cores])
    return list_of_scriptlists

def fill_cores(big_scriptlist):
    """
    Calls scripts equal to the number of cores and no more. When those
    scripts are finished, calls new ones, to continually fill up all
    the cores.
    """
    core_filling_lists = make_slices(big_scriptlist)
    for i in core_filling_lists:
        run_scripts(i)

def make_outfile(outfilename):
    with open(outfilename,'w') as f:
        f.write('#k0 k3 mean std delta\n')

def main(k0range,k3range):
    outfilename = make_outfile_name()
    make_outfile(outfilename)
    print "Beginning phasemap!"
    print "k0 = {}-{}, k3 = {}-{}.".format(k0min,k0max,k3min,k3max)
    print "Bottom boundary = {}".format(initial_boundary)
    print "Top boundary = {}".format(final_boundary)
    print "Working..."
    print "..."
    big_scriptlist = make_scripts(k0range,k3range,outfilename)
    fill_cores(big_scriptlist)
    print "All done!"
    print "Removing tuning files..."
    for i in big_scriptlist:
        os.remove(i)
    print "Finished cleanup. Good luck with your simulations!"
    
                
# Run the main program given input parameters
#if __name__ == "__main__":
#    main(k0range,k3range)
