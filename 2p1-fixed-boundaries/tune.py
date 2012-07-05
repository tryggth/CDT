#!/usr/bin/env python2
# Change the above line depending on your system

# tune.py
# Authors: Christian Anderson (original author)
#          Jonah Miller (jonah.miller@colorado.edu) updater 

# Given an alpha, and a k0, the system tunes k3 to ensure a critical
# surface. Outputs to a tuning data file.

# Import statements
#-----------------------------------------------------------------------------
from subprocess import Popen, PIPE, call
import time
import numpy
import sys
#-----------------------------------------------------------------------------

# Global variables
#-----------------------------------------------------------------------------
# The time we wait for SBCL to start up
waittime = 30

# The script written to "tuning.data."
script = """(load "cdt2p1.lisp")

(setf NUM-SWEEPS 50000)
(setf SAVE-EVERY-N-SWEEPS 100)

(initialize-t-slices-with-v-volume :num-time-slices 64
				   :target-volume (* 8 1024)
				   :spatial-topology "S2"
				   :boundary-conditions "OPEN"
                                   :initial-spatial-geometry {}
                                   :final-spatial-geometry {})

(set-k0-k3-alpha {} {} -1)

(generate-data-console)"""

# The name of the output file
outfile = 'tuning.data'

# The name of the tuning script
scriptname = 'tuning.script'

# Range for k3
k3min = 0.6
k3max = 1.2
k3range = numpy.arange(k3min,k3max,0.001)

# The implimentation of lisp in use on a given computer
lispcommand = ["sbcl","--dynamic-space-size","1024","--script",scriptname]

#k0 = 1.0
#-----------------------------------------------------------------------------

def main (k0,initial_geometry='"tetra.txt"',final_geometry='"tetra.txt"'):
    # Initial information
    with open(outfile,'w') as f:
        f.write("Data from tuning: k0=%d\n" % k0)

    # Generate tuned values
    for k3 in k3range:
        # Generate the script to run to test parameters
        with open('tuning.script','w') as f:
            f.write(script.format(initial_geometry,final_geometry,k0,k3))
        print k3
        # Run the script
        proc = Popen(lispcommand, stdout=PIPE)
        # proc = call(lispcommand, stdout=PIPE)
        time.sleep(waittime) # Wait a little bit for the system to start up
        proc.terminate() # Kill the procedure
        output = proc.stdout.read() # The to-console information produced
        print output # For anyone watching at home.
        output = output.split('\n') # Make a list
        # Only look at the part of the list we care about.
        output = filter(lambda x: x[0:5] == "start", output) 
        # Measure the standard deviation of the 3-volume. We want small values,
        # which implies a stabilized system.
        std = numpy.std(map(lambda x: int(x.split(' ')[14][:-1]), output))
        with open(outfile,'a') as f:
            if output:
                f.write("k3: %f, std: %f\n" % (k3, std))
            else:
                f.write("k3: %f, std: error\n" % (k3))
                
# Run the main program given input parameters
if __name__ == "__main__":
    num_commands = len(sys.argv) # Number of commands
    k0input = float(sys.argv[1]) # The choice for k0

    # If the number of commands is 2, run main for k0.
    if num_commands == 2:
        main(k0input)
    # Otherwise take initial and final geometry arguments
    else:
        main(k0input,sys.argv[2],sys.argv[3])
    
