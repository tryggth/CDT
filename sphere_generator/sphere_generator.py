#!/usr/bin/env python2

"""
sphere_generator.py

Time-stamp: <2013-03-13 13:39:24 (jonah)>

Author: Jonah Miller (jonah.maxwell.miller@gmail.com)

This file is the top-level interface for the sphere_generator program
which generates spheres of (as close as possible to) uniform curvature
for a given surface area by monte carlo simulation.

Call it by command line with:

./sphere_generator.py args

args can be a number of things. To run a simulation with all defaults
and a target surface are of sa, use:

./sphere_generator.py sa

You can also control surface area with the flag "--target-area". For
example, to set the surface area to 20, use:

./sphere_generator.py --target-area 20

In general, you control the simulation by calling it through the
command line with various flags. The value you want to set by calling
the flag comes after the flag. For instance, to set the surface area to 20
and the target standard deviation to 6:

./sphere_generator.py --target-std 6 --target-area 20

or

./sphere_generator.py 20 --target-std 6

Surface area is the only command line argument you can use without a
flag. If you use it without a flag, it must be the first argument.

The flags available are:

--select       (Sets the fitness function used for the simulation. The
                options are std and area. If you use "area", then the
                simulation only selects for spheres of a given target
                surface area. It ignores standard deviation. If you
                choose "std", the simulation selects for spheres of a given
                standard deviation and a given surface area. Defaults to
                "area".)

--target-area  (sets the surface area of the sphere. No default.
                You MUST set this value.)
                
--target-std   (sets the target standard deviation of curvature.
                Defaults to 0.)
                
--area-damping (a number between 0 and 1. How hard the simulation tries to
                stay at the target surface area. Defaults to 0.8.)
                
--std-damping  (a number between 0 and 1. How hard the simulation tries to
                stay at the target standard deviation of curvature.
                Defaults to 0.8.)
                
--initial      (Integer. The initial sweep number. Really only makes sense
                for continuing simulations. Defaults to 0.)
                
--final        (Integer. The final sweep number you want to go to.
                If --initial is set to 0, this is the total number of
                sweeps performed by the simulation. Defaults to 0.)
                
--save         (Integer. How often you save. If --save is set to 5,
                you save after every 5 sweeps. Defaults to 1.)
                
--many         (You don't need a value after this flag. It sets the
                simulation to save to a new file every time it saves
                the sphere. This is not the default behaviour.)
                
--one          (You don't need a value after this flag. It sets the
                simulation to save to a single file every time it saves.
                The simulation also writes to a progress file so that if
                it crashes, you can resume with no trouble.)

--micro        (You don't need a value after this flag. It controls when the
                simulation stops. If the flag is set, and the --one
                flag is set, the simulation stops when the sphere is
                close enough to flab based on the conditions for a
                micrscopically optimal--see the rest of the
                documentation--if the flag is set and the --many flag
                is set, the simulation runs until the final sweep but
                only saves when a sphere is close enough to the
                micriscopically optimal conditions. "close enough" is
                defined by the --v5 and --v6 flags. If the flag is not
                set, the simulation just runs until the final sweep.)

--exact        (This flag sets the simulation to save spheres only if they
                have the same surface area as the target area (plus or
                minus some margin of error which you set with the
                number after the flag). The simulation checks to see
                if it should save a sphere every save_every_n_sweeps
                sweeps, but only saves the sphere if it has the
                correct surface area. You set the margin of error with
                a number you give after the flag. If you set this
                flag, the --micro flag does nothing. The --final flag
                no longer counts the total number of sweeps. Rather,
                it counts the number of spheres you want to generate
                by the end of the simulation. With this flag set, the
                initial-sweep flag does nothing, and the current sweep
                in the file name represents the number of spheres
                saved so far in the simulation.)

--v5           (The margin of error for microscopically optimal spheres
                for vertices of order 5. A sphere is "close enough" to
                micrscopically optimal if
                12-v5 <= the number of vertices of order 5 <= 12+v5,
                and
                n-12-v6 <= the number of vertices of order 6 <= n-12+v6,
                where n is the total number of vertices, and the order of
                a vertex is the total number triangles attached to it.
                The default value is the target area (integer) divided by 10.)

--v6           (The margin of error for microscopically optimal spheres
                for vertices of order 6. A sphere is "close enough" to
                micrscopically optimal if
                12-v5 <= the number of vertices of order 5 <= 12+v5,
                and
                n-12-v6 <= the number of vertices of order 6 <= n-12+v6,
                where n is the total number of vertices, and the order of
                a vertex is the total number triangles attached to it.
                The default value is the target area (integer) divided by 10.)

--file         (File is a special flag. The value after it should be a
                filename ending in either ".boundary2p1" or "boundaryprg2p1".
                The simulation will load a sphere from that file. If the
                file is of type ".boundary2p1", it contains a list of
                triangles, each defined by 3 vertex numbers. This is all
                information needed to a load a sphere from file. All other
                command line flags still need to be called if you want to
                use them. If the file is of type ".boundaryprg2p1", the
                program will load the sphere file corresponding to that
                file name AND the parameters saved in the progress file.
                Other command line arguments don't do anything at all
                in this case.)

As alluded to in the tags, a sphere_generator.py call will produce one
of several types of outputs. If the --one flag is used, and the
--micro flag is not, the program will save to three files. One, ending
in ".boundary2p1", will contain exactly the format required to feed to
CDT/2p1-fixed-boundaries. It is of the form

((v1 v2 v3) (v1 v3 v4) ... )

where v1,v2,v3,v4,... are vertex numbers. Each space-separated list of
3 vertices is a triangle (i.e., a space-like 2-simplex). Another file
ends in ".boundaryprg2p1". If you select for standard deviation, it is
of the following format

target-area area-damping target-std std-damping save-every-n-sweeps
current-sweep/final-sweep

If you select for area, it is of the following format

target-area area-damping save-every-n-sweeps
current-sweep/final-sweep

The last file prints human-readable statistics on the sphere. It gives
surface area, standard deviation of curvature, and the number of
vertices of each order. It ends in ".boundarystatistics2p1"

If you select for standard deviation, the file name format will be

S2_TA0<target-area>_STD0<std>_f0<final-sweep>_started<start-date-and-time>.suffix

If you select for area, it will be

S2_TA0<target-area>_f0<final-sweep>_started<start-date-and-time>.suffix

If the --many flag is used, and the --micro flag is not, the program
will save a new sphere file every time it saves to file, so there's no
need for a progress file. There will still be statistics files. In
this case, if you select for curvature, the file names are of the
following form:

S2_TA0<target-area>_STD0<std>_io<current-sweep>_f0<final-sweep>_started<start-date-and-time>.suffix

If you select for standard deviation, it will be:

S2_TA0<target-area>_io<current-sweep>_f0<final-sweep>_started<start-date-and-time>.suffix

If the --micro flag is used, the files output will be the same, but the
formats will be

S2_TAO<target-area>_STD0<std>_M-OPTIMAL_V5D0<order_5_damping>_V6D0<order_6_damping>_started<start-time-and-date>.boundary2p1.suffix

and

S2_TAO<target-area>_STD0<std>_io<current-sweep>_f0<final-sweep>_M-OPTIMAL_V5D0<order_5_damping>_V6D0<order_6_damping>_started<start-time-and-date>.boundary2p1.suffix

respectively if you use select for curvature. If you use select for
area, they will be

S2_TAO<target-area>_M-OPTIMAL_V5D0<order_5_damping>_V6D0<order_6_damping>_started<start-time-and-date>.boundary2p1.suffix

and

S2_TAO<target-area>_io<current-sweep>_f0<final-sweep>_M-OPTIMAL_V5D0<order_5_damping>_V6D0<order_6_damping>_started<start-time-and-date>.boundary2p1.suffix

respectively. In this case, the program will either save to 1 file and
stop when the convergence conditions are met, or the program will run
final-sweep sweeps and only save a file when the file meets the
convergence conditions.

You can find some examples (with the .example suffix) in the folder
"output." By default, the files are saved to the pwd (which should
just be the sphere_generator folder), but you could change this
behavior by poking around in the code a little.
"""

### Dependencies
#-------------------------------------------------------------------------
import numpy as np
import scipy as sp
import random
import datetime
import sys
# Class data structures we need
import simplex_ancestors as sa
import simplex_descendants as sd
import state_manipulation as sm
import utilities as ut
import error_checking
import initialization as init
import state_tracking as st
import moves
import monte_carlo
import output
#-------------------------------------------------------------------------

## DEFAULTS
##----------------------------------------------------------------------
default_area_damping = 0.8
default_std_damping = 0.8
default_initial_sweep = 0
default_final_sweep = 0 # By default, the simulation makes a random
                        # sphere and does not perform metropolis at
                        # all.
default_save_every_n_sweeps = 1
default_target_std = 0 # By default we want a perfect sphere
default_algorithm = monte_carlo.select_for_area
##----------------------------------------------------------------------


## Data-gathering functions that require different inputs
##----------------------------------------------------------------------
sweep_based_functions = [output.gather_data_to_1_file,
                         output.gather_data_to_n_files]
microscopically_optimal_functions = [output.stop_at_microscopically_optimal,
                                     output.save_many_microscopically_optimal]
##----------------------------------------------------------------------


## CLASSES
## ----------------------------------------------------------------------
## These classes just hold parameters passed between the functions
## herein. It's really just an encapsulation thing, so that these
## functions can be changed later with less effort. I'm being very
## lazy with these functions because they're used in the UI but
## nowhere else, and they won't affect the actual simulation.
class parameters:
    "Holds parameters acquired through command line or progress files"
    def __init__(self,filename,target_area,area_damping_strength,
                 target_std,std_damping_strength,current_sweep,
                 final_sweep,save_every_n_sweeps,
                 v5damping, v6damping,
                 algorithm,
                 gather_data_function):
        self.filename = filename
        self.target_area = target_area
        self.area_damping_strength = area_damping_strength
        self.target_std = target_std
        self.std_damping_strength = std_damping_strength
        self.current_sweep = current_sweep
        self.final_sweep = final_sweep
        self.save_every_n_sweeps = save_every_n_sweeps
        self.v5damping = v5damping
        self.v6damping = v6damping
        self.algorithm = algorithm
        self.gather_data_function = gather_data_function

    def __str__(self):
        "For typecast into string. For debuggin."
        outstring = "filename: {}\n".format(self.filename)
        outstring += "Target area: {}\n".format(self.target_area)
        outstring += "Area Damping: {}\n".format(self.area_damping_strength)
        outstring += "Target std: {}\n".format(self.target_std)
        outstring += "std Damping: {}\n".format(self.std_damping_strength)
        outstring += "Current Sweep {}\n".format(self.current_sweep)
        outstring += "Final Sweep: {}\n".format(self.final_sweep)
        outstring += "Save every N sweeps: {}\n".format(self.save_every_n_sweeps)
        outstring += "Order 5 damping: {}\n".format(self.v5damping)
        outstring += "Order 6 damping: {}\n".format(self.v6damping)
        outstring += "Algorithm: {}\n".format(self.algorithm)
        outstring += "Gather data function: {}\n".format(self.gather_data_function)

        return outstring
        
        

## ----------------------------------------------------------------------


def get_progress_file_info(filename):
    """
    Reads a progress file and returns the information saved.

    Right now, only works for select_for_curvature
    """
    # Open the file and read it
    with open(filename,'r') as f:
        data = f.readlines()
    # Eat whitespace
    data = [s.lstrip().rstrip() for s in data]
    # Split the strings in each line 
    data[0] = data[0].split(' ')
    data[1] = data[1].split("/")
    # Now, based on position, we can start getting parameters
    current_sweep,final_sweep = [int(i) for i in data[1]]
    target_area = int(data[0][0])
    area_damping = float(data[0][1])
    save_every_n_sweeps = int(data[0][-1])
    if len(data[0]) == 3:
        target_std = default_target_std
        std_damping = default_std_damping
    else:
        target_std,std_damping = [float(s) for s in data[0][-2:-1]]
    # Gather data function is, of course, gather_data_to_1_file
    gather_data_function = output.gather_data_to_1_file
    # The filename we want to output is the same file but with
    # *.boundary2p1 at the end.
    filename = filename.rstrip(output.tracking_suffix)+output.output_suffix
    # Return the parameters as a tuple
    params =  parameters(filename,
                         target_area,area_damping,
                         target_std,std_damping,
                         current_sweep,final_sweep,
                         save_every_n_sweeps,
                         gather_data_function)
    return params

    

def parse_command_line_arguments(command_line_arguments):
    """
    Parses command line arguments
    """
    # First determine if a we're loading from a file
    filename = False
    if "--file" in command_line_arguments:
        index = command_line_arguments.index("--file")+1
        filename = command_line_arguments[index]
        if output.tracking_suffix in filename:
            # Assumes simulation in progress. So, if final_sweep ==
            # current_sweep, the simulation will load and then
            # immediately end.
            return get_progress_file_info(filename)
        if not (output.output_suffix in filename):
            raise ValueError("Can only load from *.boundaryprg2p1 or "
                             +"*.boundary2p1 files!")
        # If filename is of type *.boundary2p1, we assume its okay and
        # load from it. None of the other command line arguments
        # change.

    if "--select" in command_line_arguments:
        index = command_line_arguments.index("--select")+1
        if command_line_arguments[index] == "std":
            algorithm = monte_carlo.select_for_curvature
        elif command_line_arguments[index] == "area":
            algorithm = monte_carlo.select_for_area
        else:
            algorithm = default_algorithm
    else:
        algorithm = default_algorithm            

    if "--target-area" in command_line_arguments:
        index = command_line_arguments.index("--target-area")+1
        target_area = int(eval(command_line_arguments[index]))
    else:
        target_area = int(eval(command_line_arguments[0]))

    if "--target-std" in command_line_arguments:
        index = command_line_arguments.index("--target-std")+1
        target_std = float(eval(command_line_arguments[index]))
    else:
        target_std = default_target_std

    if "--area-damping" in command_line_arguments:
        index = command_line_arguments.index("--area-damping")+1
        area_damping_strength = float(eval(command_line_arguments[index]))
    else:
        area_damping_strength = default_area_damping
    if not 0 <= area_damping_strength <= 1:
        raise ValueError("Damping must be between 0 and 1.")

    if "--std-damping" in command_line_arguments:
        index = command_line_arguments.index("--std-damping")+1
        std_damping_strength = float(eval(command_line_arguments[index]))
    else:
        std_damping_strength = default_std_damping
    if not 0 <= area_damping_strength <= 1:
        raise ValueError("Damping must be between 0 and 1.")

    if "--initial" in command_line_arguments:
        index = command_line_arguments.index("--initial")+1
        initial_sweep = int(eval(command_line_arguments[index]))
    else:
        initial_sweep = default_initial_sweep

    if "--final" in command_line_arguments:
        index = command_line_arguments.index("--final")+1
        final_sweep = int(eval(command_line_arguments[index]))
    else:
        final_sweep = default_final_sweep

    if "--save" in command_line_arguments:
        index = command_line_arguments.index("--save")+1
        save_every_n_sweeps = int(eval(command_line_arguments[index]))
    else:
        save_every_n_sweeps = default_save_every_n_sweeps
    if save_every_n_sweeps < 1:
        raise ValueError("You must save at least every 1 sweeps!")

    if "--v5" in command_line_arguments:
        index = command_line_arguments.index("--v5")+1
        v5damping = int(eval(command_line_arguments[index]))
    else:
        v5damping = target_area/10

    if "--v6" in command_line_arguments:
        index = command_line_arguments.index("--v6")+1
        v6damping = int(eval(command_line_arguments[index]))
    else:
        v6damping = target_area/10
        
    if "--many" in command_line_arguments:
        if "--one" in command_line_arguments or "--exact" in command_line_arguments:
            raise ValueError("Contradictory input!")
        if "--micro" in command_line_arguments:
            gather_data_function = output.save_many_microscopically_optimal
        else:
            gather_data_function = output.gather_data_to_n_files
    elif "--one" in command_line_arguments:
        if "--many" in command_line_arguments or "--exact" in command_line_arguments:
            raise ValueError("Condtradictory input!")
        if "--micro" in command_line_arguments:
            gather_data_function = output.stop_at_microscopically_optimal
        else:
            gather_data_function = output.gather_data_to_1_file
    elif "--exact" in command_line_arguments:
        if "--many" in command_line_arguments or "--one" in command_line_arguments:
            raise ValueError("Contradictory input!")
        gather_data_function = output.generate_n_exact_spheres
        index = command_line_arguments.index("--exact")+1
        # In this case, v5damping is fitness_damping, as defined
        # in generate_n_exact_spheres
        v5damping = int(eval(command_line_arguments[index]))
    else:
        if "--micro" in command_line_arguments:
            gather_data_function = output.stop_at_microscopically_optimal
        else:
            gather_data_function = output.gather_data_to_1_file

    # return a class with all the info we need
    params = parameters(filename, target_area, area_damping_strength,
                        target_std, std_damping_strength,
                        initial_sweep, final_sweep,
                        save_every_n_sweeps,
                        v5damping, v6damping,
                        algorithm,
                        gather_data_function)
    return params
              

def start_simulation(parameters):
    """
    Takes an instance of the parameters class and starts a simulation
    based on those parameters. Uses the metropolis algorithm class
    type algorithm specified.
    """
    # If there's a sphere file to load, we need to load it. Otherwise,
    # we need to initialize the sphere.
    if parameters.filename:
        triangle_ids = init.load_sphere_from_file(parameters.filename)
    else:
        triangle_ids = init.initialize_sphere(parameters.target_area)

    # With the sphere initialized, we need to make an initialized
    # metropolis algorithm class instance.
    algorithm = parameters.algorithm
    if algorithm == monte_carlo.select_for_curvature:
        metro = algorithm(parameters.target_area,parameters.target_std,
                          parameters.area_damping_strength,
                          parameters.std_damping_strength)
    elif algorithm == monte_carlo.select_for_area:
        metro = algorithm(parameters.target_area,
                          parameters.area_damping_strength)
    else: #raise an error
        raise TypeError("The algorithm "+algorithm+" is not valid.")

    # Now gather data!
    parameters.gather_data_function(metro,
                                    parameters.final_sweep,
                                    parameters.v5damping,
                                    parameters.v6damping,
                                    parameters.current_sweep,
                                    parameters.save_every_n_sweeps)
                                    
    
def sphere_generator(command_line_arguments):
    """
    Takes command line arguments as input. Runs the sphere generator
    program.
    """
    start_simulation(parse_command_line_arguments(command_line_arguments))


# If this program is called from the command line, run it based on the
# command line arguments.
if __name__ == "__main__":
    sphere_generator(sys.argv[1:])
