"""
monte_carlo.py

Time-stamp: <2012-10-18 11:29:04 (jonah)>

Author: Jonah Miller (jonah.maxwell.miller@gmail.com)

This file is part of the sphere_generator program which generates
spheres of (as close as possible to) uniform curvature for a given
surface area by monte carlo simulation.

This file contains the loops for the Metropolis-Hastings algorithm.
"""


### Dependencies
#-------------------------------------------------------------------------
import numpy as np
import scipy as sp
import random
import datetime
# Class data structures we need
import simplex_ancestors as sa
import simplex_descendants as sd
import state_manipulation as sm
import utilities as ut
import error_checking
import initialization
import state_tracking as st
import moves
#-------------------------------------------------------------------------



### CONSTANTS
###----------------------------------------------------------------------
output_suffix = ".boundary2p1" # The file ending for files generated.
output_prefix = "S2_" # The file prefix for files generated.
convergence_name = "_M-OPTIMAL"
###----------------------------------------------------------------------

####---------------------------------------------------------------------####
#                               Classes                                     #
####---------------------------------------------------------------------####
class metropolis:
    """
    The metropolis class is an ancestor class and never meant to be
    initialized.

    It contains methods to help descendant classes keep
    the surface area of the sphere constant.
    """
    def area_damping(self,surface_area):
        """
        A function that is closer to 1 the closer the number of
        triangles is to target_area.

        Both target_area and area_damping_strength are constants
        defined during initialization of descendant classes.
        """
        # The difference between the number of triangles we have and
        # the number of triangles we want.
        delta_triangles = abs(surface_area - self.target_area)
        # An exponential function we can treat as a probability.
        return np.exp(-abs(self.area_damping_strength)*delta_triangles)

    def loop_once(self):
        """
        Runs one iteration of the metropolis loop.
        """
        # Move data for a random move on a random simplex
        mdata = moves.try_random(moves.list_of_try_functions)

        # If a move is topologically acceptable, check the fitness
        # function, and possibly accept it. Otherwise, do nothing.
        if mdata:
            do_move = self.confirm_move(mdata)
            if do_move:
                moves.apply_move(mdata)

    def sweep(self,n=1,print_current_state=False):
        """
        Perform n*target_area iterations of the metropolis
        algorithm. Then if you asked for it, prints status.
        """
        for i in range(n):
            for j in range(self.target_area):
                self.loop_once()
            if print_current_state:
                print i
                print self.current_state
            self.current_sweep += 1

    def print_current_state(self):
        "Prints the current state. Just syntactic sugar."
        print self.current_state

    def reset_current_sweep(self, new_current_sweep=0):
        "Resets the current sweep to new_current_sweep."
        self.current_sweep = new_current_sweep

    def get_current_sweep(self):
        "Returns the current sweep."
        return self.current_sweep

#-------------------------------------------------------------------------


#-------------------------------------------------------------------------
class select_for_curvature(metropolis):
    """
    The select_for_curvature class is a descendant of metropolis. It
    IS meant to be initialized. It contains all methods necessary to
    run a monte carlo simulation to select for a given standard
    deviation of the local curvature of the sphere.
    """
    def __init__(self, target_area, target_std,
                 area_damping_strength = 0.8,
                 std_damping_strength = 0.8):
        """
        "Initializes a select_for_curvature class instance. This
        basically means defining the constants the instance needs.
        """
        self.area_damping_strength = area_damping_strength
        self.std_damping_strength = std_damping_strength
        self.target_area = target_area
        self.target_std = target_std
        self.current_state = st.sphere()
        self.current_sweep = 0

    def get_area_damping_strength(self):
        "Return area damping strength."
        return self.area_damping_strength

    def get_std_damping_strength(self):
        "Return damping for the standard deviation."
        return self.std_damping_strength

    def get_target_area(self):
        "Return target area."
        return self.target_area

    def get_target_std(self):
        "Return target standard deviation."
        return self.target_std

    def std_damping(self,standard_deviation):
        """
        A function is closer to 1 the closer the standard deviation of
        the curvature is to self.target_std.
        """
        # The difference between the standard deviation of the
        # curvature and the target standard deviation.
        delta_std = abs(standard_deviation - self.target_std)
        return np.exp(-abs(self.std_damping_strength)*delta_std)

    def fitness_function(self,damping_parameters):
        """
        How close our simulated sphere is to what we want. The ratio
        of fitness_function_new/fitness_function_old determines
        whether or not we accept a monte carlo move.
        """
        surface_area,standard_deviation = damping_parameters
        return self.area_damping(surface_area)\
            * self.std_damping(standard_deviation)

    def confirm_move(self,move_data):
        """
        Takes in the output of a try-move function and determines
        whether or not the move should be accepted. Uses self.fitness
        function.
        """
        # Shorter names
        new_sa = move_data.predicted_surface_area()
        new_std = move_data.predicted_curvature_std_dev()
        current_damping_parameters = [self.current_state.surface_area(),
                                      self.current_state.curvature_std()]

        # Probability of the new sphere
        p_new = self.fitness_function([new_sa,new_std])
        # Probability of the old sphere
        p_old = self.fitness_function(current_damping_parameters)

        # Probability of move acceptance (if >1, truncates to 1)
        p_acceptance = p_new/float(p_old)
        
        # If p_acceptance >= 1, accept the move. Otherwise accept the
        # move with probability p_acceptance.
        if random.random() < p_acceptance:
            return True
        else:
            return False

    def make_file_name_v1(self, final_sweep):
        """
        Makes a file name that matches the type of things the class
        optimizes for. Single file name without current sweep
        information.
        """
        now = datetime.datetime.today()
        nowstr = str(now).split(' ')[0] + '_' + str(now).split(' ')[1]
        outstring = output_prefix + "TA0" + str(self.target_area) \
            + "_STD0" + str(self.target_std) + "_f0" + str(final_sweep) \
            + "_started" + nowstr
        return outstring

    def make_file_name_v2(self, final_sweep,starttime):
        """
        Makes a file name that matches the type of things the class
        optimizes for. Contains current sweep information.
        """
        outstring = output_prefix + "TA0" + str(self.target_area) \
            + "_STD0" + str(self.target_std) + "_i0" \
            + str(self.current_sweep) \
            + "_f0" + str(final_sweep) \
            + "_started" + str(starttime)
        return outstring

    def make_file_name_v3(self,order_5_damping,order_6_damping):
        """
        Makes a file name that matches the type of things the class
        optimizes for. Contains damping information.
        """
        now = datetime.datetime.today()
        nowstr = str(now).split(' ')[0] + '_' + str(now).split(' ')[1]
        outstring = output_prefix \
            + "TA0" + str(self.target_area) \
            + "_STD0" + str(self.target_std)\
            + convergence_name\
            + "_V5D0"+str(order_5_damping)\
            + "_V6D0"+str(order_6_damping)\
            + "_started"+nowstr
        return outstring

    def make_file_name_v4(self, final_sweep,starttime,
                          order_5_damping,order_6_damping):
        """
        Makes a file name that matches the type of things the class
        optimizes for. Contains damping information and current sweep
        information.
        """
        outstring = output_prefix \
            + "TA0" + str(self.target_area) \
            + "_STD0" + str(self.target_std)\
            + "_i0" + str(self.current_sweep) \
            + "_f0" + str(final_sweep) \
            + convergence_name\
            + "_V5D0"+str(order_5_damping)\
            + "_V6D0"+str(order_6_damping)\
            + "_started"+starttime
        return outstring


    def make_progress_file_output(self,final_sweep,save_every_n_sweeps):
        "Returns the text for a progress file."
        outstring = "{} {} {} {} {}\n".format(self.target_area,
                                             self.area_damping_strength,
                                             self.target_std,
                                             self.std_damping_strength,
                                             save_every_n_sweeps)
        outstring += "{}/{}\n".format(self.current_sweep,final_sweep)
        return outstring

    def make_statistics_file_output(self,vertex_count):
        """
        The statistics file contains statistics on the simulated
        sphere in a human readable format.

        vertex_count should be a state_tracking.vertex_count instance,
        or you could get strange (not necessarily bad) behavior.
        """
        outstring = "# Surface-Area\tCurvature-STD\n"
        outstring+="{}\t{}\n".format(self.current_state.surface_area(),
                                     self.current_state.curvature_std())
        outstring+="# Vertex order information\n"
        outstring+=str(vertex_count)
        return outstring

#-------------------------------------------------------------------------
