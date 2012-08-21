"""
parallel_processing_tools.py

Author: Jonah Miller (jonah.maxwell.miller@gmail.com)

This file contains a few functions useful for parallel
processing. It's used by map_phase_space_parallelized.py,
extract_action.py, and others.
"""

# Import statements, if any
#-----------------------------------------------------------------------------
import multiprocessing
import subprocess
#-----------------------------------------------------------------------------


# Functions
#-----------------------------------------------------------------------------
def make_slices(big_scriptlist):
    """
    Partitions big_scriptlist into smaller lists of script names with
    length less than or equal to the number of cores of the
    computer. Useful for not running too many processes at once.
    """
    num_cores = multiprocessing.cpu_count()
    list_of_scriptlists = [] # This will be our output.
    incrementlist = range(0,len(big_scriptlist),num_cores) # How we increment.
    for i in incrementlist:
        list_of_scriptlists.append(big_scriptlist[i:i+num_cores])
    return list_of_scriptlists

def check_processes(process_list):
    """
    Takes an input list of running processes. When they all finish, return true.
    """
    running = 1 # 0 when the subprocesses are all done
    while running:
        for proc in process_list:
            proc.poll()
            if proc.returncode == 1:
                raise RuntimeError("Program " +
                                   "number " +
                                   "{}".format(process_list.index(proc)) +
                                   " failed.")
        running = bool(sum([int(proc.returncode) for proc in process_list]))
    return True

def start_processes(program_calls):
    """
    Calls program_call to measure the spacetime action for each call
    in program_calls. Returns a list of process objects that can be
    checked on.
    """
    processes = [subprocess.Popen(c,stdout=subprocess.PIPE) \
                     for c in program_calls]
    return processes
#-----------------------------------------------------------------------------
