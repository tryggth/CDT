#!/bin/bash

# Author: Jonah Miller (jonah.maxwell.miller@gmail.com)
# Time-stamp: <2012-10-31 22:22:41 (jonah)>

# This is a short script to write a script for each file given as
# input and submit it to the janus que.


# Your Variables
# ----------------------------------------------------------------------
# Where the program runs
WORKING_DIRECTORY=/projects/jonah/CDT/2p1-fixed-boundaries 
# Where output files go. Output works a little differently. Each job
# gets its own output folder. But this is the root directory.
OUTPUT_ROOT=/lustre/janus_scratch/jonah
# Shell you want to use. Unless you have a good reason, leave this alone.
SHEBANG='#!/bin/bash'
# The script command. Defaults to SBCL
SCRIPTCOMMAND="sbcl --dynamic-space-size 1220 --script"

# Torque Variables
# ----------------------------------------------------------------------
# Controls for Torque. Don't set the output directory. The script does
# that in a clever way for you.

TORQUE_PREFIX="#PBS " # Prefix. Useful later.

# Set this to be the number of elements of the torque array. You
# shouldn't have to do this, but something is buggy.
NUM_TORQUE_COMMANDS=11

# Put commands for Torque here. 
TORQUE[0]="-l walltime=00:08:01:00"
TORQUE[1]="-l mem=2gb"
TORQUE[2]="-l vmem=2gb"
TORQUE[3]="-l procs=1"
TORQUE[4]="-j oe"
TORQUE[5]="-m a"
TORQUE[6]="-q janus-long"
TORQUE[7]="-V"
TORQUE[8]="-d $WORKING_DIRECTORY"
TORQUE[9]="-w $WORKING_DIRECTORY"
TORQUE[10]="-M jonah.maxwell.miller@gmail.com"

# There's a torque command for handling output directories too, but
# the script handles that itself.

# Variables we want the script to keep (they're important).
PROJECT_HOME=/projects/jonah
SBCL=${PROJECT_HOME}/sbcl-1.0.58
SBCL_HOME=${SBCL}/lib/sbcl
SBCLBIN=${SBCL}/bin
PATH=.:$SBCLBIN:$PATH
# ----------------------------------------------------------------------

# Now, as we begin, we need to know what directory all these files are in.
SCRIPTDIR=$(pwd)
# Log file
MYLOGFILE="$(date "+%F-%H-%M-%S-%N").qsub.log"

# Now we build the scripts.
for i in $@; do # only works if the filenames have no spaces!
    # Print to the new script
    echo $SHEBANG > $i.sh
    echo '' >> $i.sh

    # Torque commands
    for (( j=0; j<$NUM_TORQUE_COMMANDS; j++ )); do
	echo $TORQUE_PREFIX ${TORQUE[$j]} >> $i.sh
    done
    # Defines the output directory
    echo "#PBS -o $OUTPUT_ROOT" >> $i.sh 

    # Ensure work done in the directory where the program binaries are.
    echo '' >> $i.sh
    echo "cd $WORKING_DIRECTORY" >> $i.sh
    echo '' >> $i.sh

    # Tell the script what to do.
    echo $SCRIPTCOMMAND $SCRIPTDIR/$i >> $i.sh

    # Submit the script.
    qsub $i.sh >> $SCRIPTDIR/$MYLOGFILE
done