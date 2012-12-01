#!/bin/bash

# Author: Jonah Miller (jonah.maxwell.miller@gmail.com)
# Time-stamp: <2012-11-29 15:38:12 (jonah)>

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
# Sleeptime. Time between jobs started. Used to ensure unique file names.
SLEEPTIME="1s"
# processors per node. If you're requesting identical jobs.
PPN="12"


# Torque Variables
# ----------------------------------------------------------------------
# Controls for Torque. Don't set the output directory. The script does
# that in a clever way for you.

TORQUE_PREFIX="#PBS " # Prefix. Useful later.

# Set this to be the number of elements of the torque array. You
# shouldn't have to do this, but something is buggy.
NUM_TORQUE_COMMANDS=9

# Put commands for Torque here. 
TORQUE[0]="-l walltime=02:00:10:00"
TORQUE[1]="-l nodes=1:ppn=$PPN"
TORQUE[2]="-j oe"
TORQUE[3]="-m abe"
TORQUE[4]="-q janus-long"
TORQUE[5]="-V"
TORQUE[6]="-d $WORKING_DIRECTORY"
TORQUE[7]="-w $WORKING_DIRECTORY"
TORQUE[8]="-M jonah.maxwell.miller@gmail.com"

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
#SCRIPTDIR=$(pwd)
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

    # Tell the script what to do.
    echo '' >> $i.sh
    echo "cd $WORKING_DIRECTORY" >> $i.sh
    echo '' >> $i.sh

    # Tell it to generate the lisp script
    printf 'LISPSCRIPT=$(/home/jonah/utils/make_post_thermalization_script.sh %s $PBS_JOBID)\n\n' $(/home/jonah/utils/realpath $i) >> $i.sh

    # Tell it to run the lispscript $PPN times.nn
    for (( j=0; j<$PPN; j++)); do
	printf '%s ' $SCRIPTCOMMAND >> $i.sh
	printf '$LISPSCRIPT &\n' >> $i.sh
	echo "sleep $SLEEPTIME" >> $i.sh
    done

    echo '' >> $i.sh
    echo 'wait' >> $i.sh
    echo '' >> $i.sh

done
