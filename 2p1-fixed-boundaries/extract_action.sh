#!/usr/bin/env sh


# extract_action.sh
# Author: Jonah Miller (jonah.maxwell.miller@gmail.com)
# 
# This program takes a *.3sx2p1 file and returns the value of the
# un-normalized probability amplitude. In theory it should work with
# any version of the CDT code originally written by Rajesh Kommu on
# any platform. However, some changes to the source code are
# necessary. Primarily, the home directory of the CDT code will need
# to be changed. The CDT program itself must also have a module for
# calling the action: action_output.lisp.


# The home directory where the CDT stuff is
CDTHOME=~/CDT/2p1-fixed-boundaries

# The name of the script we'll use with sbcl
SCRIPT=extract_action.script.lisp

# Tell the user things are in the works
echo "Extracting the un-normalized probability amplitude."
echo "Please wait."

# Get the current directory
MYDIR=$(pwd)

# Get the filenames
PATHS=$(realpath $@)

# Change to the CDT home directory
cd $CDTHOME

# Start making the sbcl script
# import commands
echo '(load "cdt2p1.lisp")' > $SCRIPT
echo '(load "action_output.lisp")' >> $SCRIPT

# The list containing all file names
echo '(defparameter *filenamelist* nil)' >> $SCRIPT
for i in $PATHS; do
    echo "(push \"$i\" *filenamelist*)" >> $SCRIPT
done

# Get lisp to calculate mean and standard deviations
echo '(mean-and-std-dev-probability-amplitude *filenamelist* t)' >> $SCRIPT

# Call the script and saves the output string
OUTPUT=$(sbcl --dynamic-space-size 2000 --script $SCRIPT)

# Cleanup
rm $SCRIPT

# Return to the original directory and output
cd $MYDIR
echo $OUTPUT
