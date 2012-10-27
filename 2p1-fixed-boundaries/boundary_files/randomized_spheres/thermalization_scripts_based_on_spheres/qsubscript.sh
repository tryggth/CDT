#!/bin/bash

# Lines starting with #PBS are treated by bash as comments, but interpreted by qsub
# as arguments.  For more details about usage of these arguments see "man qsub"

# Script names
THERM_SCRIPTS="AUTO_S2_OPEN_T028_V030850_B_S2_TA0100_i07_f010_started2012-10-24_13:11:26.110506.script.lisp AUTO_S2_OPEN_T028_V030850_B_S2_TA0300_i05_f010_started2012-10-24_13:11:57.802288.script.lisp AUTO_S2_OPEN_T028_V030850_B_S2_TA0100_i010_f010_started2012-10-24_13:11:26.110506.script.lisp AUTO_S2_OPEN_T028_V030850_B_S2_TA0100_i01_f010_started2012-10-24_13:11:26.110506.script.lisp AUTO_S2_OPEN_T028_V030850_B_S2_TA0100_i02_f010_started2012-10-24_13:11:26.110506.script.lisp AUTO_S2_OPEN_T028_V030850_B_S2_TA0100_i03_f010_started2012-10-24_13:11:26.110506.script.lisp AUTO_S2_OPEN_T028_V030850_B_S2_TA0100_i04_f010_started2012-10-24_13:11:26.110506.script.lisp AUTO_S2_OPEN_T028_V030850_B_S2_TA0100_i05_f010_started2012-10-24_13:11:26.110506.script.lisp AUTO_S2_OPEN_T028_V030850_B_S2_TA0100_i06_f010_started2012-10-24_13:11:26.110506.script.lisp AUTO_S2_OPEN_T028_V030850_B_S2_TA0100_i08_f010_started2012-10-24_13:11:26.110506.script.lisp AUTO_S2_OPEN_T028_V030850_B_S2_TA0100_i09_f010_started2012-10-24_13:11:26.110506.script.lisp AUTO_S2_OPEN_T028_V030850_B_S2_TA0300_i010_f010_started2012-10-24_13:11:57.802288.script.lisp AUTO_S2_OPEN_T028_V030850_B_S2_TA0300_i01_f010_started2012-10-24_13:11:57.802288.script.lisp AUTO_S2_OPEN_T028_V030850_B_S2_TA0300_i02_f010_started2012-10-24_13:11:57.802288.script.lisp AUTO_S2_OPEN_T028_V030850_B_S2_TA0300_i03_f010_started2012-10-24_13:11:57.802288.script.lisp AUTO_S2_OPEN_T028_V030850_B_S2_TA0300_i04_f010_started2012-10-24_13:11:57.802288.script.lisp AUTO_S2_OPEN_T028_V030850_B_S2_TA0300_i06_f010_started2012-10-24_13:11:57.802288.script.lisp AUTO_S2_OPEN_T028_V030850_B_S2_TA0300_i07_f010_started2012-10-24_13:11:57.802288.script.lisp AUTO_S2_OPEN_T028_V030850_B_S2_TA0300_i08_f010_started2012-10-24_13:11:57.802288.script.lisp AUTO_S2_OPEN_T028_V030850_B_S2_TA0300_i09_f010_started2012-10-24_13:11:57.802288.script.lisp AUTO_S2_OPEN_T028_V030850_B_S2_TA0500_i010_f010_started2012-10-24_13:12:38.321109.script.lisp AUTO_S2_OPEN_T028_V030850_B_S2_TA0500_i01_f010_started2012-10-24_13:12:38.321109.script.lisp AUTO_S2_OPEN_T028_V030850_B_S2_TA0500_i02_f010_started2012-10-24_13:12:38.321109.script.lisp AUTO_S2_OPEN_T028_V030850_B_S2_TA0500_i03_f010_started2012-10-24_13:12:38.321109.script.lisp AUTO_S2_OPEN_T028_V030850_B_S2_TA0500_i04_f010_started2012-10-24_13:12:38.321109.script.lisp AUTO_S2_OPEN_T028_V030850_B_S2_TA0500_i05_f010_started2012-10-24_13:12:38.321109.script.lisp AUTO_S2_OPEN_T028_V030850_B_S2_TA0500_i06_f010_started2012-10-24_13:12:38.321109.script.lisp AUTO_S2_OPEN_T028_V030850_B_S2_TA0500_i07_f010_started2012-10-24_13:12:38.321109.script.lisp AUTO_S2_OPEN_T028_V030850_B_S2_TA0500_i08_f010_started2012-10-24_13:12:38.321109.script.lisp AUTO_S2_OPEN_T028_V030850_B_S2_TA0500_i09_f010_started2012-10-24_13:12:38.321109.script.lisp qsubscript.sh"

## Local environment data
##----------------------------------------------------------------------
# Dot kit
. /curc/tools/utils/dkinit
# SBCL
export SBCL=${HOME}/sbcl-1.0.58
export SBCL_HOME=${SBCL}/lib/sbcl
export SBCLBIN=${SBCL}/bin
export PATH=.:$SBCLBIN:$PATH

# Use commands
use Torque
##----------------------------------------------------------------------


# Torque commands
##----------------------------------------------------------------------
# Set a walltime for the job. The time format is HH:MM:SS

# Run for 24 hours:
# #PBS -l walltime=24:00:00

# 

##----------------------------------------------------------------------

MYDIR=/projects/jonah/CDT/2p1-fixed-boundaries/

# Job submission
for i in $THERM_SCRIPTS; do
    qsub -r n -j -m abe -o /lustre/janus_scratch/jonah/$i.log -w $MYDIR -c shutdown -q janus-long -v SCBCLBIN,PATH sbcl --dynamic-space-size 2000 --script $i &;
done

EOF
