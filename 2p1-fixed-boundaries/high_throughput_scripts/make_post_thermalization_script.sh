#!/bin/bash

# other stuff
realpath='/home/jonah/utils/realpath'

# Variables for the job
# ----------------------------------------------------------------------
SIMTIME=172800 # in seconds -- approximately 2 days
NUMSWEEPS=2000
# ----------------------------------------------------------------------

# Local environment data
# ----------------------------------------------------------------------
OUTDIR=/lustre/janus_scratch/jonah/
MYDIR=/projects/jonah/CDT/2p1-fixed-boundaries
infile=$($realpath $1)

cd $MYDIR # only work in the CDT directory.

echo '(load "cdt2p1")' > $1.script.lisp
echo '' >> $1.script.lisp
echo "(setf NUM-SWEEPS $NUMSWEEPS)" >> $1.script.lisp
echo '(setf SAVE-EVERY-N-SWEEPS 500)' >> $1.script.lisp
echo '' >> $1.script.lisp
printf '(defparameter *filename* "%s")\n' $infile >> $1.script.lisp
echo '' >> $1.script.lisp
printf '(setf *output-directory* "%s")' $OUTDIR >> $1.script.lisp
echo '' >> $1.script.lisp
echo '(with-open-file (f *filename*) (load-spacetime-from-file f))' >> $1.script.lisp
echo '' >> $1.script.lisp
(printf '(generate-data-in-time-v2 %s 1 "%s")\n' $SIMTIME $2) >> $1.script.lisp

echo $1.script.lisp
# ----------------------------------------------------------------------
