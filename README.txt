Tools and code for Causal Dynamical Triangulations

Github location: https://github.com/ucdavis/CDT

* dataviz.nb - Mathematica file to convert CDT output files into nice movies
* SurfaceOfResolution.m - Mathematica Algorithms for dataviz.nb
* gr.nb - Symbolically compute the Einstein tensor
* 2p1-fixed-boundaries-old - modified by David Kamensky to support open/fixed boundary conditions. Inherets bugs from older code. (LISP)
* 2p1-fixed-boundaries - Updated version of 2p1 to incorporate David Kamensky's modifications. Merger performed by Jonah Miller. (LISP)
* 2p1 - 2+1 spacetime code (LISP)
* 3p1 - 3+1 spacetime code (LISP)
* HL - Horava-Lifshitz code (LISP)
* Newton - 3+1 code with mass to test the Newtonian limit, written in Clojure, a Lisp running on the JVM (Clojure)
* runlisp.sh - A shell script to run Lisp jobs on Sun Grid Engine using 'qsub runlisp.sh'

2p1 and 3p1 code by Rajesh Kommu
HL (Horava-Lifshitz) code by Christian Anderson

TODO: Complete Newton
TODO: Update Christian's code for fixed boundary conditions
DONE: Merge 2p1 and 2p1-fixed-boundaries. The updated code is called 2p1-fixed-boundaries. Works for open or periodic boundary conditions.
DONE: Fixed Christian's code to run on Sun Grid Engine
