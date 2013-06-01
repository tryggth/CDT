Tools and code for Causal Dynamical Triangulations

See http://arxiv.org/abs/hep-th/0105267 for background

Github location: https://github.com/ucdavis/CDT

* dataviz.nb - Mathematica file to convert CDT output files into nice movies
* SurfaceOfResolution.m - Mathematica Algorithms for dataviz.nb
* gr.nb - Symbolically compute the Einstein tensor
* 2p1-fixed-boundaries-old - modified by David Kamensky to support open/fixed boundary conditions. Inherets bugs from older code. (LISP)
* 2p1-fixed-boundaries - Updated version of 2p1 to incorporate David Kamensky's modifications. Merger performed by Jonah Miller. (LISP)
* 2p1 - 2+1 spacetime code (LISP)
* 3p1 - 3+1 spacetime code (LISP)
* HL - Horava-Lifshitz code (LISP)
* Newton - 3+1 code with mass to test the Newtonian limit, and general rewrite of 3p1 code (Julia)
* runlisp.sh - A shell script to run Lisp jobs on Sun Grid Engine using 'qsub runlisp.sh'

2p1 and 3p1 code by Rajesh Kommu
HL (Horava-Lifshitz) code by Christian Anderson
2p1-fixed-boundaries by Jonah Miller

TODO: Complete Newton
TODO: Update Christian's code for fixed boundary conditions
DONE: Merge 2p1 and 2p1-fixed-boundaries. The updated code is called 2p1-fixed-boundaries. Works for open or periodic boundary conditions.
DONE: Fixed Christian's code to run on Sun Grid Engine
