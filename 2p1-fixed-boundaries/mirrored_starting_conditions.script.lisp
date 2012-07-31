;;;; mirrored_starting_conditions.script.lisp
;;;; Jonah Miller (jonah.maxwell.miller@gmail.com)
;;;; Date: July 30, 2012

;;;; This is a script to test whether or not we can control where the
;;;; spatial extent of the simulation with fixed boundaries appears by
;;;; initialization.

(load "cdt2p1.lisp")

(setf NUM-SWEEPS 10000)

(with-open-file 
    (infile "~/CDT/2p1-fixed-boundaries/S2-OPEN-T064-V030851-1.0-0.7577-0.02--1-000000001-000050000-on-dewitt-started2012-07-27-17-04-43.time-inverted.3sx2p1")
  (load-spacetime-from-file infile))

(generate-spacetime-and-movie-data)
