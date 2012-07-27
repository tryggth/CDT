;;;; 120625_test_first_run.lisp
;;;; Author: Jonah Miller (jonah.miller@colorado.edu)
;;;; Date: June 25, 2012

;;;; This is a script to test the behaviour of David Kemansky's code
;;;; now that I know the code works.

;; load the simulation environment
(load "cdt2p1.lisp")

;; Number of sweeps
(setf NUM-SWEEPS 10000)

;; Set the topology
(setf STOPOLOGY "S2")

;; Initialize a minimal triangulation
(initialize-s2-triangulation 64 "OPEN" "icosahedron.txt" "icosahedron.txt")
