;; An example script for using the HL CDT code
;; See the README for more information

(defpackage "CDT-HL"
  (:use "COMMON-LISP"))

(in-package "CDT-HL")

;; Loading files
(load "../utilities.lisp")
(load "globals.lisp")
(load "simplex.lisp")
(load "moves.lisp")
(load "initialization.lisp")
(load "montecarlo.lisp")

;; Defining the topology
(set-t-slices-with-v-volume :num-time-slices 64
			    :target-volume(* 8 1024)
			    :spatial-topology"S2"
			    :boundary-conditions"PERIODIC")

;; Setting *eps* to the minimum interesting value (for the sake of tuning)
(setf *eps* .02)

;; Optional last argument invokes tuning
(set-k0-k3-alpha-lambda-mu 1.0 0.0 -1.0 1.0 3.0 t)

;; Setting *eps* to the actual value of interest
(setf *eps* .05)

;; Initializing
(initialize)

;; Setting some parameters for my data collection
(setf SAVE-EVERY-N-SWEEPS 100)
(setf NUM-SWEEPS 50000)

;; Performing data collection
(generate-movie-data-console-v2)

