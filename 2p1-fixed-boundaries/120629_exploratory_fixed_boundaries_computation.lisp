;;;; 120629_exploratory_fixed_boundaries_computation.lisp
;;;; Jonah Miller (jonah.miller@colorado.edu)
;;;; Date: June 25, 2012

;;;; This is a script to test the behaviour of my bug fix and update
;;;; of David Kamensky's 2p1-fixed-boundaries code to work with Rajesh
;;;; Kommu's updated data structures and bug fixes.

;; Load the simulation environment
(load "cdt2p1.lisp")

;; Set the number of sweeps to 2000
(setf NUM-SWEEPS 2000)

;; Initialize a small spacetime for exploratory calculations
(initialize-t-slices-with-v-volume :num-time-slices          64
				   :target-volume            10000
				   :spatial-topology         "s2"
				   :boundary-conditions      "open"
				   :initial-spatial-geometry "tetra.txt"
				   :final-spatial-geometry   "tetra.txt")

;; Initialize the coupling constants to a set that should get us
;; (hopefully) into the de-Sitter phase of the phase diagram.

;; Hopefully one of the following ordered pairs should work:
;; (k0,  k3)
;; (1.0, 0.753)
;; (2.0, 0.936)
;; (3.0, 1.108)

;; alpha is set to -1 for consistency. Other values are possible.
(set-k0-k3-alpha 2.0 0.936 -1)

;; Take data!
;; (generate-spacetime-and-movie-data)
