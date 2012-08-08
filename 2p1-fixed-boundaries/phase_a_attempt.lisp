;;;; default_script.lisp
;;;; Jonah Miller (jonah.miller@colorado.edu)
;;;; Date: July 12, 2012

;;;; This is a script to test the behaviour of my bug fix and update
;;;; of David Kamensky's 2p1-fixed-boundaries code to work with Rajesh
;;;; Kommu's updated data structures and bug fixes.

;; Load the simulation environment
(load "cdt2p1.lisp")

;; Set the number of sweeps to 10000
(setf NUM-SWEEPS 50000)

;; Initialize a small spacetime for exploratory calculations
(initialize-t-slices-with-v-volume :num-time-slices          64
				   :target-volume            30850
				   :spatial-topology         "s2"
				   :boundary-conditions      "open"
                                   :initial-spatial-geometry "tetra.txt"
				   :final-spatial-geometry   "tetra.txt")

;; Initialize the coupling constants to a set that should get us
;; (hopefully) get us onto the critical surface in the k0-k3 phase diagram

;; alpha is set to -1 for consistency. Other values are possible.
(set-k0-k3-alpha 4.0 1.332 -1)

;; Take data!
(generate-data)
