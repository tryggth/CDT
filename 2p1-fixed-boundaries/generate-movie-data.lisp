(load "cdt2p1.lisp")
(reset-spacetime)

(setf NUM-SWEEPS 100000)

(initialize-T-slices-with-V-volume :num-time-slices     64
				   :target-volume       80000
				   :spatial-topology    "s2"
				   :boundary-conditions "open"
				   :initial-spatial-geometry "tetra.txt"
				   :final-spatial-geometry   "tetra.txt")

(set-k-litL-alpha 0.20 5.00 -1.0)

;(generate-data-console)
(generate-spacetime-and-movie-data)