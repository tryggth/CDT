(load "cdt2p1.lisp")
(initialize-t-slices-with-v-volume :num-time-slices 64
				   :target-volume (* 64 1024)
				   :spatial-topology "S2"
				   :boundary-conditions "PERIODIC")
(setf k0 0.25)
(setf k3 0.66)
(setf NUM-SWEEPS 500000)
(generate-movie-data)