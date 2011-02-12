(load "cdt2p1.lisp")
(initialize-t-slices-with-v-volume :num-time-slices 16
				   :target-volume (* 32 1024)
				   :spatial-topology "S2"
				   :boundary-conditions "PERIODIC")
(setf k0 0.25)
(setf k3 0.66)
(setf NUM-SWEEPS 1000000)
(generate-data)