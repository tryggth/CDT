(load "cdt2p1.lisp")
(initialize-t-slices-with-v-volume :num-time-slices 32
				   :target-volume (* 128 1024)
				   :spatial-topology "S2"
				   :boundary-conditions "PERIODIC")
(setf k0 1.0)
(setf k3 0.78)
(setf NUM-SWEEPS 250000)
(generate-data)