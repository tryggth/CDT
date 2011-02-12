(declaim (optimize (speed 3)
		   (compilation-speed 0)
		   (safety 0)
		   (debug 0)))

(load "cdt2p1.lisp")
(initialize-t-slices-with-v-volume :num-time-slices 8
				   :target-volume 8192
				   :spatial-topology "S2"
				   :boundary-conditions "PERIODIC")
(setf k0 0.0)
(setf k3 0.63)
(setf NUM-SWEEPS 300000)
(generate-data)