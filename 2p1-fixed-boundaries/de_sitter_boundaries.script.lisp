(load "cdt2p1.lisp")

(setf NUM-SWEEPS 50000)

(initialize-t-slices-with-v-volume :num-time-slices 28
				   :target-volume 30850
				   :spatial-topology "S2"
				   :boundary-conditions "OPEN"
				   :initial-spatial-geometry "TS2-V64-k1.0-kkk0.7577200487786184d0-bottom.boundary"
				   :final-spatial-geometry "TS2-V64-k1.0-kkk0.7577200487786184d0-top.boundary")

(set-k0-k3-alpha 1.0 0.75772 -1)

(generate-spacetime-and-movie-data)
