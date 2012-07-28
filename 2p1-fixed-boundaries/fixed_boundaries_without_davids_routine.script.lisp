(load "cdt2p1.lisp")

;; Sets sweeps
(setf NUM-SWEEPS 50000)
(setf SAVE-EVERY-N-SWEEPS 500)

;; Initialize spacetime with fixed boundary conditions, but don't
;; replace the initiali or final slices. Leave them as tetrahedra.
(initialize-t-slices-with-v-volume :num-time-slices 64
				   :target-volume 30850
				   :spatial-topology "S2"
				   :boundary-conditions "OPEN")

;; Set coupling constants
(set-k0-k3-alpha 1.0 0.7577 -1)

;; Gather data
(generate-spacetime-and-movie-data)
