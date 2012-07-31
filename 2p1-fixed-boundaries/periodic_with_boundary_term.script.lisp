;;;; periodic_with_boundary_term.script.lisp
;;;; Jonah Miller (jonah.maxwell.miller@gmail.com)
;;;; Date: July 30, 2012

;;;; Script with periodic boundary conditions, but such that the
;;;; boundary term is kept track of the entire time. A test of whether
;;;; or not the action is correct.

(load "cdt2p1.lisp")

(setf NUM-SWEEPS 50000)

; Initialize the spacetime with periodic boundary conditions
(initialize-t-slices-with-v-volume :num-time-slices 64
				   :target-volume 30850
				   :spatial-topology "S2"
				   :boundary-conditions "PERIODIC")

; Tell the simulation that the boundary is NOT periodic.
(setf BCTYPE "OPEN")

; Set the B-vector accordingly
(set-b-vector (count-spacelike-links-at-time 0)
	      (count-simplices-in-sandwich-of-type (1- NUM-T) 0 2)
	      (count-simplices-in-sandwich-of-type (1- NUM-T) 0 1)
	      (count-spacelike-links-at-time 0)
	      (count-simplices-in-sandwich-of-type 0 1 2)
	      (count-simplices-in-sandwich-of-type 0 1 3))

; Set the coupling constants
(set-k0-k3-alpha 1.0 0.75772 -1)

(generate-spacetime-and-movie-data)
