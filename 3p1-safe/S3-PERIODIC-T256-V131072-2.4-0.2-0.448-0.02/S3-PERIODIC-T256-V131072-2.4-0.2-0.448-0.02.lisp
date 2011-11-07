(declaim (optimize (safety 0) (speed 3)))

(load "../../utilities.lisp")
(load "../globals.lisp")
(load "../simplex.lisp")
(load "../moves.lisp")
(load "../initialization.lisp")
(load "../montecarlo.lisp")

(initialize-t-slices-with-v-volume :num-time-slices 256
				   :target-volume (* 128 1024)
				   :boundary-conditions "periodic")
(set-kappa0-delta-kappa4 2.4 0.2 0.448)

(setf NUM-SWEEPS 100000)
(setf SAVE-EVERY-N-SWEEPS 1000)
(format t "starting run ~A~%" (generate-filename))
(generate-spacetime-and-movie-data)
(format t "finished run ~A~%~%" (generate-filename))
(setf NUM-SWEEPS 400000)
(format t "starting run ~A~%" (generate-filename 100001))
(generate-data-v3 100001)
(format t "finished run ~A~%~%" (generate-filename 100001))