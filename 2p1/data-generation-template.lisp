(proclaim '(optimize (speed 3)
		   (compilation-speed 0)
		   (safety 0)
		   (debug 0)))

(load "../../utilities.lisp")
(load "../globals.lisp")
(load "../simplex.lisp")
(load "../moves.lisp")
(load "../initialization.lisp")
(load "../montecarlo.lisp")

;;(initialize-t-slices-with-v-volume :num-time-slices 64
;;				   :target-volume (* 80 1024)
;;				   :spatial-topology "S2"
;;				   :boundary-conditions "PERIODIC")

;;(set-k0-k3-alpha 1.0 0.78 -1.0)
(setf SAVE-EVERY-N-SWEEPS 100)
(setf NUM-SWEEPS 50000)
(with-open-file (infile "S2-PERIODIC-T064-V081921-1.0-0.78-0.02--1.0-000173351-000250000-on-kusch-started2010-04-23-16-19-48.3sx2p1")
  (load-spacetime-from-file infile))
(format t "finished loading data at ~A~%" (cdt-now-str))
(generate-data-v2 250001)