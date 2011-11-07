(defpackage "CDT-HL"
  (:use "COMMON-LISP"))

(in-package "CDT-HL")

(load "../../utilities.lisp")
(load "../globals.lisp")
(load "../simplex.lisp")
(load "../moves.lisp")
(load "../initialization.lisp")
(load "../montecarlo.lisp")

(set-t-slices-with-v-volume :num-time-slices 64
			    :target-volume(* 100 1024)
			    :spatial-topology "S2"
			    :boundary-conditions "PERIODIC")

(set-k0-k3-alpha-lambda-mu 1.0 0.0 -1.0 0.75 0.5 t)

(setf *eps* .02)

(initialize)

(setf SAVE-EVERY-N-SWEEPS 100)
(setf NUM-SWEEPS 150000)

(format t "data generation begun at ~A~%" (cdt-now-str))
(generate-data-v2)
(format t "data generation completed at ~A~%" (cdt-now-str))
