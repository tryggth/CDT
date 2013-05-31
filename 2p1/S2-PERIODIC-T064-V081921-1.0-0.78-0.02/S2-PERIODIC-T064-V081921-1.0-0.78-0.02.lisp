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
			
;; Fresh simulation section
(initialize-t-slices-with-v-volume 	:num-time-slices 64 
									:target-volume (* 80 1024) 
									:spatial-topology "S2" 
									:boundary-conditions "PERIODIC")
									
(set-k0-k3-alpha 1.0 0.78 -1.0)
(setf SAVE-EVERY-N-SWEEPS 100)
(setf NUM-SWEEPS 50000)
(generate-data)
;; end of fresh simulation section

;; resume simulation section
;; (setf SAVE-EVERY-N-SWEEPS 100)
;; (setf NUM-SWEEPS 50000)
;; (with-open-file (infile "???.3sx2p1") (load-spacetime-from-file infile))
;; (format t "finished loading data at ~A~%" (cdt-now-str))
;; (generate-data-v2 250001)