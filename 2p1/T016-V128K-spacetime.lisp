(load "cdt2p1.lisp")
(initialize-t-slices-with-v-volume :num-time-slices 16
				   :target-volume (* 128 1024)
				   :spatial-topology "S2"
				   :boundary-conditions "PERIODIC")
(setf k0 0.25)
(setf k3 0.66)
(setf NUM-SWEEPS 500000)
(generate-data)
;(with-open-file (infile "S2-PERIODIC-T016-V131224-1.0-0.78-0.02-000000001-000200000.3sx2p1" 
;			:direction :input)
;  (format t "opened the file, about to load~%")
;  (load-spacetime-from-file infile))
;(format t "finished loading the data~%")
;(generate-data 200001)
