(load "cdt2p1.lisp")
(initialize-t-slices-with-v-volume :num-time-slices 64
				   :target-volume (* 128 1024)
				   :spatial-topology "S2"
				   :boundary-conditions "PERIODIC")
(setf k0 0.25)
(setf k3 0.66)
(setf NUM-SWEEPS 200000)
;;(with-open-file (infile "S2-PERIODIC-T064-V524288-0.25-0.66-0.02-000000001-000250000.3sx2p1" 
;;			:direction :input)
;;  (format t "opened the file~%")
;;  (load-spacetime-from-file infile))
(generate-data)