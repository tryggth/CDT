(proclaim '(optimize (speed 3)
		   (compilation-speed 0)
		   (safety 0)
		   (debug 0)))

(load "cdt2p1.lisp")

(initialize-T-slices-with-V-volume :num-time-slices 64 
				   :spatial-topology 'S2
				   :boundary-conditions 'periodic
				   :target-volume (* 80 1024))
(set-k0-k3-alpha 1.0 0.78 -1.0)
(setf NUM-SWEEPS 100000)
(generate-movie-data)

;(with-open-file (infile 
;		 "pbc-T8_V8492_eps0.02_kz2.0_kt0.985_sweeps10001to20000.data2p1"
;		 :direction :input)
;  (load-spacetime-from-file infile))
;(generate-data 20001)
