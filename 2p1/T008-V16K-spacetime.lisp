(declaim (optimize (speed 3)
		   (compilation-speed 0)
		   (safety 0)
		   (debug 0)))

(load "cdt2p1.lisp")
;;(initialize-t-slices-with-v-volume :num-time-slices 8
;;				   :target-volume (* 16 1024)
;;				   :spatial-topology "S2"
;;				   :boundary-conditions "PERIODIC")
;;(setf k0 0.50)
;;(setf k3 0.70)
(setf NUM-SWEEPS 400000)

(with-open-file (infile "S2-PERIODIC-T008-V016387-0.5-0.7-0.02-000000001-000100000.3sx2p1"
			:direction :input)
  (load-spacetime-from-file infile))
(format t "starting data generation for ~A at ~A~%" (generate-filename 100001) (cdt-now-str))
(generate-data 100001)
(format t "finished data generation for ~A at ~A~%" (generate-filename 100001) (cdt-now-str))