(load "cdt2p1.lisp")

;;function to run a small-volume test case
(defun test (new-k new-litL topology)
  (format t "~%::::::: spatial topology of ~A :::::::~%" topology)

  (reset-spacetime)
  (setf NUM-SWEEPS 5000)
  (initialize-T-slices-with-V-volume :num-time-slices     8
				     :target-volume       8192
				     :spatial-topology    topology
				     :boundary-conditions "PERIODIC")
  (set-k-litL-alpha new-k new-litL -1.0)
  (generate-data-console))

;; function to test a set of k/lambda parameters on both S2 and T2
(defun run-test-on-both-topologies (new-k new-litL)

  (format t "~%~%******************** TEST ******************************~%~%")
  (format t "======= PARAMETERS =======~%")
  (format t "   k     = ~$~%" new-k)
  (format t "   litL  = ~$~%" new-litL)
  (format t "   alpha = -1.0~%")
  (format t "==========================~%")

  ;;test S2 with set of parameters -------------------------------------------
  (test new-k new-litL "S2")

  ;;test T2 with same set of parameters ---------------------------------------
  (test new-k new-litL "T2"))



;; run some tests with above function

(run-test-on-both-topologies 0.2 5.10) ; should be stable for S2
;;(run-test-on-both-topologies 1.0 5.10) ; k too large (for S2)
;;(run-test-on-both-topologies 0.2 5.50) ; litL too large (for S2)

(quit)