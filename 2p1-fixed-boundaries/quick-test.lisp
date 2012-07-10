;;quick test to make sure that the program does not crash

(load "cdt2p1.lisp")

;;function to run a small-volume test case
(defun test (new-k new-litL topology boundary-conditions)
  (format t "~%::::::: spatial topology of ~A :::::::~%" topology)

  (reset-spacetime)
  (setf NUM-SWEEPS 1000)
  (initialize-T-slices-with-V-volume :num-time-slices     8
				     :target-volume       8192
				     :spatial-topology    topology
				     :boundary-conditions boundary-conditions
				     :initial-spatial-geometry "tetra.txt"
				     :final-spatial-geometry   "tetra.txt")
  (set-k-litL-alpha new-k new-litL -1.0)
  (generate-data-console))

;;a set of stable parameters for s2, at the given volume
(test 0.20 5.00 "s2" "open")

(quit)