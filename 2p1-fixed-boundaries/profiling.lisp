;;;; profiling.lisp
;;;; Author: Jonah Miller (jonah.maxwell.miller@gmail.com)

;;;; The puprose of this script is to profile the cdt2p1-fixed
;;;; boundaries algorithm and test how quickly each piece goes. Uses
;;;; sbcl statistical profiler.

(in-package :cl-user)
(require :sb-sprof)
(declaim (optimize speed))
(load "cdt2p1")

(setf NUM-SWEEPS 50000) ; Not useful. Just for convenience really.
(defparameter *run-time* (* 60 5) "walltime to run the tests.") ; in seconds

;; Initialize a typical spacetime
(initialize-t-slices-with-v-volume :num-time-slices 28
				   :target-volume 30850
				   :spatial-topology "S2"
				   :boundary-conditions "OPEN")

; set typical runtime conditions
(set-k0-k3-alpha 1.0 0.75772 -1)

;;;; CPU profiling

;;;; Take up to 1000 samples of the sweep algorithm and give a flat
;;;; table report at the end. When the loop ends, if the max samples
;;;; haven't been generated, too bad. A sample count will be printed
;;;; after each iteration.
(sb-sprof:with-profiling (:max-samples 10000
			  :report :flat
			  :loop nil)
  (generate-data-in-time *run-time*))

;;; Record call counts for functions defined on symbols in the CL-USER
;;; package
(sb-sprof:profile-call-counts "CL-USER")

;;;; Allocation profiling
(sb-sprof:with-profiling (:max-samples 1000
			  :mode :alloc
			  :report :flat)
  (generate-data-in-time *run-time*))
