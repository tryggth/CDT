;;;; default_script.lisp
;;;; Jonah Miller (jonah.miller@colorado.edu)
;;;; Date: July 12, 2012

;;;; This is a script to test the behaviour of my bug fix and update
;;;; of David Kamensky's 2p1-fixed-boundaries code to work with Rajesh
;;;; Kommu's updated data structures and bug fixes.

;; Load the simulation environment
(load "cdt2p1.lisp")

;; Calculate the number of sweeps we need per computer thread we use.
; bloch kusch cherenkov glaser moessbauer born kendall friedman taylor waals
(defparameter *num-cores* 10 "Number of cores we will run our code on.")
(defparameter *ensemble-size* 1000 "Number of spacetimes we want.")
(defparameter *spacetimes/core* (/ *ensemble-size* *num-cores*)
  "Spacetimes per core.")
;; Save every 500 sweeps
(setf SAVE-EVERY-N-SWEEPS 500)
; Sweeps per core
(setf NUM-SWEEPS (* SAVE-EVERY-N-SWEEPS *spacetimes/core*))


;; Initialize a small spacetime for exploratory calculations
(with-open-file (infile "S2-OPEN-T064-V030850-1.0-0.75772-0.02--1-000000001-000037500-000050000-on-moessbauer-start2012-07-19-16-50-06-curr2012-07-20-19-35-46.3sx2p1")
  (load-spacetime-from-file infile))

;; Take data!
(generate-data-v2)
