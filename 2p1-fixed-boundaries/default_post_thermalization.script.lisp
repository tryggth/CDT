;;;; default_post_thermalization.script.lisp

;;;; Author: Jonah Miller (jonah.maxwell.miller@gmail.com)

;;;; Loads a thermalized file and calculates how many sweeps are
;;;; necessary given parallelization.

(load "cdt2p1")

;;;; Constants
(defvar *filename* "S2-OPEN-T064-V030853-1.0-0.75772-0.02--1-000000001-000050000-on-jonah-ftdubs-started2012-08-05-00-17-14.3sx2p1"
  "The file we will use.")
(defvar *num-cores* 28 "The number of cores we can run the simulation on.")

(setf SAVE-EVERY-N-SWEEPS 500) ; How different each element of the
			       ; ensemble should be.
(defvar *ensemble-size* 1000 "How many spacetimes we want total.")

(defvar *total-sweeps* (* *ensemble-size* SAVE-EVERY-N-SWEEPS) 
  "How many sweeps we need to get the ensemble we want.")
(defvar *sweeps-per-core* (ceiling (/ *total-sweeps* *num-cores*))
  "The number of sweeps used per core.")

(setf NUM-SWEEPS *sweeps-per-core*)

(with-open-file (f *filename*) (load-spacetime-from-file f))

(generate-data-v2)
