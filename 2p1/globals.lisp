;+---------------------------------------------------------------------------------------------------------+
;| cdt-2+1-globals.lisp --- all the parameters that might need to be accessed from multiple files          |
;+---------------------------------------------------------------------------------------------------------+
; we use the same random state during testing to verify that the bugs are being fixed.
#+ :sbcl (with-open-file (rndstt "../cdt-random-state-004.rndsbcl" :direction :input)
	   (setf *random-state* (read rndstt)))
#+ :ccl (with-open-file (rndstt "../cdt-random-state-001.rndccl" :direction :input)
	  (setf *random-state* (read rndstt)))


;; comment out the following three lines if mt19937 is not installed on your system
(require 'asdf)
(asdf:oos 'asdf:load-op :mt19937)
(setf mt19937:*random-state* (mt19937:make-random-state t))

(setf *random-state* (make-random-state t))

(defun reload-random-state ()
  #+ :sbcl (with-open-file (rndstt "../cdt-random-state-004.rndsbcl" :direction :input)
	     (setf *random-state* (read rndstt)))
  #+ :ccl (with-open-file (rndstt "../cdt-random-state-001.rndccl" :direction :input)
	    (setf *random-state* (read rndstt)))
  )

(defparameter *LAST-USED-2SXID* 0)
(defparameter *LAST-USED-3SXID* 0)
(defparameter *LAST-USED-POINT* 0)

(defmacro next-pt ()
  `(incf *LAST-USED-POINT*))
(defmacro set-last-used-pt (pt)
  `(setf *LAST-USED-POINT* ,pt))
(defmacro next-2simplex-id ()
  `(incf *LAST-USED-2SXID*))
(defmacro next-3simplex-id ()
  `(incf *LAST-USED-3SXID*))

(defun 2simplex->id-equality (2sx1 2sx2)
  (set-equal? (fourth 2sx1) (fourth 2sx2)))
(defun 2simplex->id-hashfn (2sx)
  (sxhash (sort (copy-list (fourth 2sx)) #'<)))

#+sbcl
(sb-ext:define-hash-table-test 2simplex->id-equality 2simplex->id-hashfn)

#+sbcl
(defparameter *2SIMPLEX->ID* (make-hash-table :test '2simplex->id-equality)) 

#+ccl
(defparameter *2SIMPLEX->ID* (make-hash-table :test '2simplex->id-equality 
					      :hash-function '2simplex->id-hashfn))

(defparameter *ID->2SIMPLEX* (make-hash-table))

(defparameter *ID->3SIMPLEX* (make-hash-table :test 'equal))

(defconstant 26MTYPE 0 "move type (2,6)")
(defconstant 23MTYPE 1 "move type (2,3)")
(defconstant 44MTYPE 2 "move type (4,4)")
(defconstant 32MTYPE 3 "move type (3,2)")
(defconstant 62MTYPE 4 "move type (6,2)")

(defparameter ATTEMPTED-MOVES (list 0 0 0 0 0) "number of attempted moves for each move type")
(defparameter SUCCESSFUL-MOVES (list 0 0 0 0 0) "number of successful moves for each move type")

(defun reset-move-counts ()
  (for (n 0 4)
       (setf (nth n ATTEMPTED-MOVES) 0 (nth n SUCCESSFUL-MOVES) 0)))

(defparameter N0 0 "number of points")
(defparameter N1-SL 0 "number of spacelike links")
(defparameter N1-TL 0 "number of timelike links")
(defparameter N2-SL 0 "number of spacelike triangles")
(defparameter N2-TL 0 "number of timelike triangles")
(defparameter N3-TL-31 0 "number of (1,3) + (3,1) timelike tetrahedra")
(defparameter N3-TL-22 0 "number of (2,2) timelike tetrahedra")
(defmacro N3 ()
  "total number of timelike 3simplices (tetrahedra)"
  `(+ N3-TL-31 N3-TL-22))

(defun set-f-vector (v1 v2 v3 v4 v5 v6 v7)
  (setf N0 v1 N1-SL v2 N1-TL v3 N2-SL v4 N2-TL v5 N3-TL-31 v6 N3-TL-22 v7))
(defun update-f-vector (dv)
  (incf N0 (nth 0 dv))
  (incf N1-SL (nth 1 dv))
  (incf N1-TL (nth 2 dv))
  (incf N2-SL (nth 3 dv))
  (incf N2-TL (nth 4 dv))
  (incf N3-TL-31 (nth 5 dv))
  (incf N3-TL-22 (nth 6 dv)))

(defparameter DF26 '(1 3 2 2 6 4 0))
(defparameter DF62 '(-1 -3 -2 -2 -6 -4 0))
(defparameter DF44 '(0 0 0 0 0 0 0))
(defparameter DF23 '(0 0 1 0 2 0 1))
(defparameter DF32 '(0 0 -1 0 -2 0 -1))

(defparameter CURRENT-MOVE-IDENTIFIER "UNKNOWN")
(defparameter CURRENT-MOVE-NUMBER 0)
(defparameter STOPOLOGY "unknown" "spatial slice topology --- S2 or T2")
(defparameter BCTYPE "unknown" "boundary conditions --- PERIODIC or OPEN")
(defparameter SAVE-EVERY-N-SWEEPS 100 "save every 100 sweeps by default")
(defparameter NUM-T 666666 "number of time slices --- set to a non-zero value so (mod ts NUM-T) works")
(defparameter N-INIT 0 "initial volume of spacetime; we try to keep the volume close to this number")
(defparameter NUM-SWEEPS 0 "number of sweeps for which the simulation is to run")
(defparameter k0 0.0)
(defparameter k3 0.0)
(defparameter eps 0.02)


#+ :sbcl 
(defparameter RNDEXT ".rndsbcl" "used for storing the random state under SBCL compiler")
#+ :ccl 
(defparameter RNDEXT ".rndccl" "used for storing the random state under CCL compiler")
(defparameter 3SXEXT ".3sx2p1" "used for storing the parameters and 3simplex information")
(defparameter 2SXEXT ".2sx2p1" "used for storing the timlike 2simplex (triangle) information")
(defparameter PRGEXT ".prg2p1" "used for keeping track of the progress of a simulation run")
(defparameter MOVEXT ".mov2p1" "used for storing the movie data information")

;; STOPOLOGY-BCTYPE-NUMT-NINIT-k0-k3-eps-startsweep-endsweep
(defun generate-filename (&optional (start-sweep 1) (end-sweep (+ start-sweep NUM-SWEEPS -1)))
  (format nil "~A-~A-T~3,'0d-V~6,'0d-~A-~A-~A-~9,'0d-~9,'0d" 
	  STOPOLOGY BCTYPE NUM-T N-INIT k0 k3 eps start-sweep end-sweep))

(defvar action nil)

(defun action-S1xS2 (num-0 num-3)
  (+ (- (* k3 num-3) (* k0 num-0)) (* eps (abs (- num-3 N-INIT)))))

;; 5, 6, 7, 8
;;------------ t=1
;; 1, 2, 3, 4
;;------------ t=0

;(defparameter S2-1/2-31 '((1 2 3 5) (2 3 4 6) (3 4 1 7) (4 1 2 8)))
;(defparameter S2-1/2-22 '((1 2 5 8) (2 3 5 6) (3 1 5 7) (3 4 6 7) (4 2 6 8) (4 1 7 8)))
;(defparameter S2-1/2-13 '((1 5 7 8) (2 5 6 8) (3 5 6 7) (4 6 7 8)))
;(defparameter S2-SPATIAL-TRIANGLES-ON-T0 '((1 2 3) (2 3 4) (3 4 1) (4 1 2)))
;(defparameter S2-SPATIAL-TRIANGLES-ON-T1 '((5 6 7) (6 7 8) (7 8 5) (8 5 6)))
(defparameter S2-1/2-31 '((1 2 3 13) (1 3 4 17) (1 4 5 21) (1 5 6 25) 
			  (1 6 2 29) (2 3 10 14) (2 6 9 30) (2 9 10 31) 
			  (3 4 11 18) (3 10 11 15) (4 5 7 22) (4 7 11 19) 
			  (5 6 8 26) (5 7 8 23) (6 8 9 27) (12 7 8 24) 
			  (12 8 9 28) (12 9 10 32) (12 10 11 16) (12 11 7 20)))

(defparameter S2-1/2-22 '((2 3 13 14) (1 3 13 17) (3 4 17 18) (3 11 15 18) 
			  (3 10 14 15) (1 4 17 21) (4 5 21 22) (4 7 19 22) 
			  (4 11 18 19) (1 5 21 25) (5 6 25 26) (5 8 23 26) 
			  (5 7 22 23) (1 6 25 29) (6 2 29 30) (6 9 27 30) 
			  (6 8 26 27) (10 11 15 16) (11 12 16 20) (11 7 19 20)
			  (7 12 20 24) (7 8 23 24) (8 12 24 28) (8 9 27 28) 
			  (9 12 28 32) (9 10 31 32) (2 10 31 14) (2 9 30 31) 
			  (1 2 29 13) (10 12 32 16)))

(defparameter S2-1/2-13 '((1 33 13 17) (1 33 17 21) (1 33 21 25) (1 33 25 29) 
			  (1 33 29 13) (2 34 13 14) (2 34 14 31) (2 34 30 31) 
			  (2 34 29 30) (2 34 29 13) (3 35 13 17) (3 35 17 18) 
			  (3 35 18 15) (3 35 15 14) (3 35 14 13) (4 36 17 21) 
			  (4 36 21 22) (4 36 22 19) (4 36 19 18) (4 36 18 17)
			  (5 37 21 25) (5 37 25 26) (5 37 26 23) (5 37 23 22) 
			  (5 37 22 21) (6 38 25 29) (6 38 29 30) (6 38 30 27) 
			  (6 38 27 26) (6 38 26 25) (7 39 19 22) (7 39 22 23) 
			  (7 39 23 24) (7 39 24 20) (7 39 20 19) (8 40 23 26) 
			  (8 40 26 27) (8 40 27 28) (8 40 28 24) (8 40 24 23) 
			  (9 41 27 30) (9 41 30 31) (9 41 31 32) (9 41 32 28) 
			  (9 41 28 27) (10 42 31 14) (10 42 14 15) (10 42 15 16) 
			  (10 42 16 32) (10 42 32 31) (11 43 15 18) (11 43 18 19) 
			  (11 43 19 20) (11 43 20 16) (11 43 16 15) (12 44 16 20) 
			  (12 44 20 24) (12 44 24 28) (12 44 28 32) (12 44 32 16)))
