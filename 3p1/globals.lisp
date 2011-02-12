;+---------------------------------------------------------------------------------------------------------+
;| globals.lisp --- all the parameters that might need to be accessed from multiple files          |
;+---------------------------------------------------------------------------------------------------------+
; we use the same random state during testing to verify that the bugs are being fixed.
#+ :sbcl (with-open-file (rndstt "../cdt-random-state-004.rndsbcl" :direction :input)
	   (setf *random-state* (read rndstt)))
#+ :ccl (with-open-file (rndstt "../cdt-random-state-001.rndccl" :direction :input)
	  (setf *random-state* (read rndstt)))

; uncomment the following code for "production" runs
;(setf *random-state* (make-random-state t)) ; a different sequence of random numbers each time

(defun reload-random-state ()
  #+ :sbcl (with-open-file (rndstt "../cdt-random-state-004.rndsbcl" :direction :input)
	     (setf *random-state* (read rndstt)))
  #+ :ccl (with-open-file (rndstt "../cdt-random-state-001.rndccl" :direction :input)
	    (setf *random-state* (read rndstt)))
  )

(defparameter *LAST-USED-3SXID* 0)
(defparameter *LAST-USED-4SXID* 0)
(defparameter *LAST-USED-POINT* 0)

(defmacro next-pt ()
  `(incf *LAST-USED-POINT*))
(defmacro set-last-used-pt (pt)
  `(setf *LAST-USED-POINT* ,pt))
(defmacro next-3simplex-id ()
  `(incf *LAST-USED-3SXID*))
(defmacro next-4simplex-id ()
  `(incf *LAST-USED-4SXID*))

;; 3simplex = (type tmlo tmhi (p0 p1 p2 p3))
(defun 3simplex->id-equality (3sx1 3sx2)
  (set-equal? (fourth 3sx1) (fourth 3sx2)))
(defun 3simplex->id-hashfn (3sx)
  (sxhash (sort (copy-list (fourth 3sx)) #'<)))

#+sbcl
(sb-ext:define-hash-table-test 3simplex->id-equality 3simplex->id-hashfn)

#+sbcl
(defparameter *3SIMPLEX->ID* (make-hash-table :test '3simplex->id-equality)) 

#+ccl
(defparameter *3SIMPLEX->ID* (make-hash-table :test '3simplex->id-equality 
					      :hash-function '3simplex->id-hashfn))

(defparameter *ID->3SIMPLEX* (make-hash-table))

(defparameter *ID->4SIMPLEX* (make-hash-table :test 'equal))

(defconstant 28MTYPE 0 "move type (2,8)")
(defconstant 46MTYPE 1 "move type (4,6)")
(defconstant 24v1MTYPE 2 "move type (2,4) 1st version")
(defconstant 24v2MTYPE 3 "move type (2,4) 2nd version")
(defconstant 33v1MTYPE 4 "move type (3,3) 1st version")
(defconstant 33v2MTYPE 5 "move type (3,3) 2nd version")
(defconstant 42v2MTYPE 6 "move type (4,2) 2nd version")
(defconstant 42v1MTYPE 7 "move type (4,2) 1st version")
(defconstant 64MTYPE 8 "move type (6,4)")
(defconstant 82MTYPE 9 "move type (8,2)")

(defparameter ATTEMPTED-MOVES (list 0 0 0 0 0 0 0 0 0 0) "number of attempted moves for each move type")
(defparameter SUCCESSFUL-MOVES (list 0 0 0 0 0 0 0 0 0 0) "number of successful moves for each move type")

(defun reset-move-counts ()
  (for (n 0 9)
       (setf (nth n ATTEMPTED-MOVES) 0 (nth n SUCCESSFUL-MOVES) 0)))

(defvar 28DF '(1 4 2 6 8 3 12 0 6 0) "the delta f vector for the (2,8) move")
(defvar 24DF '(0 0 1 0 4 0 2 3 0 2) "the delta f vector for the (2,4) move")
(defvar 46DF '(0 1 0 2 2 1 4 0 2 0) "the delta f vector for the (4,6) move")
(defvar 33DF '(0 0 0 0 0 0 0 0 0 0) "the delta f vector for the (3,3) move")
(defvar 64DF '(0 -1 0 -2 -2 -1 -4 0 -2 0) "the delta f vector for the (6,4) move")
(defvar 42DF '(0 0 -1 0 -4 0 -2 -3 0 -2) "the delta f vector for the (4,2) move")
(defvar 82DF '(-1 -4 -2 -6 -8 -3 -12 0 -6 0) "the delta f vector for the (8,2) move")

(defparameter N0 0 "number of points")
(defparameter N1-SL 0 "number of spacelike links")
(defparameter N1-TL 0 "number of timelike links")
(defparameter N2-SL 0 "number of spacelike triangles")
(defparameter N2-TL 0 "number of timelike triangles")
(defparameter N3-TL-31 0 "number of (1,3) + (3,1) timelike tetrahedra")
(defparameter N3-TL-22 0 "number of (2,2) timelike tetrahedra")
(defparameter N3-SL 0 "number of spacelike tetrahedra")
(defparameter N4-TL-41 0 "number of (1,4) + (4,1) timelike pentachora")
(defparameter N4-TL-32 0 "number of (2,3) + (3,2) timelike pentachora")
(defmacro N4 ()
  "total number of timelike 4simplices (pentachora)"
  `(+ N4-TL-41 N4-TL-32))
(defmacro N3 ()
  "total number of 3simplices (tetrahedra)"
  `(+ N3-SL N3-TL-31 N3-TL-22))

(defun set-f-vector (v1 v2 v3 v4 v5 v6 v7 v8 v9 v10)
  (setf N0 v1 N1-SL v2 N1-TL v3 N2-SL v4 N2-TL v5 
	N3-TL-31 v6 N3-TL-22 v7 N3-SL v8 N4-TL-41 v9 N4-TL-32 v10))

(defun update-f-vector (dv)
  (incf N0 (nth 0 dv))
  (incf N1-SL (nth 1 dv))
  (incf N1-TL (nth 2 dv))
  (incf N2-SL (nth 3 dv))
  (incf N2-TL (nth 4 dv))
  (incf N3-TL-31 (nth 5 dv))
  (incf N3-TL-22 (nth 6 dv))
  (incf N3-SL (nth 7 dv))
  (incf N4-TL-41 (nth 8 dv))
  (incf N4-TL-32 (nth 9 dv)))

(defparameter CURRENT-MOVE-NUMBER 0)
(defparameter CURRENT-MOVE-IDENTIFIER "UNKNOWN")
(defparameter STOPOLOGY "S2" "spatial slice topology --- S2 or T2")
(defparameter BCTYPE "UNKNOWN")

(defparameter NUM-T 666666 "number of time slices --- set to a non-zero value so (mod ts NUM-T) works")
(defparameter N-INIT 0 "initial volume of spacetime; we try to keep the volume close to this number")
(defparameter NUM-SWEEPS 0 "number of sweeps for which the simulation is to run")
(defparameter k4 0.0)
(defparameter k2 0.0)
(defparameter alpha 0.0)
(defparameter eps 0.02)


#+ :sbcl 
(defparameter RNDEXT ".rndsbcl" "used for storing the random state under SBCL compiler")
#+ :ccl 
(defparameter RNDEXT ".rndccl" "used for storing the random state under CCL compiler")
(defparameter 4SXEXT ".4sx3p1" "used for storing the timlike 4simplex (pentachora) information")

;; 6 7 8 9 10
;;-------------- t = 1
;; 1 2 3 4 5
;;-------------- t = 0
(defparameter S3-1/2-41 '((1 2 3 4 6) (2 3 4 5 7) (3 4 5 1 8) (4 5 1 2 9) (5 1 2 3 10)))
(defparameter S3-1/2-32 '((1 2 3 6 10) (2 3 4 6 7) (3 4 1 6 8) (4 1 2 6 9) (3 4 5 7 8)
			  (4 5 2 7 9) (5 2 3 7 10) (4 5 1 8 9) (5 1 3 8 10) (5 1 2 9 10)))
(defparameter S3-1/2-23 '((1 2 6 9 10) (2 3 6 7 10) (3 1 6 8 10) (3 4 6 7 8) (4 2 6 7 9)
			  (4 1 6 8 9) (4 5 7 8 9) (5 3 7 8 10) (5 2 7 9 10) (5 1 8 9 10)))
(defparameter S3-1/2-14 '((1 6 8 9 10) (2 6 7 9 10) (3 6 7 8 10) (4 6 7 8 9) (5 7 8 9 10)))
