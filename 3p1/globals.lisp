(declaim (optimize (speed 3)
		   (compilation-speed 0)
		   (debug 0)
		   (safety 0)))
;; globals.lisp --- all the parameters that might need to be accessed from 
;; multiple files

(setf *random-state* (make-random-state t))

(defparameter *LAST-USED-4SXID* 0)
(defparameter *RECYCLED-4SX-IDS* '())
(defparameter *LAST-USED-POINT* 0)
(defparameter *LAST-USED-S3SXID* 0)

(defmacro next-pt ()
  `(incf *LAST-USED-POINT*))
(defmacro set-last-used-pt (pt)
  `(setf *LAST-USED-POINT* ,pt))
(defmacro next-s3simplex-id ()
  `(incf *LAST-USED-S3SXID*))
(defmacro next-4simplex-id ()
  `(if (null *RECYCLED-4SX-IDS*)
       (incf *LAST-USED-4SXID*)
       (pop *RECYCLED-4SX-IDS*)))
(defmacro recycle-4simplex-id (sxid)
  `(push ,sxid *RECYCLED-4SX-IDS*))

;;------------------------------------------------------------------------------
;; timelike subsimplices have the form (type tmlo (p0 p1 ...))
(defun tlsubsx->id-hashfn (tlsx)
  (sxhash (sort (copy-list (third tlsx)) #'<)))
(defun tlsubsx->id-equality (tlsx1 tlsx2)
  (set-equal? (third tlsx1) (third tlsx2)))
(sb-ext:define-hash-table-test tlsubsx->id-equality tlsubsx->id-hashfn)
(defparameter *TL3SIMPLEX->ID* (make-hash-table :test 'tlsubsx->id-equality))
(defparameter *TL2SIMPLEX->ID* (make-hash-table :test 'tlsubsx->id-equality))
(defparameter *TL1SIMPLEX->ID* (make-hash-table :test 'tlsubsx->id-equality))

;; spacelike subsimplices have the form (tslice (p0 p1 ...))
(defun slsubsx->id-hashfn (slsx)
  (sxhash (sort (copy-list (second slsx)) #'<)))
(defun slsubsx->id-equality (slsx1 slsx2)
  (set-equal? (second slsx1) (second slsx2)))
(sb-ext:define-hash-table-test slsubsx->id-equality slsubsx->id-hashfn)
(defparameter *SL3SIMPLEX->ID* (make-hash-table :test 'slsubsx->id-equality))
(defparameter *SL2SIMPLEX->ID* (make-hash-table :test 'slsubsx->id-equality))
(defparameter *SL1SIMPLEX->ID* (make-hash-table :test 'slsubsx->id-equality))
;;-----------------------------------------------------------------------------
(defparameter *ID->SPATIAL-3SIMPLEX* (make-hash-table))
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

(defparameter ATTEMPTED-MOVES (list 1 1 1 1 1 1 1 1 1 1) 
  "number of attempted moves for each move type")
(defparameter SUCCESSFUL-MOVES (list 1 1 1 1 1 1 1 1 1 1) 
  "number of successful moves for each move type")

(defun reset-move-counts ()
  (for (n 0 9)
       (setf (nth n ATTEMPTED-MOVES) 0 (nth n SUCCESSFUL-MOVES) 0)))

(defun accept-ratios ()
  (format nil "[~A ~A ~A ~A ~A ~A ~A ~A ~A ~A]"
	  (* 100.0 (/ (nth 0 SUCCESSFUL-MOVES) (nth 0 ATTEMPTED-MOVES)))
	  (* 100.0 (/ (nth 1 SUCCESSFUL-MOVES) (nth 1 ATTEMPTED-MOVES)))
	  (* 100.0 (/ (nth 2 SUCCESSFUL-MOVES) (nth 2 ATTEMPTED-MOVES)))
	  (* 100.0 (/ (nth 3 SUCCESSFUL-MOVES) (nth 3 ATTEMPTED-MOVES)))
	  (* 100.0 (/ (nth 4 SUCCESSFUL-MOVES) (nth 4 ATTEMPTED-MOVES)))
	  (* 100.0 (/ (nth 5 SUCCESSFUL-MOVES) (nth 5 ATTEMPTED-MOVES)))
	  (* 100.0 (/ (nth 6 SUCCESSFUL-MOVES) (nth 6 ATTEMPTED-MOVES)))
	  (* 100.0 (/ (nth 7 SUCCESSFUL-MOVES) (nth 7 ATTEMPTED-MOVES)))
	  (* 100.0 (/ (nth 8 SUCCESSFUL-MOVES) (nth 8 ATTEMPTED-MOVES)))
	  (* 100.0 (/ (nth 9 SUCCESSFUL-MOVES) (nth 9 ATTEMPTED-MOVES)))))

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
(defmacro N2 ()
  "total number of 2-simplices (triangles)"
  `(+ N2-TL N2-SL))

(defun set-f-vector (v1 v2 v3 v4 v5 v6 v7 v8 v9 v10)
  (setf N0 v1 N1-SL v2 N1-TL v3 N2-SL v4 N2-TL v5 
	N3-SL v6 N3-TL-31 v7 N3-TL-22 v8 N4-TL-41 v9 N4-TL-32 v10))

(defun update-f-vector (dv)
  (incf N0 (nth 0 dv))
  (incf N1-SL (nth 1 dv))
  (incf N1-TL (nth 2 dv))
  (incf N2-SL (nth 3 dv))
  (incf N2-TL (nth 4 dv))
  (incf N3-SL (nth 5 dv))
  (incf N3-TL-31 (nth 6 dv))
  (incf N3-TL-22 (nth 7 dv))
  (incf N4-TL-41 (nth 8 dv))
  (incf N4-TL-32 (nth 9 dv)))

(defvar DF28 '(1 4 2 6 8 3 12 0 6 0) 
  "the delta f vector for the (2,8) move")
(defvar DF46 '(0 1 0 2 2 1 4 0 2 0) 
  "the delta f vector for the (4,6) move")
(defvar DF24 '(0 0 1 0 4 0 2 3 0 2) 
  "the delta f vector for the (2,4) move [same for v1 and v2]")
(defvar DF33 '(0 0 0 0 0 0 0 0 0 0) 
  "the delta f vector for the (3,3) move [same for v1 and v2]")
(defvar DF42 '(0 0 -1 0 -4 0 -2 -3 0 -2) 
  "the delta f vector for the (4,2) move [same for v1 and v2]")
(defvar DF64 '(0 -1 0 -2 -2 -1 -4 0 -2 0) 
  "the delta f vector for the (6,4) move")
(defvar DF82 '(-1 -4 -2 -6 -8 -3 -12 0 -6 0) 
  "the delta f vector for the (8,2) move")

(defparameter CURRENT-MOVE-NUMBER 0)
(defparameter CURRENT-MOVE-IDENTIFIER "UNKNOWN")
(defparameter STOPOLOGY "S3" "spatial slice topology --- always S3 in 3+1")
(defparameter BCTYPE "UNKNOWN" "boundary conditions --- PERIODIC or OPEN")

(defparameter SAVE-EVERY-N-SWEEPS 10)

(defparameter NUM-T 666666 
  "number of time slices --- set to a non-zero value so (mod ts NUM-T) works")
(defparameter N-INIT 0 
  "initial volume of spacetime; we try to keep the volume close to this number")
(defparameter NUM-SWEEPS 0 
  "number of sweeps for which the simulation is to run")
(defparameter KAPPA-0 0.0)
(defparameter DELTA 0.0)
(defparameter KAPPA-4 0.0)
(defparameter EPS 0.02)

(defparameter SIM-START-TIME (cdt-now-str) 
  "set again inside the generate methods for more accurate value")

(defvar 28MARKER 0.0)
(defvar 46MARKER 0.0)
(defvar 24v1MARKER 0.0)
(defvar 24v2MARKER 0.0)
(defvar 33v1MARKER 0.0)
(defvar 33v2MARKER 0.0)
(defvar 42v2MARKER 0.0)
(defvar 42v1MARKER 0.0)
(defvar 64MARKER 0.0)
(defvar 82MARKER 0.0)

(defvar action nil)

(defun damping (num41 num32)
  (* EPS (abs (- (+ num41 num32) N-INIT))))

(defun initialize-move-markers ()
  (setf 28MARKER 5.0) ;; 0
  (setf 46MARKER (+ 28MARKER 5.0)) ;;1
  (setf 24v1MARKER (+ 46MARKER 5.0));;2
  (setf 24v2MARKER (+ 24v1MARKER 5.0));;3
  (setf 33v1MARKER (+ 24v2MARKER 5.0));;4
  (setf 33v2MARKER (+ 33v1MARKER 5.0));;5
  (setf 42v2MARKER (+ 33v2MARKER 5.0));;6
  (setf 42v1MARKER (+ 42v2MARKER 5.0));;7
  (setf 64MARKER (+ 42v1MARKER 5.0));;8
  (setf 82MARKER (+ 64MARKER 5.0)));;9

(defun update-move-markers ()
  (let ((num-successful-moves (sum SUCCESSFUL-MOVES)))
    (setf 28MARKER (float (/ num-successful-moves 
			     (nth 28MTYPE SUCCESSFUL-MOVES))))
    (setf 46MARKER (float (/ num-successful-moves 
			     (nth 46MTYPE SUCCESSFUL-MOVES))))
    (setf 24v1MARKER (float (/ num-successful-moves 
			       (nth 24v1MTYPE SUCCESSFUL-MOVES))))
    (setf 24v2MARKER (float (/ num-successful-moves 
			       (nth 24v2MTYPE SUCCESSFUL-MOVES))))
    (setf 33v1MARKER (float (/ num-successful-moves 
			       (nth 33v1MTYPE SUCCESSFUL-MOVES))))
    (setf 33v2MARKER (float (/ num-successful-moves 
			       (nth 33v2MTYPE SUCCESSFUL-MOVES))))
    (setf 42v2MARKER (float (/ num-successful-moves 
			       (nth 42v2MTYPE SUCCESSFUL-MOVES))))
    (setf 42v1MARKER (float (/ num-successful-moves 
			       (nth 42v1MTYPE SUCCESSFUL-MOVES))))
    (setf 64MARKER (float (/ num-successful-moves 
			     (nth 64MTYPE SUCCESSFUL-MOVES))))
    (setf 82MARKER (float (/ num-successful-moves 
			     (nth 82MTYPE SUCCESSFUL-MOVES))))))

(defun select-move ()
  (let ((rndval (random 82MARKER))
	(mtype -1))
    (if (< rndval 28MARKER)
	(setf mtype 0)
	(if (< rndval 46MARKER)
	    (setf mtype 1)
	    (if (< rndval 24v1MARKER)
		(setf mtype 2)
		(if (< rndval 24v2MARKER)
		    (setf mtype 3)
		    (if (< rndval 33v1MARKER)
			(setf mtype 4)
			(if (< rndval 33v2MARKER)
			    (setf mtype 5)
			    (if (< rndval 42v2MARKER)
				(setf mtype 6)
				(if (< rndval 42v1MARKER)
				    (setf mtype 7)
				    (if (< rndval 64MARKER)
					(setf mtype 8)
					(setf mtype 9))))))))))))

(defun set-kappa0-delta-kappa4 (k0 D k4)
  (setf KAPPA-0 k0 DELTA D KAPPA-4 k4)
  (initialize-move-markers)
  (let ((c0 (* -1.0 (+ KAPPA-0 (* 6.0 DELTA))))
	(c41 (+ KAPPA-4 (* 2.0 DELTA)))
	(c32 (+ KAPPA-4 DELTA)))
    (setf (symbol-function 'action)
	  #'(lambda (n0 n41 n32)
	      (+ (* c0 n0) (* c41 n41) (* c32 n32))))))

(defparameter RNDEXT ".rndsbcl" 
  "used for storing the random state under SBCL compiler")
(defparameter 4SXEXT ".4sx3p1" 
  "used for storing the timlike 4simplex (pentachora) information")
(defparameter PRGEXT ".prg3p1" 
  "used for keeping track of the progress of a simulation run")
(defparameter MOVEXT ".mov3p1" 
  "used for storing the movie data information")
(defparameter S3SXEXT ".s3sx3p1" 
  "used for storing the spatial 3-simplex information")

;; STOPOLOGY-BCTYPE-NUMT-NINIT-KAPPA0-DELTA-KAPPA4-EPS-startsweep-endsweep-hostname-currenttime
(defun generate-filename (&optional (start-sweep 1) 
			  (end-sweep (+ start-sweep NUM-SWEEPS -1)))
  (format nil "~A-~A-T~3,'0d-V~6,'0d-~A-~A-~A-~A-~9,'0d-~9,'0d-on-~A-at-~A" 
	  STOPOLOGY BCTYPE NUM-T N-INIT KAPPA-0 DELTA KAPPA-4 EPS 
	  start-sweep end-sweep (hostname) (cdt-now-str)))

(defun generate-filename-v2 (&optional (ssweep 1) (csweep 0) 
			     (esweep (+ ssweep NUM-SWEEPS -1)))
  (format nil "~A-~A-T~3,'0d-V~6,'0d-~A-~A-~A-~A-~9,'0d-~9,'0d-~9,'0d-on-~A-start~A-curr~A" 
	  STOPOLOGY BCTYPE NUM-T N-INIT KAPPA-0 DELTA KAPPA-4 EPS 
	  ssweep csweep esweep (hostname) SIM-START-TIME (cdt-now-str)))

;; 6 7 8 9 10
;;-------------- t = 1
;; 1 2 3 4 5
;;-------------- t = 0
(defparameter S3-1/2-41 
  '((1 2 3 4 6) (2 3 4 5 7) (3 4 5 1 8) (4 5 1 2 9) (5 1 2 3 10)))
(defparameter S3-1/2-32 
  '((1 2 3 6 10) (2 3 4 6 7) (3 4 1 6 8) (4 1 2 6 9) (3 4 5 7 8)
    (4 5 2 7 9) (5 2 3 7 10) (4 5 1 8 9) (5 1 3 8 10) (5 1 2 9 10)))
(defparameter S3-1/2-23 
  '((1 2 6 9 10) (2 3 6 7 10) (3 1 6 8 10) (3 4 6 7 8) (4 2 6 7 9)
    (4 1 6 8 9) (4 5 7 8 9) (5 3 7 8 10) (5 2 7 9 10) (5 1 8 9 10)))
(defparameter S3-1/2-14 
  '((1 6 8 9 10) (2 6 7 9 10) (3 6 7 8 10) (4 6 7 8 9) (5 7 8 9 10)))

(defparameter N0-PER-SLICE 5 
  "number of points per spatial slice")
(defparameter N1-SL-PER-SLICE 10 
  "number of spacelike links per spatial slice")
(defparameter N1-TL-PER-SLICE 20 
  "number of timelike links per spatial slice")
(defparameter N2-SL-PER-SLICE 10 
  "number of spacelike triangles per spatial slice")
(defparameter N2-TL-PER-SLICE 60 
  "number of timelike triangles per spatial slice")
(defparameter N3-TL-31-PER-SLICE 20 
  "number of (3,1) timelike tetrahedra per spatial slice")
(defparameter N3-TL-13-PER-SLICE 20 
  "number of (1,3) timelike tetrahedra per spatial slice")
(defparameter N3-TL-22-PER-SLICE 30 
  "number of (2,2) timelike tetrahedra per spatial slice")
(defparameter N3-SL-PER-SLICE 5 
  "number of spacelike tetrahedra per spatial slice")
(defparameter N4-TL-41-PER-SLICE 5 
  "number of (1,4) timelike pentachora per spatial slice")
(defparameter N4-TL-14-PER-SLICE 5 
  "number of (4,1) timelike pentachora per spatial slice")
(defparameter N4-TL-32-PER-SLICE 10 
  "number of (2,3) timelike pentachora per spatial slice")
(defparameter N4-TL-23-PER-SLICE 10 
  "number of (3,2) timelike pentachora per spatial slice")

;;------------------------------------------------------------------------------
;; a vertex is (point tmslice)
;(defun vertex-hashfn (vtx) (sxhash vtx))
;(defun vertex-equality (vtx1 vtx2) (equal vtx1 vtx2))
;(sb-ext:define-hash-table-test vertex-equality vertex-hashfn)
;(defparameter *VERTEX-GRAPH* (make-hash-table :test 'vertex-equality))
