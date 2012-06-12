;+---------------------------------------------------------------------------------------------------------+
;| cdt-2+1-globals.lisp --- all the parameters that might need to be accessed from multiple files          |
;+---------------------------------------------------------------------------------------------------------+
; we use the same random state during testing to verify that the bugs are being fixed.
#+ :sbcl (with-open-file (rndstt "../cdt-random-state-004.rndsbcl" :direction :input)
	   (setf *random-state* (read rndstt)))
#+ :ccl (with-open-file (rndstt "../cdt-random-state-001.rndccl" :direction :input)
	  (setf *random-state* (read rndstt)))

;; comment the following line to use a fixed seed from above
(setf *random-state* (make-random-state t))

(defun reload-random-state ()
  #+ :sbcl (with-open-file (rndstt "../cdt-random-state-004.rndsbcl" :direction :input)
	     (setf *random-state* (read rndstt)))
  #+ :ccl (with-open-file (rndstt "../cdt-random-state-001.rndccl" :direction :input)
	    (setf *random-state* (read rndstt)))
  )

(defparameter *LAST-USED-2SXID* 0)
(defparameter *LAST-USED-3SXID* 0)
(defparameter *RECYCLED-3SX-IDS* '())
(defparameter *LAST-USED-POINT* 0)
(defparameter *LAST-USED-S2SXID* 0)

(defmacro next-pt ()
  `(incf *LAST-USED-POINT*))
(defmacro set-last-used-pt (pt)
  `(setf *LAST-USED-POINT* ,pt))
(defmacro next-spatial-2simplex-id ()
  `(incf *LAST-USED-S2SXID*))
(defmacro next-2simplex-id ()
  `(incf *LAST-USED-2SXID*))
(defmacro next-3simplex-id ()
  `(if (null *RECYCLED-3SX-IDS*)
       (incf *LAST-USED-3SXID*)
       (pop *RECYCLED-3SX-IDS*)))
(defmacro recycle-3simplex-id (sxid)
  `(push ,sxid *RECYCLED-3SX-IDS*))

(defun 2simplex->id-equality (2sx1 2sx2)
  (set-equal? (fourth 2sx1) (fourth 2sx2)))
(defun 2simplex->id-hashfn (2sx)
  (sxhash (sort (copy-list (fourth 2sx)) #'<)))

;;macro to determine the euler characteristic of the spatial slices
;;assumes that the only two available are s2 and t2
(defmacro euler-char ()
  `(if (string= STOPOLOGY "S2") 2 1)) 
;; rkommu 2011-05-03 the second number above (chi for torus) should be 0

#+sbcl
(sb-ext:define-hash-table-test 2simplex->id-equality 2simplex->id-hashfn)

(defparameter *ID->SPATIAL-2SIMPLEX* (make-hash-table))

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

(defconstant ROOT2 (sqrt 2.0))
(defconstant KAPPA (/ (acos (/ 1 3)) pi))
(defconstant 6ROOT2 (* 6.0 ROOT2))
(defconstant 3KAPPAMINUS1 (- (* 3 KAPPA) 1))

(defparameter ATTEMPTED-MOVES (list 1 1 1 1 1) "number of attempted moves for each move type")
(defparameter SUCCESSFUL-MOVES (list 1 1 1 1 1) "number of successful moves for each move type")

(defun reset-move-counts ()
  (for (n 0 4)
       (setf (nth n ATTEMPTED-MOVES) 1 (nth n SUCCESSFUL-MOVES) 1)))

(defun accept-ratios ()
  (format nil "[~A ~A ~A ~A ~A]"
	  (* 100.0 (/ (nth 0 SUCCESSFUL-MOVES) (nth 0 ATTEMPTED-MOVES)))
	  (* 100.0 (/ (nth 1 SUCCESSFUL-MOVES) (nth 1 ATTEMPTED-MOVES)))
	  (* 100.0 (/ (nth 2 SUCCESSFUL-MOVES) (nth 2 ATTEMPTED-MOVES)))
	  (* 100.0 (/ (nth 3 SUCCESSFUL-MOVES) (nth 3 ATTEMPTED-MOVES)))
	  (* 100.0 (/ (nth 4 SUCCESSFUL-MOVES) (nth 4 ATTEMPTED-MOVES)))))

(defmacro percent-tv ()
  `(* 100.0 (/ (abs (- (N3) N-INIT)) N-INIT)))

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

;;make the deltas of the f-vector macros, so that
;;they can vary depending on whether the simplex on which
;;the move is performed is at a boundary

;;first, define macros to determine helpful things about
;;the position of a particular simplex
(defmacro in-upper-sandwich (sxid) 
  `(and (string= BCTYPE "OPEN") (= (3sx-tmhi (get-3simplex ,sxid))) NUM-T))
(defmacro in-lower-sandwich (sxid)
  `(and (string= BCTYPE "OPEN") (= (3sx-tmlo (get-3simplex ,sxid))) 0))
(defmacro in-either-boundary-sandwich (sxid)
  `(or (in-upper-sandwich ,sxid) (in-lower-sandwich ,sxid)))
(defmacro has-face-on-boundary (sxid)
  `(let* ((sx (get-3simplex ,sxid))  ;for fixed boundaries, this is unnecessary
	  (ty (3sx-type sx))
	  (th (3sx-tmhi sx))
	  (tl (3sx-tmlo sx)))
     (and (string= BCTYPE "OPEN")
	  (or (and (= ty 1)
		   (or (= th 0) (= th NUM-T)))
	      (and (= ty 3)
		   (or (= tl 0) (= tl NUM-T)))))))

;;use the above macros to help determine when special
;;DF's need to be applied
(defmacro DF26 (sxid)
  `(if (has-face-on-boundary ,sxid) ;for fixed boundaries, this is unnecessary
       (list 1 3 1 2 3 2 0)
       (list 1 3 2 2 6 4 0)))
(defmacro DF62 (sxid)
  `(if (has-face-on-boundary ,sxid) ;for fixed boundaries, this is unnecessary
       (list -1 -3 -1 -2 -3 -2 0)
       (list -1 -3 -2 -2 -6 -4 0)))
(defparameter DF44 (list 0 0 0 0 0 0 0))
(defparameter DF23 (list 0 0  1 0  2 0  1))
(defparameter DF32 (list 0 0 -1 0 -2 0 -1))


;;define some quantities of geometrical objects on the boundaries
(defparameter N1-SL-boundary 0)
(defparameter N3-22-boundary 0)
(defparameter N3-31-boundary 0)

;;the above three quantities compose the "b-vector"
(defun set-b-vector (n1 n2 n3)
  (setf N1-SL-boundary n1
	N3-22-boundary n2
	N3-31-boundary n3))
(defun update-b-vector (dv)
  (incf N1-SL-boundary (first  dv))
  (incf N3-22-boundary (second dv))
  (incf N3-31-boundary (third  dv)))

;;the deltas of the b-vector are dependent on where a move occurs
(defmacro DB26 (sxid)
  `(if (has-face-on-boundary ,sxid) ;for fixed boundaries, this is unnecessary
       (list 3 0 2) 
       (list 0 0 0)))
(defmacro DB62 (sxid)
  `(if (has-face-on-boundary ,sxid) ;for fixed boundaries, this is unnecessary
       (list -3 0 -2) 
       (list  0 0  0)))
(defparameter DB44 (list 0 0 0))
(defmacro DB23 (sxid)
  `(if (in-either-boundary-sandwich ,sxid)
       (list 0 1 0)
       (list 0 0 0)))
(defmacro DB32 (sxid)
  `(if (in-either-boundary-sandwich ,sxid)
       (list 0 -1 0)
       (list 0  0 0)))


(defparameter CURRENT-MOVE-IDENTIFIER "UNKNOWN")
(defparameter CURRENT-MOVE-NUMBER 0)
(defparameter STOPOLOGY "unknown" "spatial slice topology --- S2 or T2")
(defparameter BCTYPE "unknown" "boundary conditions --- PERIODIC or OPEN")
(defparameter SAVE-EVERY-N-SWEEPS 10 "save every 10 sweeps by default")
(defparameter NUM-T 666666 "number of time slices --- set to a non-zero value so (mod ts NUM-T) works")
(defparameter N-INIT 0 "initial volume of spacetime; we try to keep the volume close to this number")
(defparameter NUM-SWEEPS 0 "number of sweeps for which the simulation is to run")
(defparameter k0 0.0)
(defparameter k3 0.0)
(defparameter eps 0.02)
(defparameter SIM-START-TIME (cdt-now-str) "set again inside the generate methods for more accurate value")


#+ :sbcl 
(defparameter RNDEXT ".rndsbcl" "used for storing the random state under SBCL compiler")
#+ :ccl 
(defparameter RNDEXT ".rndccl" "used for storing the random state under CCL compiler")

(defparameter 3SXEXT ".3sx2p1" "used for storing the parameters and 3simplex information")
(defparameter PRGEXT ".prg2p1" "used for keeping track of the progress of a simulation run")
(defparameter MOVEXT ".mov2p1" "used for storing the movie data information")

(defvar action nil)

;;(defun action-S1xS2 (num-0 num-3)
;;  (+ (- (* k3 num-3) (* k0 num-0)) (* eps (abs (- num-3 N-INIT)))))

;;; wrsqrt is the "wick rotated" sqrt function. Basically wrsqrt(x) = -i*sqrt(-x) when x < 0 and not 
;;; i*sqrt(-x). So wrsqrt(-1) = -i
(defmacro wrsqrt (val)
  `(if (< ,val 0)
       (* -1 *cmplx-i* (sqrt (* -1 ,val)))
       (sqrt ,val)))

(defparameter *a* 1.0)
(defparameter *alpha* -1.0)
(defparameter *cmplx-i* #C(0.0 1.0)) ;; complex number i
(defparameter *minus-cmplx-i* #C(0.0 -1.0)) ;; complex number -i
(defparameter *k* 1.0)
(defparameter *litL* 1.0)

;; STOPOLOGY-BCTYPE-NUMT-NINIT-k0-k3-eps-alpha-startsweep-endsweep-hostname-currenttime
(defun generate-filename (&optional (start-sweep 1) (end-sweep (+ start-sweep NUM-SWEEPS -1)))
  (format nil "~A-~A-T~3,'0d-V~6,'0d-~A-~A-~A-~A-~9,'0d-~9,'0d-on-~A-started~A" 
	  STOPOLOGY BCTYPE NUM-T N-INIT k0 k3 eps *alpha* start-sweep end-sweep (hostname) (cdt-now-str)))

;; STOPOLOGY-BCTYPE-NUMT-NINIT-k0-k3-eps-alpha-startsweep-currsweep-endsweep-hostname-starttime-currenttime
(defun generate-filename-v2 (&optional (ssweep 1) (csweep 0) (esweep (+ ssweep NUM-SWEEPS -1)))
  (format nil "~A-~A-T~3,'0d-V~6,'0d-~A-~A-~A-~A-~9,'0d-~9,'0d-~9,'0d-on-~A-start~A-curr~A" 
	  STOPOLOGY BCTYPE NUM-T N-INIT k0 k3 eps *alpha* ssweep csweep esweep 
	  (hostname) SIM-START-TIME (cdt-now-str)))

(defvar 26MARKER 0.0)
(defvar 23MARKER 0.0)
(defvar 44MARKER 0.0)
(defvar 32MARKER 0.0)
(defvar 62MARKER 0.0)

;;the dihedral angle of an equilateral tetrahedron
(defparameter *theta* (acos 1/3))
(defparameter *sxvol* (/ 1 (* 6 (sqrt 2))))

;;define action as the _euclidean_ action, to be used directly as the weight in the
;;partition function

;;also, assumes alpha = -1, a = 1

(defun action (num1-sl num1-tl num3-22  num3-31
	       ;;some additional arguments to be passed for
	       ;;open boundary conditions
	       num1-sl-boundary
	       num3-22-boundary
	       num3-31-boundary)

  ;;define some variables to ease the expression-writing
  (let* ((num3  (+ num3-22 num3-31))
	 (num1  (+ num1-sl num1-tl)))

    ;;expression for the euclidean action
    (+ (* *litL* *sxvol* num3) 
       (* (- *k*) (- (* 2 pi (- num1 num1-sl-boundary)) 
		     (* (- (* 6 num3) (* 3 num3-31-boundary) num3-22-boundary) *theta*)))
       
       ;;boundary term:
       (* (- *k*) (- (* pi num1-sl-boundary) (* *theta* (+ (* 3 num3-31-boundary) num3-22-boundary)))))))


(defun damping (num3)
  (* eps (abs (- num3 N-INIT))))

(defun initialize-move-markers ()
  (setf 26MARKER 5.0)
  (setf 23MARKER (+ 26MARKER 5.0))
  (setf 44MARKER (+ 23MARKER 5.0))
  (setf 32MARKER (+ 44MARKER 5.0))
  (setf 62MARKER (+ 32MARKER 5.0)))

(defun update-move-markers ()
  (let ((num-successful-moves (apply #'+ SUCCESSFUL-MOVES)))
    (setf 26MARKER (float (/ num-successful-moves (nth 26MTYPE SUCCESSFUL-MOVES))))
    (setf 23MARKER (+ 26MARKER(float (/ num-successful-moves (nth 23MTYPE SUCCESSFUL-MOVES)))))
    (setf 44MARKER (+ 23MARKER(float (/ num-successful-moves (nth 44MTYPE SUCCESSFUL-MOVES)))))
    (setf 32MARKER (+ 44MARKER(float (/ num-successful-moves (nth 32MTYPE SUCCESSFUL-MOVES)))))
    (setf 62MARKER (+ 32MARKER(float (/ num-successful-moves (nth 62MTYPE SUCCESSFUL-MOVES)))))))

(defun select-move ()
  (let ((rndval (random 62MARKER))
	(mtype -1))
    (if (< rndval 26MARKER)
	(setf mtype 0)
	(if (< rndval 23MARKER)
	    (setf mtype 1)
	    (if (< rndval 44MARKER)
		(setf mtype 2)
		(if (< rndval 32MARKER)
		    (setf mtype 3)
		    (setf mtype 4)))))
    mtype))

(defun set-k0-k3-alpha (kay0 kay3 alpha)
  (setf k0 kay0 k3 kay3 *alpha* alpha)
  (setf *k* (/ k0 (* 2 *a* pi)))
  (setf *litL* (* (- k3 (* 2 *a* pi *k* 3KAPPAMINUS1)) (/ 6ROOT2 (* *a* *a* *a*))))
  (initialize-move-markers))

;;function analogous to set-k0-k3-alpha, but for k, litL, and alpha
;;(does not bother to set k0 and k3)
(defun set-k-litL-alpha (new-k new-litL new-alpha)
  (setf *k*     new-k
        *litL*  new-litL
        *alpha* new-alpha)
	
  ;;also set k0, k3 for compatibility with rest of code (saving/loading, filename generation, etc)
  (setf k0 (* *k* (* 2 *a* pi))
	k3 (+ (/ (* *litL* *a* *a* *a*) 6ROOT2) (* 2 *a* pi *k* 3KAPPAMINUS1)))

  (initialize-move-markers))



;;----------------------------------------------------------------------------------------------------------
;; initialization data
;;----------------------------------------------------------------------------------------------------------

;; 5, 6, 7, 8
;;------------ t=1
;; 1, 2, 3, 4
;;------------ t=0
(defparameter N0-PER-SLICE 4)
(defparameter N1-SL-PER-SLICE 6)
(defparameter N1-TL-PER-SLICE 12)
(defparameter N2-SL-PER-SLICE 4)
(defparameter N2-TL-PER-SLICE 24)
(defparameter N3-TL-13-PER-SLICE 4)
(defparameter N3-TL-22-PER-SLICE 6)
(defparameter N3-TL-31-PER-SLICE 4) ; here (3,1) means just (3,1), not (3,1)+(1,3)
(defparameter S2-1/2-31 '((1 2 3 5) (2 3 4 6) (3 4 1 7) (4 1 2 8)))
(defparameter S2-1/2-22 '((1 2 5 8) (2 3 5 6) (3 1 5 7) (3 4 6 7) (4 2 6 8) (4 1 7 8)))
(defparameter S2-1/2-13 '((1 5 7 8) (2 5 6 8) (3 5 6 7) (4 6 7 8)))

;; 13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44
;;------------------------------------------------------------------------------------------------- t=1
;; 1,2,3,4,5,6,7,8,9,10,11,12
;;------------------------------------------------------------------------------------------------- t=0
;;(defparameter N0-PER-2-SLICES 44)

;;(defparameter S2-1/2-31 '((1 2 3 13) (1 3 4 17) (1 4 5 21) (1 5 6 25) 
;;			  (1 6 2 29) (2 3 10 14) (2 6 9 30) (2 9 10 31) 
;;			  (3 4 11 18) (3 10 11 15) (4 5 7 22) (4 7 11 19) 
;;			  (5 6 8 26) (5 7 8 23) (6 8 9 27) (12 7 8 24) 
;;			  (12 8 9 28) (12 9 10 32) (12 10 11 16) (12 11 7 20)))
;;
;;(defparameter S2-1/2-22 '((2 3 13 14) (1 3 13 17) (3 4 17 18) (3 11 15 18) 
;;			  (3 10 14 15) (1 4 17 21) (4 5 21 22) (4 7 19 22) 
;;			  (4 11 18 19) (1 5 21 25) (5 6 25 26) (5 8 23 26) 
;;			  (5 7 22 23) (1 6 25 29) (6 2 29 30) (6 9 27 30) 
;;			  (6 8 26 27) (10 11 15 16) (11 12 16 20) (11 7 19 20)
;;			  (7 12 20 24) (7 8 23 24) (8 12 24 28) (8 9 27 28) 
;;			  (9 12 28 32) (9 10 31 32) (2 10 31 14) (2 9 30 31) 
;;			  (1 2 29 13) (10 12 32 16)))
;;
;;(defparameter S2-1/2-13 '((1 33 13 17) (1 33 17 21) (1 33 21 25) (1 33 25 29) 
;;			  (1 33 29 13) (2 34 13 14) (2 34 14 31) (2 34 30 31) 
;;			  (2 34 29 30) (2 34 29 13) (3 35 13 17) (3 35 17 18) 
;;			  (3 35 18 15) (3 35 15 14) (3 35 14 13) (4 36 17 21) 
;;			  (4 36 21 22) (4 36 22 19) (4 36 19 18) (4 36 18 17)
;;			  (5 37 21 25) (5 37 25 26) (5 37 26 23) (5 37 23 22) 
;;			  (5 37 22 21) (6 38 25 29) (6 38 29 30) (6 38 30 27) 
;;			  (6 38 27 26) (6 38 26 25) (7 39 19 22) (7 39 22 23) 
;;			  (7 39 23 24) (7 39 24 20) (7 39 20 19) (8 40 23 26) 
;;			  (8 40 26 27) (8 40 27 28) (8 40 28 24) (8 40 24 23) 
;;			  (9 41 27 30) (9 41 30 31) (9 41 31 32) (9 41 32 28) 
;;			  (9 41 28 27) (10 42 31 14) (10 42 14 15) (10 42 15 16) 
;;			  (10 42 16 32) (10 42 32 31) (11 43 15 18) (11 43 18 19) 
;;			  (11 43 19 20) (11 43 20 16) (11 43 16 15) (12 44 16 20) 
;;			  (12 44 20 24) (12 44 24 28) (12 44 28 32) (12 44 32 16)))


(defun reset-spacetime ()
  "the state of the simulation after (reset-spacetime) is identical 
to the state after (load \"cdt2p1.lisp\")"
  ;; clear the hash tables
  (clrhash *2SIMPLEX->ID*)
  (clrhash *ID->2SIMPLEX*)
  (clrhash *ID->3SIMPLEX*)
  ;; reset the counters
  (setf *LAST-USED-2SXID* 0)
  (setf *LAST-USED-3SXID* 0)
  (setf *RECYCLED-3SX-IDS* '())
  (setf *LAST-USED-POINT* 0)
  ;; reset the ''bulk'' variables
  (setf N0 0)
  (setf N1-SL 0)
  (setf N1-TL 0)
  (setf N2-SL 0)
  (setf N2-TL 0)
  (setf N3-TL-31 0)
  (setf N3-TL-22 0)
  ;; reset the parameters
  (setf k0 0.0)
  (setf k3 0.0)
  (setf eps 0.02)
  (setf *a* 1.0)
  (setf *alpha* -1.0)
  (setf *k* 1.0)
  (setf *litL* 1.0)
  (setf CURRENT-MOVE-IDENTIFIER "UNKNOWN")
  (setf CURRENT-MOVE-NUMBER 0)
  (setf STOPOLOGY "unknown")
  (setf BCTYPE "unknown")
  (setf SAVE-EVERY-N-SWEEPS 10)
  (setf NUM-T 666666)
  (setf N-INIT 0)
  (setf NUM-SWEEPS 0)
  ;; reset the move markers
  (setf 26MARKER 0.0)
  (setf 23MARKER 0.0)
  (setf 44MARKER 0.0)
  (setf 32MARKER 0.0)
  (setf 62MARKER 0.0))