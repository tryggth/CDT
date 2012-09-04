;; cdt-2+1-globals.lisp --- all the parameters that might need to be accessed 
;; from multiple files

;; Authors:
;; ------- Rajesh Kommu
;; ------- David Kamensky
;; ------- Jonah Miller (jonah.maxwell.miller@gmail.com)

(setf *random-state* (make-random-state t))

;;; We must re-use id numbers for various simplices to avoid the
;;; numbers getting ridiculously large.
(defparameter *LAST-USED-3SXID* 0)
(defparameter *RECYCLED-3SX-IDS* '())
(defparameter *LAST-USED-POINT* 0)

(defmacro next-pt ()
  `(incf *LAST-USED-POINT*))
(defmacro set-last-used-pt (pt)
  `(setf *LAST-USED-POINT* ,pt))
(defmacro next-3simplex-id ()
  `(if (null *RECYCLED-3SX-IDS*)
       (incf *LAST-USED-3SXID*)
       (pop *RECYCLED-3SX-IDS*)))
(defmacro recycle-3simplex-id (sxid)
  `(push ,sxid *RECYCLED-3SX-IDS*))



;;---------------------------------------------------------------------------
;; timelike subsimplices have the form (type tmlo (p0 p1 ...))
(defun tlsubsx->id-hashfn (tlsx)
  (sxhash (sort (copy-list (third tlsx)) #'<)))
(defun tlsubsx->id-equality (tlsx1 tlsx2)
  (and (= (first tlsx1) (first tlsx2))
       (= (second tlsx1) (second tlsx2))
       (set-equal? (third tlsx1) (third tlsx2))))
(sb-ext:define-hash-table-test tlsubsx->id-equality tlsubsx->id-hashfn)
(defparameter *TL2SIMPLEX->ID* (make-hash-table :test 'tlsubsx->id-equality))
(defparameter *TL1SIMPLEX->ID* (make-hash-table :test 'tlsubsx->id-equality))
;; spacelike subsimplices have the form (tslice (p0 p1 ...))
(defun slsubsx->id-hashfn (slsx)
  (sxhash (sort (copy-list (second slsx)) #'<)))
(defun slsubsx->id-equality (slsx1 slsx2)
  (and (= (first slsx1) (first slsx2)) 
       (set-equal? (second slsx1) (second slsx2))))
(sb-ext:define-hash-table-test slsubsx->id-equality slsubsx->id-hashfn)
(defparameter *SL2SIMPLEX->ID* (make-hash-table :test 'slsubsx->id-equality))
(defparameter *SL1SIMPLEX->ID* (make-hash-table :test 'slsubsx->id-equality))
;;---------------------------------------------------------------------------
(defparameter *ID->SPATIAL-2SIMPLEX* (make-hash-table))
(defparameter *ID->3SIMPLEX* (make-hash-table :test 'equal))

(defconstant 26MTYPE 0 "move type (2,6)")
(defconstant 23MTYPE 1 "move type (2,3)")
(defconstant 44MTYPE 2 "move type (4,4)")
(defconstant 32MTYPE 3 "move type (3,2)")
(defconstant 62MTYPE 4 "move type (6,2)")

(defparameter ATTEMPTED-MOVES (list 1 1 1 1 1) 
  "number of attempted moves for each move type")
(defparameter SUCCESSFUL-MOVES (list 1 1 1 1 1) 
  "number of successful moves for each move type")

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

(defparameter CURRENT-MOVE-IDENTIFIER "UNKNOWN")
(defparameter CURRENT-MOVE-NUMBER 0)
(defparameter STOPOLOGY "unknown" "spatial slice topology --- S2 or T2")
(defparameter BCTYPE "unknown" "boundary conditions --- PERIODIC or OPEN")
(defparameter *merge-faces* nil "Whether or not a periodic simulation 
                                 cares about the boundary.")
(defparameter SAVE-EVERY-N-SWEEPS 10 "save every 10 sweeps by default")
(defparameter NUM-T 666666 
  "number of time slices --- set to a non-zero value so (mod ts NUM-T) works")
(defparameter N-INIT 0 
  "initial volume of spacetime; 
   we try to keep the volume close to this number")
(defparameter NUM-SWEEPS 0 
  "number of sweeps for which the simulation is to run")
(defparameter SIM-START-TIME (cdt-now-str) 
  "set again inside the generate methods for more accurate value")
(defparameter 3SXEXT ".3sx2p1" 
  "used for storing the parameters and 3simplex information")
(defparameter PRGEXT ".prg2p1" 
  "used for keeping track of the progress of a simulation run")
(defparameter MOVEXT ".mov2p1" 
  "used for storing the movie data information")
(defparameter S2SXEXT ".s2sx2p1" 
  "used for storing the spatial 2-simplex information")

;; pi is a builtin constant
(defparameter *k0* 0.0)   ; Coupling constant. Related to k and lambda
(defparameter *k3* 0.0)   ; Coupling constant. Related to k and lambda
(defparameter *eps* 0.02) ; Damping parameter. For Metropolis algorithm
(defparameter *a* 1.0)    ; Space-like edge-length
;; length-squared of a time-like edge is (* *alpha* *a*). Since
;; *alpha* is negative, the action has been Wick Rotated.
(defparameter *alpha* -1.0)
(defparameter *k* 1.0) ; Coupling constant. Changed at runtime.
;; Small lambda in Ambjorn and Loll. It's (/ *k* G), where G is
;; Newton's constant.
(defparameter *litL* 1.0)       
(defparameter *i* #C(0.0 1.0))   ; complex number i
(defparameter *-i* #C(0.0 -1.0)) ; complex number -i
(defparameter *2/i* (/ 2 *i*))
(defparameter *pi/i* (/ pi *i*))
(defparameter *2pi/i* (* *2/i* pi))
(defparameter *3/i* (/ 3 *i*))
;; JM : It's bad lisp style to use constants without *constant*. 
;; TODO: fix these constant names.
(defparameter ROOT2 (sqrt 2.0))
;; Real part of dihedral angle for a 3-simplex around a time-like bone
;; for both types of simplices, magically. Assumes (= *alpha* -1).
;; JM : This is somewhat deprecated. Kept only to prevent breakage.
;; For clarity, I use the formulas in Ambjorn and Loll.
(defparameter KAPPA (/ (acos (/ 1 3)) pi)) 
(defparameter 6ROOT2 (* 6.0 ROOT2))
;; Other useful dihedral angle.
(defparameter 3KAPPAMINUS1 (- (* 3 KAPPA) 1)) 

;;; wrsqrt is the "wick rotated" sqrt function. Basically 
;;; wrsqrt(x) = -i*sqrt(-x) when x < 0 and not i*sqrt(-x). 
;;; So wrsqrt(-1) = -i
(defmacro wrsqrt (val)
  `(if (< ,val 0)
       (* -1 ,*i* (sqrt (* -1 ,val)))
       (sqrt ,val)))

;; Declare we'll use the action later.
(defvar action)

(defvar 26MARKER 0.0)
(defvar 23MARKER 0.0)
(defvar 44MARKER 0.0)
(defvar 32MARKER 0.0)
(defvar 62MARKER 0.0)

;; The damping function is used in accept-move?. It effectively widens
;; the critical surface in phase space.
(defun damping (num3)
  (* *eps* (abs (- num3 N-INIT))))

(defun initialize-move-markers ()
  (setf 26MARKER 5.0)
  (setf 23MARKER (+ 26MARKER 5.0))
  (setf 44MARKER (+ 23MARKER 5.0))
  (setf 32MARKER (+ 44MARKER 5.0))
  (setf 62MARKER (+ 32MARKER 5.0)))

(defun update-move-markers ()
  (let ((num-successful-moves (apply #'+ SUCCESSFUL-MOVES)))
    (setf 26MARKER (float (/ num-successful-moves 
			     (nth 26MTYPE SUCCESSFUL-MOVES))))
    (setf 23MARKER (+ 26MARKER(float (/ num-successful-moves 
					(nth 23MTYPE SUCCESSFUL-MOVES)))))
    (setf 44MARKER (+ 23MARKER(float (/ num-successful-moves 
					(nth 44MTYPE SUCCESSFUL-MOVES)))))
    (setf 32MARKER (+ 44MARKER(float (/ num-successful-moves 
					(nth 32MTYPE SUCCESSFUL-MOVES)))))
    (setf 62MARKER (+ 32MARKER(float (/ num-successful-moves 
					(nth 62MTYPE SUCCESSFUL-MOVES)))))))


;; JM: I've generated a select-move function for debugging. Name the
;; one you want to use "select-move" and name the other one something
;; else.
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

;; The debugging version
;; Annoying. Use at your own risk.
(defun select-move-debug () 
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
    (let ((movename (case mtype 
		       (0 "2->6")
		       (1 "2->3")
		       (2 "4->4")
		       (3 "3->2")
		       (4 "6->2")
		       (otherwise "dunno."))))
       (format t "~%Attempted move: ~a.~%" movename)
       mtype)))

;; For printing a move type in a readable format:
(defun printmove (mtype)
  (let ((movename (case mtype 
		    (0 "2->6")
		    (1 "2->3")
		    (2 "4->4")
		    (3 "3->2")
		    (4 "6->2")
		    (otherwise "dunno."))))
    (format t "~%Accepted move: ~a.~%" movename)))
    


;;;------------------------------------------------------------------------  
;;; initialization data 
;;;------------------------------------------------------------------------  


;;; JM: I am unconvinced the initialization data should be in the
;;; globals section... it seems like it should live in
;;; "initialization.lisp." I guess I'll leave it here for now.

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
; here (3,1) means just (3,1), not (3,1)+(1,3)
(defparameter N3-TL-31-PER-SLICE 4) 
(defparameter S2-1/2-31 '((1 2 3 5) (2 3 4 6) (3 4 1 7) (4 1 2 8)))
(defparameter S2-1/2-22 '((1 2 5 8) (2 3 5 6) (3 1 5 7) (3 4 6 7) 
			  (4 2 6 8) (4 1 7 8)))
(defparameter S2-1/2-13 '((1 5 7 8) (2 5 6 8) (3 5 6 7) (4 6 7 8)))

;; JM: Style rules broken here because I can't figure out how to fix
;; this comment without ruining it.

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

