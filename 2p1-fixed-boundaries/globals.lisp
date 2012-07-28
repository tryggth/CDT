;; cdt-2+1-globals.lisp --- all the parameters that might need to be accessed 
;; from multiple files

(setf *random-state* (make-random-state t))

(defparameter *LAST-USED-3SXID* 0)
(defparameter *RECYCLED-3SX-IDS* '())
(defparameter *LAST-USED-POINT* 0)
(defparameter *LAST-USED-S2SXID* 0)

(defmacro next-pt ()
  `(incf *LAST-USED-POINT*))
(defmacro set-last-used-pt (pt)
  `(setf *LAST-USED-POINT* ,pt))
(defmacro next-s2simplex-id ()
  `(incf *LAST-USED-S2SXID*))
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


;;; Book-keeping for the action terms. We input the f-vector and
;;; b-vector into the action to determine whehter or not to accept a
;;; move. It is worth noting that we ALSO only accept a move depending
;;; on whether or not it changes the boundary. Changes to the boundary
;;; are not allowed.

;; Some macros to determine whether or not a simplex is on a boundary
;;first, define macros to determine helpful things about
;;the position of a particular simplex
(defmacro in-upper-sandwich (sxid) 
  `(and (string= BCTYPE "OPEN") (= (3sx-tmhi (get-3simplex ,sxid)) NUM-T)))
(defmacro in-lower-sandwich (sxid)
  `(and (string= BCTYPE "OPEN") (= (3sx-tmlo (get-3simplex ,sxid)) 0)))
(defmacro in-either-boundary-sandwich (sxid)
  `(or (in-upper-sandwich ,sxid) (in-lower-sandwich ,sxid)))
(defmacro has-face-on-boundary (sxid)
  `(let* ((sx (get-3simplex ,sxid))  ;for fixed boundaries, this is
				     ;unnecessary
	  (ty (3sx-type sx))
	  (th (3sx-tmhi sx))
	  (tl (3sx-tmlo sx)))
     (and (string= BCTYPE "OPEN")
	  (or (and (= ty 1)
		   (or (= th 0) (= th NUM-T)))
	      (and (= ty 3)
		   (or (= tl 0) (= tl NUM-T)))))))

;; The f-vector manipulation functions
(defun set-f-vector (v1 v2 v3 v4 v5 v6 v7)
  (setf N0 v1 N1-SL v2 N1-TL v3 N2-SL v4 N2-TL v5 N3-TL-31 v6 N3-TL-22 v7))
(defun update-f-vector (dv)
  (incf N0       (nth 0 dv))
  (incf N1-SL    (nth 1 dv))
  (incf N1-TL    (nth 2 dv))
  (incf N2-SL    (nth 3 dv))
  (incf N2-TL    (nth 4 dv))
  (incf N3-TL-31 (nth 5 dv))
  (incf N3-TL-22 (nth 6 dv)))

;; Changes in the f-vector depending on moves.
(defparameter DF26 '(1 3 2 2 6 4 0))
(defparameter DF62 '(-1 -3 -2 -2 -6 -4 0))
(defparameter DF44 '(0 0 0 0 0 0 0))
(defparameter DF23 '(0 0 1 0 2 0 1))
(defparameter DF32 '(0 0 -1 0 -2 0 -1))

(defun f-vector ()
  "Print the f-vector."
  (format nil "[~A ~A ~A ~A ~A ~A ~A] : ~A ~%"
	  N0
	  N1-SL
	  N1-TL
	  N2-SL
	  N2-TL
	  N3-TL-31
	  N3-TL-22
	  (N3)))


;;; In addition to the f-vector, which keeps track of total information
;;; (not bulk information), we also have the b-vector, which keeps
;;; track of boundary information only. This is useful for
;;; fixed-boundary systems. Below is all the information the system
;;; needs to keep track of the boundary (not including hash tables).

;;; The b-vector keeps track of the top boundary (at t=t_final) and
;;; the bottom boundary (at t=t_initial) seperately. This is because
;;; these terms are accounted for in the action with opposite sign.

;;; Note: the following definitions are not in the order they appear
;;; in the b-vector!

;; Define some quantities of geometrical objects on the boundaries
(defparameter *N1-SL-TOP* 0) ; Number of spacelike links on the
			     ; boundary. These are the bones used to
			     ; calculate the extrinsic curvature on
			     ; the boundary. For the top boundary.
(defparameter *N1-SL-BOT* 0) ; For the bottom boundary

(defparameter *N3-22-TOP* 0) ; Number of (2,2)-simplexes that have 2
			     ; vertexes on the boundary. Useful for
			     ; extrinsic curvature calculations. For
			     ; the top boundary.
(defparameter *N3-22-BOT* 0) ; For the bottom boundary.

(defparameter *N3-31-TOP* 0) ; Number of (3,1)- or (1,3)-simplexes
			     ; that have at least one vertex on the
			     ; boundary. We don't need them for
			     ; extrinsic curvature calculations, but
			     ; we need to compare them to N3-31
			     ; (total) to figure out the number of
			     ; (3,1)-simplices in the bulk. For the
			     ; top boundary.
(defparameter *N3-31-BOT* 0) ; For the bottom boundary.

;; The above three quantities compose the "b-vector," which keeps
;; track of boundary information. The following two functions
;; manipulate the b-vector.

;; Extract the order in the b-vector from these definitions.
(defun set-b-vector (n1 n2 n3 n4 n5 n6)
  "Input: *N1-SL-TOP*, *N1-SL-BOT*, *N3-22-TOP*, *N3-22-BOT*, 
*N3-31-TOP*, *N3-31-BOT*."
  (setf *N1-SL-TOP* n1
	*N3-22-TOP* n2
	*N3-31-TOP* n3
	*N1-SL-BOT* n4
	*N3-22-BOT* n5
	*N3-31-BOT* n6))
(defun update-b-vector (dv)
  "Input list: (DELTA-N1-SL-TOP DELTA-N1-SL-BOT
 DELTA-N3-22-TOP DELTA-N3-22-BOT DELTA-N3-31-TOP DELTA-N3-31-BOT)"
  (incf *N1-SL-TOP* (first  dv))
  (incf *N3-22-TOP* (second dv))
  (incf *N3-31-TOP* (third  dv))
  (incf *N1-SL-BOT* (fourth dv))
  (incf *N3-22-BOT* (fifth  dv))
  (incf *N3-31-BOT* (sixth  dv)))

;; As with the f-vector, the change in the b-vector depends on the
;; move that acts on the spacetime. Unlike the f-vector, the change in
;; the b-vector ALSO depends on where the move occurs. If the move
;; acts on a 2-simplex that is NOT on the the boundary, the b-vector
;; does not change. The only move allowed to act on the boundary is
;; the 23-move and its inverse, the 32-move. Thus these are the only
;; moves which have a chance of changing the b-vector.
(defmacro DB23 (sxid) ; Change in b-vector due to 23-move.
  `(cond ((in-upper-sandwich ,sxid) (list 0 1 0 0 0 0))
	 ((in-lower-sandwich ,sxid) (list 0 0 0 0 1 0))
	 (t                         (list 0 0 0 0 0 0))))

;; change in b-vector due to 32-move. Inverse of DB23
(defmacro DB32 (sxid) 
  `(cond ((in-upper-sandwich ,sxid) (list 0 -1 0 0 0  0))
	 ((in-lower-sandwich ,sxid) (list 0  0 0 0 -1 0))
	 (t                         (list 0  0 0 0  0 0))))

;; The following moves can't affect the boundary
(defparameter DB44 '(0 0 0 0 0 0)) ; Change in b-vector due to a
				   ; 44-move. 44 is its own inverse.
(defparameter DB26 '(0 0 0 0 0 0)) ; Change in b-vector due to a
				   ; 26-move.
(defparameter DB62 '(0 0 0 0 0 0)) ; Change in b-vector due to a
				   ; 62-move. Would be the inverse of
				   ; DB26 if there were any change.

(defun b-vector()
  "Print the b-vector."
  (format t "[~A, ~A, ~A, ~A, ~A, ~A]~%"
	  *N1-SL-TOP*
	  *N3-22-TOP*
	  *N3-31-TOP*
	  *N1-SL-BOT*
	  *N3-22-BOT*
	  *N3-31-BOT*))

(defparameter CURRENT-MOVE-IDENTIFIER "UNKNOWN")
(defparameter CURRENT-MOVE-NUMBER 0)
(defparameter STOPOLOGY "unknown" "spatial slice topology --- S2 or T2")
(defparameter BCTYPE "unknown" "boundary conditions --- PERIODIC or OPEN")
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

;; STOPOLOGY-BCTYPE-NUMT-NINIT-k0-k3-eps-alpha-startsweep-endsweep-hostname-currenttime
(defun generate-filename (&optional (start-sweep 1) 
			    (end-sweep (+ start-sweep NUM-SWEEPS -1)))
  (format nil "~A-~A-T~3,'0d-V~6,'0d-~A-~A-~A-~A-~9,'0d-~9,'0d-on-~A-started~A" 
	  STOPOLOGY BCTYPE NUM-T N-INIT *k0* *k3* *eps* *alpha* start-sweep 
	  end-sweep (hostname) (cdt-now-str)))

;; STOPOLOGY-BCTYPE-NUMT-NINIT-k0-k3-eps-alpha-startsweep-currsweep-endsweep-hostname-starttime-currenttime
(defun generate-filename-v2 (&optional (ssweep 1) (csweep 0) 
			       (esweep (+ ssweep NUM-SWEEPS -1)))
  (format nil "~A-~A-T~3,'0d-V~6,'0d-~A-~A-~A-~A-~9,'0d-~9,'0d-~9,'0d-on-~A-start~A-curr~A" 
	  STOPOLOGY BCTYPE NUM-T N-INIT *k0* *k3* *eps* *alpha* 
	  ssweep csweep esweep 
	  (hostname) SIM-START-TIME (cdt-now-str)))

(defvar 26MARKER 0.0)
(defvar 23MARKER 0.0)
(defvar 44MARKER 0.0)
(defvar 32MARKER 0.0)
(defvar 62MARKER 0.0)

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
    


;;; The following functions will be used to construct the
;;; action. action-exposed is the form of the action, and this is
;;; where changes to the form of the action should be made. It is not,
;;; however, the action function used in the metropolis
;;; algorithm. Instead, make-action is called by set-k0-k3-alpha (or
;;; alternatively set-k-litL-alpha), which then sets the variable
;;; name, 'action, which was set with defvar above to the proper
;;; function. This saves computational time because alpha, k, and litL
;;; don't have to passed to the function every time, but are still
;;; change-able at run-time. This is a trick that would only work in a
;;; language like lisp, where functions are compiled but an
;;; interpreted environment still exists at runtime.

;; Functional form of the corrected action that uses arbitrary alpha
;; and arbitrary k and lambda. Set ret-coup to true for
;; debugging. Note that the action is purely complex. This is expected
;; after Wick rotation and the correct partition function is e^{i
;; action}. Note that if the b-vector is zero (i.e., there is no
;; boundary, i.e., we have periodic boundary conditions), the action
;; reduces to the bulk action.
(defun action-exposed (num1-sl num1-tl num3-31 num3-22 
	       ; Some additional arguments that need to be passed for
	       ; open boundary conditions. In theory, the number of
	       ; (3,1)-simplices connected to the boundary should be
	       ; looked at too. At the fixed bounary, they don't
	       ; change and have no effect on the dynamics. However at
	       ; the other boundary, which we do not keep fixed, they
	       ; have a substantial effect.
	      num1-sl-top ; spacelike links on top boundary
	      num3-22-top ; (2,2)-simplices connected to top boundary
	      num3-31-top ; (3,1)&(1,3)-simplices connected to top boundary
	      num1-sl-bot ; spacelike links on bottom boundary
	      num3-22-bot ; (2,2)-simplices connected to bottom boundary
	      num3-31-bot ; (3,1)&(1,3)-simplices connected to bottom boundary
	      alpha k litL ; Tuning parameters
	      &optional (ret-coup nil))
  "The action before setting coupling constants."
  (let* ((2alpha+1 (+ (* 2 alpha) 1)) ; self-explanatory
	 (4alpha+1 (+ (* 4 alpha) 1))
	 (4alpha+2 (+ (* 4 alpha) 2))
	 (3alpha+1 (+ (* 3 alpha) 1))
	 ; dihedral angle around spacelike bone for (2,2) simplices
	 (theta-22-sl (asin (/ (* *-i* (wrsqrt (* 8 2alpha+1))) 4alpha+1))) 
	 ; dihedral angle around spacelike bones for (3,1) simplices
	 (theta-31-sl (acos (/ *-i* (wrsqrt (* 3 4alpha+1)))))
	 ;dihedral angle around timelike bones for (2,2) simplices
	 (theta-22-tl (acos (/ -1 4alpha+1)))
	 ; dihedral angle around time-like bones for (3,1) simplices
	 (theta-31-tl (acos (/ 2alpha+1 4alpha+1)))
	 ; 3-volume of a (2,2)-simplex
	 (v3-22 (wrsqrt 4alpha+2))
	 ; 3-volume of a (3,1)- or (1,3)-simplex
	 (v3-31 (wrsqrt 3alpha+1))

	 ;; Coefficients for action assuming closed manifold.

	 ;; BULK

         ; spacelike edges term (for total angle we take the deficit from)
	 (K1SL (* *pi/i* k)) 
	  ; timelike edges term (for the total angle we take deficit from)
	 (K1TL (* (wrsqrt alpha) 2 pi k))  
	 ; coefficient for (2,2)-simplices around spacelike bone term
	 (K3SL-22 (* -1  (/ k *i*) theta-22-sl))  
	 ; coefficient for dihedral angle of (3,1)-simplices at each edge.
	 (K3SL-31 (* -1 3 (/ k *i*) theta-31-sl))  
	 ; term for dihedral angle around timelike edges of (2,2)-simplices
	 (K3TL-22 (* -4 (wrsqrt alpha) k theta-22-tl)) 
	 ; term for dihedral angle around timelike edges of (1,3)- and
	 ; (3,1)-simplices
	 (K3TL-31 (* -3 k (wrsqrt alpha) theta-31-tl)) 
	 (KV22 (* -1 (/ litL 12) v3-22)) ; volume term for (2,2)-simplices
	 (KV31 (* -1 (/ litL 12) v3-31)) ; volume term for (3,1)-simplices

	 ;; BOUNDARY

	 (B1SL (* k (/ pi *i*))) ; term for sum over spacelike edges
	 ; term for sum over dihedral angles of (3,1)-simplices at each bone.
	 (B3SL-31 (* -1 *2/i* k theta-31-sl)) 
	 ; term for sum over dihedral angles of (2,2)-simplices
	 ; attached at each bone
	 (B3SL-22 (* -1 (/ k *i*) theta-22-sl))) 

    ;; JM: Note that I keep track of both boundaries separately, even
    ;; though this isn't really necessary. This is for clarity and in
    ;; case my conception of the action is incorrect, changes are
    ;; easier.	     

    (if ret-coup ; if ret-coup, just return debugging info. Otherwise,
		 ; return the action.
;;	(list (* 2 K1SL) K1TL
;;	      (+ K3SL-31 K3TL-31 KV31) (+ (* 2 K3SL-22) K3TL-22 KV22)
;;	      K3SL-31 K3TL-31 KV31)
;;	(list theta-22-sl theta-31-sl theta-22-tl theta-31-tl
;;	      v3-22 v3-31 K1SL K1TL K3SL-22 K3TL-22 K3TL-31 KV22 KV31
;;	      B1SL B3SL-31 B3SL-22)
	(list  (+ (* 2 K3SL-22) K3TL-22 KV22) (* 2 K3SL-22) K3TL-22 KV22)

	(+ ;; BULK TERM
	 ; we need to subtract the boundary simplices from the bulk
	 (* K1SL (- (* 2 num1-sl) (+ num1-sl-top num1-sl-bot)))
	 (* K1TL num1-tl) ; there are no boundary timelike simplices 

	 ; We subtract half of the (2,2)-simplices at the boundary
	 ; because each (2,2)-simplex at the boundary contributes one
	 ; dihedral angle to the bulk sum and one dihedral angle to
	 ; the boundary sum, while (2,2)-simplices in the bulk
	 ; contribute 2 dihedral angles to the bulk sum.
	 (* K3SL-22 (- (* 2 num3-22) (+ num3-22-top num3-22-bot)))
	 ; For (3,1)-simplices attacked to the boundary, we have the
	 ; same problem as with (2,2)-simplices.
	 (* K3SL-31 (- num3-31 (+ num3-31-top num3-31-bot)))
	 ; There are no timelike bones in the boundary, so we don't
	 ; have to subtract for angles around timelike bones.
	 (* K3TL-31 num3-31) 
	 (* K3TL-22 num3-22)
	 ; volume terms
	 (* KV22 num3-22)
	 (* KV31 num3-31)

	 ;; BOUNDARY TERM
	 ;; top boundary
	 (* B1SL num1-sl-top) ; pi term for deficit angle
	 (* B3SL-31 num1-sl-top) ; (3,1)-simplex contribution to deficit angle
	 (* B3SL-22 num3-22-top) ; (2,2)-simplex contribution to deficit angle
	 ;; bottom boundary
	 (* B1SL num1-sl-bot) ; pi term for deficit angle
	 (* B3SL-31 num1-sl-bot) ; (3,1)-simplex contribution to deficit angle
	 (* B3SL-22 num3-22-bot))))) ; (2,2)-simplex contribution to
		                     ; deficit angle
	  
(defun make-action (alpha k litL)
  "Construct an action with fixed coupling constants 
for use in the simulation."
  (setf (symbol-function 'action)
	#'(lambda (num1-sl num1-tl num3-31 num3-22
		   num1-sl-top num3-22-top num3-31-top
		   num1-sl-bot num3-22-bot num3-31-bot)
	    (action-exposed num1-sl num1-tl num3-31 num3-22
			    num1-sl-top num3-22-top num3-31-top
			    num1-sl-bot num3-22-bot num3-31-bot 
			    alpha k litL))))

;;; Functions to set coupling constants and build action. Both
;;; functions are equivalent. They just take different inputs.

;; Takes k0-k3-alpha as input, sets k, lambda, alpha, k0, k3, and then
;; constructs an action and initializes moves.
(defun set-k0-k3-alpha (k0 k3 alpha)
  (setf *k0* k0 *k3* k3 *alpha* alpha)
  (setf *k* (/ *k0* (* 2 *a* pi)))
  (setf *litL* (* (- *k3* (* 2 *a* pi *k* 3KAPPAMINUS1)) 
		  (/ 6ROOT2 (* *a* *a* *a*))))
  (make-action *alpha* *k* *litL*)
  (initialize-move-markers))

;; Function analogous to set-k0-k3-alpha, but takes k, litL, and alpha
;; as inputs.
(defun set-k-litL-alpha (k litL alpha)
  (prog nil
     (setf *k* k
	   *litL*  litL
	   *alpha* alpha)
     (setf *k0* (* *k* (* 2 *a* pi))
	   *k3* (+ (/ (* *litL* *a* *a* *a*) 6ROOT2) 
		   (* 2 *a* pi *k* 3KAPPAMINUS1)))
     (make-action *alpha* *k* *litL*)
     (initialize-move-markers)))

;; JM: Deprecated. I've replaced this with the functions
;; above. However, for completeness and for debugging, I've left the
;; original function in here. A warning: IT DOES NOT CONSTRUCT AN
;; ACTION WITH THE BOUNDARY CONDITION TERM INCLUDED.
(defun set-k0-k3-alpha-deprecated (kay0 kay3 alpha)
  (setf *k0* kay0 *k3* kay3 *alpha* alpha)
  (initialize-move-markers)
  (let* ((k (/ *k0* (* 2 *a* pi)))
	 (litL (* (- *k3* (* 2 *a* pi k 3KAPPAMINUS1)) 
		  (/ 6ROOT2 (* *a* *a* *a*))))
	 (2alpha+1 (+ (* 2 *alpha*) 1))
	 (4alpha+1 (+ (* 4 *alpha*) 1))
	 (4alpha+2 (+ (* 4 *alpha*) 2))
	 (3alpha+1 (+ (* 3 *alpha*) 1))
	 (arcsin-1 (asin (/ (* *-i* (wrsqrt (* 8 2alpha+1))) 4alpha+1)))
	 (arccos-1 (acos (/ *-i* (wrsqrt (* 3 4alpha+1)))))
	 (arccos-2 (acos (/ -1 4alpha+1)))
	 (arccos-3 (acos (/ 2alpha+1 4alpha+1)))
	 (k1SL (* *2pi/i* k))
	 (k1TL (* 2 pi k (wrsqrt *alpha*)))
	 (k3TL31 (+ (* k *3/i* arccos-1) 
		    (* 3 k (wrsqrt *alpha*) arccos-3)
		    (* (/ litL 12) (wrsqrt 3alpha+1))))
	 (k3TL22 (+ (* k *2/i* arcsin-1)
		    (* 4 k (wrsqrt *alpha*) arccos-2)
		    (* (/ litL 12) (wrsqrt 4alpha+2)))))
    (setf (symbol-function 'action)
	  #'(lambda (n1SL n1TL n3TL31 n3TL22)
	      (- (+ (* k1SL n1SL) (* k1TL n1TL)) 
		 (+ (* k3TL31 n3TL31) (* k3TL22 n3TL22)))))))
  

;;;------------------------------------------------------------------------  
;;; initialization data 
;;;------------------------------------------------------------------------  


;;; JM: I am unconvinced the initialization data should be in the
;;; globals section... it seems like it should live in
;;; "initialization.lisp." I guess I'll leave it here for now.

;;; CHANGELOG: Using David Kemansky's initialization data rather than
;;; Rajesh Kommu's original initialization data.



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



;;;--------------------------------------------------------------------------
;;; RESET THE SIMULATION
;;;--------------------------------------------------------------------------

;;; JM: Resetting the simulation at runtime is tricky business. The
;;; following two functions reset the system in different ways. It's
;;; not terribly important that these functions work properly. They're
;;; mostly usefull for debuging at runtime.

;; reset-spacetime-fast is a faster way to reset the spacetime, but it
;; is dangerous because I am not certain all hash tables and important
;; variables are reset. Use at your own risk.
(defun reset-spacetime-fast ()
  "the state of the simulation after (reset-spacetime) is identical 
to the state after (load \"cdt2p1.lisp\"). Use at your own risk."
  ;; clear the hash tables
  (clrhash *TL2SIMPLEX->ID*)
  (clrhash *SL2SIMPLEX->ID*)
  (clrhash *TL1SIMPLEX->ID*)
  (clrhash *SL1SIMPLEX->ID*)
  (clrhash *ID->SPATIAL-2SIMPLEX*)
  (clrhash *ID->3SIMPLEX*)
  ;; reset the counters
  ;; (setf *LAST-USED-2SXID* 0)
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
  ;; Reset the ''boundary'' variables
  (setf *N1-SL-TOP* 0)
  (setf *N1-SL-BOT* 0)
  (setf *N3-22-TOP* 0)
  (setf *N3-22-BOT* 0)
  (setf *N3-31-TOP* 0)
  (setf *N3-31-BOT* 0)
  ;; reset the parameters
  (setf *k0* 0.0)
  (setf *k3* 0.0)
  (setf *eps* 0.02)
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
  (setf 62MARKER 0.0)
  ;; Reset attempted moves list
  (reset-move-counts))

;; reset-spacetime-slow simply reloads all modules and essentially
;; restarts the simulation. Only useful for runtime debugging.
(defun reset-spacetime-slow nil
  "Reloads all modules and restarts simulation. Use at your own risk."
  (load "cdt2p1.lisp"))

;; For debugging. Use only to quickly test the simulation while in the
;; lisp REPL.
;(defun test-initialization nil
;  "A simple initialization command for use in debugging."
;  (initialize-t-slices-with-v-volume :num-time-slices 64
;				     :target-volume   80000
;				     :spatial-topology "s2"
;				     :boundary-conditions "open"
;				     :initial-spatial-geometry "tetra.txt"
;				     :final-spatial-geometry "tetra.txt"))
