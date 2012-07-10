(load "../utilities.lisp")
(load "globals.lisp")
(load "simplex.lisp")
(load "moves.lisp")
(load "initialization.lisp")
(load "montecarlo.lisp")
;; cdt-2+1-globals.lisp --- all the parameters that might need to be accessed 
;; from multiple files

(setf *random-state* (make-random-state t))

;; IDs for the different simplices
;; CA : a 2sx or 3sx could just be hashed
(defparameter *LAST-USED-2SXID* 0)
(defparameter *LAST-USED-3SXID* 0)
(defparameter *RECYCLED-3SX-IDS* '())

(defparameter *LAST-USED-POINT* 0)
;; CA is there any problem with recycling points?

(defparameter *LAST-USED-S2SXID* 0)

(defmacro next-pt ()
  `(incf *LAST-USED-POINT*))
(defmacro set-last-used-pt (pt)
  `(setf *LAST-USED-POINT* ,pt))
(defmacro next-s2simplex-id ()
  `(incf *LAST-USED-S2SXID*))
(defmacro next-2simplex-id ()
  `(incf *LAST-USED-2SXID*))
(defmacro next-3simplex-id ()
  `(if (null *RECYCLED-3SX-IDS*)
       (incf *LAST-USED-3SXID*)
       (pop *RECYCLED-3SX-IDS*)))
(defmacro recycle-3simplex-id (sxid)
  `(push ,sxid *RECYCLED-3SX-IDS*))

;;
;; Hash tables
;;

;; CA we should guestimate a size for these

;; New in HL:

;; Point -> Triangles that contain it
(defparameter *0SIMPLEX->NS2SX* (make-hash-table))
;;Point -> Its time
(defparameter *0SIMPLEX->TIME* (make-hash-table))

;; Line -> 3simplices that contain it
(defun 1simplex->id-equality (1sx1 1sx2)
  (set-equal? 1sx1 1sx2))
(defun 1simplex->id-hashfn (1sx)
  (sxhash (sort (copy-list 1sx) #'<)))

(sb-ext:define-hash-table-test 1simplex->id-equality 1simplex->id-hashfn)

(defparameter *1SIMPLEX->N3SX22UD* (make-hash-table :test '1simplex->id-equality))
(defparameter *1SIMPLEX->POINTS* (make-hash-table :test '1simplex->id-equality))

;; A related table for debugging only

(defparameter *1SIMPLEX?* (make-hash-table :test '1simplex->id-equality))

;; CA these were already in place and are pretty bad ideas
;; IDs are probably ~sequential, so why use hash table for ID -> foo
;; In the case of 2simplices, why should they have ids?

(defun 2simplex->id-equality (2sx1 2sx2)
  (set-equal? (fourth 2sx1) (fourth 2sx2)))
(defun 2simplex->id-hashfn (2sx)
  (sxhash (sort (copy-list (fourth 2sx)) #'<)))

(sb-ext:define-hash-table-test 2simplex->id-equality 2simplex->id-hashfn)

(defparameter *2SIMPLEX->ID* (make-hash-table :test '2simplex->id-equality)) 
(defparameter *ID->2SIMPLEX* (make-hash-table))
(defparameter *ID->SPATIAL-2SIMPLEX* (make-hash-table))
(defparameter *ID->3SIMPLEX* (make-hash-table :test 'equal))

;; Arrays

;; An array that stores the number of spatial triangles in each time slice
;; Used to:
;; 1) give a quick view of spatial volume during run time
;; 2) Reject 6->2 moves where the spatial slice has only four triangles.

(defvar *SVOLS*)

;; The number of points in each time slice. 
(defvar *0SX-PER-TIME*)

(defun reset-global-state ()
  (setf *LAST-USED-2SXID* 0)
  (setf *LAST-USED-3SXID* 0)
  (setf *RECYCLED-3SX-IDS* '())
  (setf *LAST-USED-POINT* 0)
  (setf *LAST-USED-S2SXID* 0)
  (setf *0SIMPLEX->NS2SX* (make-hash-table))
  (setf *1SIMPLEX->N3SX22UD* (make-hash-table :test '1simplex->id-equality))
  (setf *1SIMPLEX->POINTS* (make-hash-table :test '1simplex->id-equality))
  (setf *2SIMPLEX->ID* (make-hash-table :test '2simplex->id-equality)) 
  (setf *ID->2SIMPLEX* (make-hash-table))
  (setf *ID->SPATIAL-2SIMPLEX* (make-hash-table))
  (setf *ID->3SIMPLEX* (make-hash-table :test 'equal))
  (setf *SVOLS* (make-array NUM-T :initial-element 0)))

;; Constants. No more magic numbers!

(defconstant 3SX31 3)
(defconstant 3SX22 2)
(defconstant 3SX13 1)

(defconstant 2SXSPATIAL 0)
(defconstant 2SX21 1)
(defconstant 2SX12 2)

(defconstant 26MTYPE 0 "move type (2,6)")
(defconstant 23MTYPE 1 "move type (2,3)")
(defconstant 44MTYPE 2 "move type (4,4)")
(defconstant 32MTYPE 3 "move type (3,2)")
(defconstant 62MTYPE 4 "move type (6,2)")

(defparameter ATTEMPTED-MOVES (list 1 1 1 1 1) 
  "number of attempted moves for each move type")
(defparameter SUCCESSFUL-MOVES (list 1 1 1 1 1) 
  "number of successful moves for each move type")
(defparameter THRASHING 0)

(defun reset-move-counts ()
  (for (n 0 4)
       (setf (nth n ATTEMPTED-MOVES) 1 (nth n SUCCESSFUL-MOVES) 1)))

(defun accept-ratios ()
  (format nil "~A [~A ~A ~A ~A ~A]"
	  ATTEMPTED-MOVES
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

(defun f-vector ()
  (format nil "[~A ~A ~A ~A ~A ~A ~A] : ~A ~%"
	  N0
	  N1-SL
	  N1-TL
	  N2-SL
	  N2-TL
	  N3-TL-31
	  N3-TL-22
	  (N3)))

;;Update f-vector from delta vector
(defun update-f-vector (dv)
  (incf N0 (nth 0 dv))
  (incf N1-SL (nth 1 dv))
  (incf N1-TL (nth 2 dv))
  (incf N2-SL (nth 3 dv))
  (incf N2-TL (nth 4 dv))
  (incf N3-TL-31 (nth 5 dv))
  (incf N3-TL-22 (nth 6 dv)))

;;Delta vectors for the different move types
(defparameter DF26 '(1 3 2 2 6 4 0))
(defparameter DF62 '(-1 -3 -2 -2 -6 -4 0))
(defparameter DF44 '(0 0 0 0 0 0 0))
(defparameter DF23 '(0 0 1 0 2 0 1))
(defparameter DF32 '(0 0 -1 0 -2 0 -1))

;;Additional run-time state
(defparameter CURRENT-MOVE-IDENTIFIER "UNKNOWN")
(defparameter CURRENT-MOVE-NUMBER 0)

;;Set at start time
(defparameter STOPOLOGY "unknown" "spatial slice topology --- S2 or T2")
(defparameter BCTYPE "unknown" "boundary conditions --- PERIODIC or OPEN")
(defparameter SAVE-EVERY-N-SWEEPS 10 "save every 10 sweeps by default")
(defparameter NUM-T 666666 
  "number of time slices --- set to a non-zero value so (mod ts NUM-T) works")
(defparameter N-INIT 0 
  "initial volume of spacetime; we try to keep the volume close to this number")
(defparameter NUM-SWEEPS 0 
  "number of sweeps for which the simulation is to run")
(defparameter SIM-START-TIME (cdt-now-str) 
  "set again inside the generate methods for more accurate value")
(defparameter TARG-VOL-IDEAL 0
  "useful to know if the program wants to analyze its own behavior")

;;File extensions
(defparameter 3SXEXT ".3sx2p1" 
  "used for storing the parameters and 3simplex information")
(defparameter PRGEXT ".prg2p1" 
  "used for keeping track of the progress of a simulation run")
(defparameter MOVEXT ".mov2p1" 
  "used for storing the movie data information")
(defparameter S2SXEXT ".s2sx2p1" 
  "used for storing the spatial 2-simplex information")

;; wrsqrt is the "wick rotated" sqrt function. 
;; Basically wrsqrt(x) = -i*sqrt(-x) when x < 0 and not 
;; i*sqrt(-x). So wrsqrt(-1) = -i
(defmacro wrsqrt (val)
  `(if (< ,val 0)
       (* -1 *i* (sqrt (* -1 ,val)))
       (sqrt ,val)))

;;Numerical parameters

;; The free physical constants
(defparameter *k0* 0.0 "Per AJL")
(defparameter *k3* 0.0 "Per AJL")
(defparameter *alpha* -1.0 "The length scaling.")
(defparameter *lambda* 0.0 "The K^2 coupling.")
(defparameter *mu* 0.0 "The R^2 coupling.")

;; CA: I think this is the length of spatial edges?
(defparameter *a* 1.0)

;; The damping constant for MC
(defparameter *eps* 0.02)

;; The derived physical constants
(defparameter *k* 1.0)
(defparameter *litL* 1.0)
(defparameter *theta22* 1.0)
(defparameter *theta31* 1.0)
(defparameter *eta* 1.0)

;; Numerical constants
(defparameter *i* #C(0.0 1.0)) ;; complex number i
(defparameter *-i* #C(0.0 -1.0)) ;; complex number -i
(defparameter *2/i* (/ 2 *i*))
(defparameter *2pi/i* (* *2/i* pi))
(defparameter *3/i* (/ 3 *i*))
(defparameter ROOT2 (sqrt 2.0))
(defparameter KAPPA (/ (acos (/ 1 3)) pi))
(defparameter 6ROOT2 (* 6.0 ROOT2))
(defparameter 3KAPPAMINUS1 (- (* 3 KAPPA) 1))

;; STOPOLOGY-BCTYPE-NUMT-NINIT-k0-k3-eps-alpha-startsweep-endsweep-hostname-currenttime
(defun generate-filename (&optional (start-sweep 1) (end-sweep (+ start-sweep NUM-SWEEPS -1)))
  (format nil "~A-~A-T~3,'0d-V~6,'0d-~A-~A-~A-~A-~A-~A-~9,'0d-~9,'0d-on-~A-started~A" 
	  STOPOLOGY BCTYPE NUM-T N-INIT *k0* *k3* *eps* *alpha* *lambda* *mu* start-sweep end-sweep (hostname) (cdt-now-str)))

;; STOPOLOGY-BCTYPE-NUMT-NINIT-k0-k3-eps-alpha-startsweep-currsweep-endsweep-hostname-starttime-currenttime
(defun generate-filename-v2 (&optional (ssweep 1) (csweep 0) (esweep (+ ssweep NUM-SWEEPS -1)))
  (format nil "~A-~A-T~3,'0d-V~6,'0d-~A-~A-~A-~A-~A-~A-~9,'0d-~9,'0d-~9,'0d-on-~A-start~A-curr~A" 
	  STOPOLOGY BCTYPE NUM-T N-INIT *k0* *k3* *eps* *alpha* *lambda* *mu* ssweep csweep esweep 
	  (hostname) SIM-START-TIME (cdt-now-str)))

;; CA: shouldn't this function read in an fvector instead?

;; refer to eqn (35) of dynamically triangulating lorentzian quantum gravity
;;(defun action-n1TL-n3TL31-n3TL22 (n1tl n3tl31 n3tl22)
;;  (+
;;   (* 2 pi invG (sqrt alpha) n1tl)
;;   (* n3tl31 (+ (* -3 invG (asinh (/ 1 (sqrt (+ (* 12 alpha) 3)))))
;;		(* -3 invG (sqrt alpha) (acos (/ (+ (* 2 alpha) 1) (+ (* 4 alpha) 1))))
;;		(* -1 (small-lambda) (sqrt (+ (* 3 alpha) 1)) 1/12)))
;;   (* n3tl22 (+ (* 2 invG (asinh (/ (sqrt (+ (* 16 alpha) 8)) (+ (* 4 alpha) 1))))
;;		(* -4 invG (sqrt alpha) (acos (/ -1 (+ (* 4 alpha) 1))))
;;		(* -1 (small-lambda) (sqrt (+ (* 4 alpha) 2)) 1/12)))
;;   (* -1 *eps* (abs (- (+ n3tl31 n3tl22) N-INIT)))))

;; CA: this code should have a more elaborate system for citing its sources
;; refer to eqn (34) of Dynamically Triangulating Lorentzian Quantum Gravity

(defun damping (num3)
  "A function to generate a damping exponent based on volume."
  (* *eps* (abs (- num3 N-INIT))))

;; When choosing a move we generate (random 62MARKER)
;; and see which of these buckets it falls into 

(defvar 26MARKER 0.0)
(defvar 23MARKER 0.0)
(defvar 44MARKER 0.0)
(defvar 32MARKER 0.0)
(defvar 62MARKER 0.0)

(defun initialize-move-markers ()
  (setf 26MARKER 5.0)
  (setf 23MARKER (+ 26MARKER 5.0))
  (setf 44MARKER (+ 23MARKER 5.0))
  (setf 32MARKER (+ 44MARKER 5.0))
  (setf 62MARKER (+ 32MARKER 5.0)))

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
;;    (format t "Randomly selected move type: ~A ~%" mtype) ;;COMP
    mtype))

;; CA: this is a very clumsy way to do this
;; ... but this function is never called
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
					(nth 62MTYPE
					     SUCCESSFUL-MOVES)))))))

(defun action-exposed (num1-sl num1-tl num3-31 num3-22 alpha k litL &optional (ret-coup nil))
  (let* ((2alpha+1 (+ (* 2 alpha) 1))
	 (4alpha+1 (+ (* 4 alpha) 1))
	 (4alpha+2 (+ (* 4 alpha) 2))
	 (3alpha+1 (+ (* 3 alpha) 1))
	 (arcsin-1 (asin (/ (* *-i* (wrsqrt (* 8 2alpha+1))) 4alpha+1)))
	 (arccos-1 (acos (/ *-i* (wrsqrt (* 3 4alpha+1)))))
	 (arccos-2 (acos (/ -1 4alpha+1)))
	 (arccos-3 (acos (/ 2alpha+1 4alpha+1)))
	 (A (* *2pi/i* k))
	 (B (* (wrsqrt alpha) 2 pi k))
	 (C (- (+ (* *3/i* arccos-1 k) (* (wrsqrt alpha) 3 arccos-3 k) (* (/ litL 12) (wrsqrt 3alpha+1)))))
	 (D (- (+ (* *2/i* arcsin-1 k) (* (wrsqrt alpha) 4 arccos-2 k) (* (/ litL 12) (wrsqrt 4alpha+2))))))
    (if ret-coup
	(list A B C D)
	(+ (* A num1-sl) (* B num1-tl) (* C num3-31) (* D num3-22)))))

(defun build-regge-action (A B C D)
  (format t "Building regge action: ~A ~A ~A ~A ~%" A B C D)
  (defun action-cached (num1-sl num1-tl num3-31 num3-22)
    (let ((out (+ (* A num1-sl) 
		  (* B num1-tl) 
		  (* C num3-31) 
		  (* D num3-22))))
      out)))

(defparameter *FIX-SELECTIONS* nil)

;; alpha = 1.8 (defparameter PATTEMPTED (list .317 .346 .136 .168 .033))
;; alpha = -1.8 (defparameter PATTEMPTED (list .324 .340 .128 .158 .050))
;; alpha = 3.0 (defparameter PATTEMPTED (list .313 .349 .130 .173 .023))
;; alpha = 3.0 no 1/N (defparameter PATTEMPTED (list .336 .294 .188 .111 .070))
(if *FIX-SELECTIONS*
    (defparameter PATTEMPTED (list .2 .2 .2 .2 .2))
    (defparameter PATTEMPTED (list .317 .346 .136 .168 .033)))

(defun update-pattempted ()
  (setf PATTEMPTED (mapcar (lambda (x) (/ x (sum ATTEMPTED-MOVES))) ATTEMPTED-MOVES)))

(defvar *r2-pluggable*)
(defvar *k2-pluggable*)
(defvar *regge-action*)

(flet ((r2-change-26 (initial-counts)
	 (zero-or *mu*
		  (- (funcall *r2-pluggable* (cons 3 (mapcar (lambda (x y) (+ x y)) initial-counts (list 1 1 1))))
		     (funcall *r2-pluggable* initial-counts))))
       (r2-change-62 (initial-counts)
	 (zero-or *mu* 
		  (- (funcall *r2-pluggable* (mapcar (lambda (x y) (+ x y)) initial-counts (list -1 -1 -1)))
		     (funcall *r2-pluggable* (cons 3 initial-counts))))))
  (defun tuned-litL (k)
    (flet ((func (x)
	     (+ (* 4 (elt PATTEMPTED 26MTYPE) 
		   (min 1 (exp (- (realpart (* *i* (action-exposed 3 2 4 0 *alpha* k x)))
				  (r2-change-26 (list 6 6 6))))))
		(* 1 (elt PATTEMPTED 23MTYPE) 
		   (min 1 (exp (realpart (* *i* (action-exposed 0 1 0 1 *alpha* k x))))))
		(* -1 (elt PATTEMPTED 32MTYPE)
		   (min 1 (exp (realpart (* *i* (action-exposed 0 -1 0 -1 *alpha* k x))))))
		(* -4 (elt PATTEMPTED 62MTYPE)
		   (min 1 (exp (- (realpart (* *i* (action-exposed -3 -2 -4 0 *alpha* k x)))
				  (r2-change-62 (list 6 6 6)))))))))
      ;; We shoot in the dark to find bounds for the root.
      ;; The honest truth: for general values of func, I have no idea where the root will be.
      (let* ((one 0)
	     (fone (func one))
	     (desired (if (> fone 0) #'< #'>))
	     (two (if (funcall desired (func 1) fone) 1 -1)))
	(until (funcall desired (func two) 0) (setf two (* two 2)))
	(if (> (func one) 0)
	    (progn
	      (assert (< (func two) 0))
	      (false-position-root #'func one two (expt 10 -5)))
	    (progn
	      (assert (> (func two) 0))
	      (false-position-root #'func two one (expt 10 -5))))))))

(flet ((objective-function (litL)
	 (let ((avg 0))
	   (set-new-litL litL)
	   (format t "Tuning litL to ~A ~%" *litL*)
	   (initialize)
	   (for (ns-therm 0 80)
		(sweep)
		(when (> (N3) (* 2 N-INIT))
		  (return-from objective-function (- N-INIT (N3)))))
	   (for (ns 0 20)
		(setf avg (/ (+ (* avg ns) (N3)) (+ ns 1)))
		(when (> (N3) (* 2 N-INIT))
		  (return-from objective-function (- N-INIT (N3)))))
	   (format t "Average vol is ~A ~%" avg)
	   (- N-INIT avg))))
  (defun runtime-tuning ()
    (format t "The value of litL before tuning : ~A ~%" *litL*)
    (let ((out 0))
      (setf out (false-position-root-near #'objective-function *litL* 200))
      (format t "The value of litL after tuning : ~A ~%" *litL*)
      out)))

(defun set-new-litL (litL)
  (setf *litL* litL)
  (setf *regge-action* (apply 'build-regge-action 
			      (action-exposed 0 0 0 0 *alpha* *k* *litL* t))))

(defun set-k0-k3-alpha-lambda-mu (kay0 kay3 alpha lambda mu &optional tune)
  "Sets *k0* *k3* *alpha* *k* *litL* and initializes move markers"
  (setf *k0* kay0 *k3* kay3 *alpha* alpha *lambda* lambda *mu* mu)

  (setf *eta* (- *alpha*))
  (setf *theta22* (acos (/ (- (* 4 *eta*) 3) (- (* 4 *eta*) 1))))
  (setf *theta31* (- (/ pi 2) (acos (/ (* 2 (sqrt (- (* 3 *eta*) 1))) 
				       (* (sqrt 3) (sqrt (- (* 4 *eta*) 1)))))))

  (initialize-move-markers)

  ;; We require alpha to be such that we can plug it into the Lorentzian 
  ;; action and recover i*Euclidean action.
  (assert (< *alpha* -.5))

  (setf *k* (/ *k0* (* 2 *a* pi)))
  (setf *r2-pluggable* (build-r2-pluggable))
  (setf *k2-pluggable* (build-k2-pluggable))
  (if tune
      (progn
	(format t "Tuning is beginning! ~%")
	(set-new-litL (tuned-litL *k*))
	(runtime-tuning)
	(setf *k3* (+ (/ (* *litL* *a* *a* *a*) 6ROOT2) (* 2 pi *a* *k* 3KAPPAMINUS1)))
	(format t "Tuning has finished. The new value of k3 is: ~A ~%" *k3*))
      (set-new-litL (* (- *k3* (* 2 *a* pi *k* 3KAPPAMINUS1)) 
		       (/ 6ROOT2 (* *a* *a* *a*)))))
  (format t "k: ~A, cc: ~A, mu: ~A, lambda: ~A ~%" *k* *litL* *mu* *lambda*))

;; initialization data

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
(defparameter N3-TL-31-PER-SLICE 4) ; here (3,1) means (3,1), not (3,1)+(1,3)
(defparameter S2-1/2-31 '((1 2 3 5) (2 3 4 6) (3 4 1 7) (4 1 2 8)))
(defparameter S2-1/2-22 '((1 2 5 8) (2 3 5 6) (3 1 5 7) (3 4 6 7) 
			  (4 2 6 8) (4 1 7 8)))
(defparameter S2-1/2-13 '((1 5 7 8) (2 5 6 8) (3 5 6 7) (4 6 7 8)))

;; spatial-2simplex
;; (time (p0 p1 p2) (n0 n1 n2))

;; 3simplex
;; (type tmlo tmhi (t0 t1 t2 t3) (n0 n1 n2 n3))
;; type = 1,2,3 tm[lo|hi] - [lo|hi] spatial time slice
;; tj = id of the 2sx that does not have pj
;; nj = id of the 3sx that does not have pj

;;
;; Code for 2 simplices
;;

;; Get a simplex by its ID.
(defmacro get-2simplex (sxid) `(gethash ,sxid *ID->2SIMPLEX*))

;; Read the specified field of a 2-simplex.
(defmacro 2sx-type (2sx) `(first ,2sx))
(defmacro 2sx-tmlo (2sx) `(second ,2sx))
(defmacro 2sx-tmhi (2sx) `(third ,2sx)) 
(defmacro 2sx-points (2sx) `(fourth ,2sx))

;; CA: we maintain two hash tables. Are both really necessary?
(defun make-2simplex (type tmlo tmhi p0 p1 p2)
  "Returns the id of the 2-simplex with the specified points. If a 2-simplex with points 
p0 p1 p2 already exists, the id of that simplex is returned"
  (let* ((2sx (list type tmlo tmhi (list p0 p1 p2)))
	 (2sxid (gethash 2sx *2SIMPLEX->ID*)))
    (unless 2sxid
;;      (format t "Made one: ~A ~A ~A ~%" p0 p1 p2)
      (setf 2sxid (next-2simplex-id))
      (setf (gethash 2sx *2SIMPLEX->ID*) 2sxid)
      (setf (gethash 2sxid *ID->2SIMPLEX*) 2sx)
      (mapc (lambda (x) (gethash-update x *1SIMPLEX?* '(1)) t) (pairs (list p0 p1 p2)))
      (when (= type 2SXSPATIAL)
	(incf (elt *SVOLS* tmlo))
;;	(mapc (lambda (x) (setf (gethash x *0SIMPLEX->TIME*) (2sx-tmlo 2sx))) (list p0 p1 p2))
	(mapc (lambda (x) (gethash-update x *0SIMPLEX->NS2SX* '(1))) (list p0 p1 p2))
	(mapc (lambda (x) (gethash-push (car x) *1SIMPLEX->POINTS* (cadr x)))
	      (list (list (list p0 p1) p2)
		    (list (list p0 p2) p1)
		    (list (list p1 p2) p0)))))
    2sxid))

(defmacro remove-2simplex (2sxid)
  `(let* ((2sx (gethash ,2sxid *ID->2SIMPLEX*))
	  (points (2sx-points 2sx))
	  (p0 (car points))
	  (p1 (cadr points))
	  (p2 (caddr points)))
     (remhash ,2sxid *ID->2SIMPLEX*)
     (remhash 2sx *2SIMPLEX->ID*)
     (mapc (lambda (x) (gethash-update x *1SIMPLEX?* '(-1))) (pairs (2sx-points 2sx)))
     (when (= (2sx-type 2sx) 2SXSPATIAL)
       (decf (elt *SVOLS* (2sx-tmlo 2sx)))
;;       (mapc (lambda (x) (remhash x *0SIMPLEX->TIME*)) (2sx-points 2sx))
       (mapc (lambda (x) (gethash-update x *0SIMPLEX->NS2SX* '(-1))) (2sx-points 2sx))
       (mapc (lambda (x) (gethash-delete (car x) *1SIMPLEX->POINTS* (cadr x)))
	     (list (list (list p0 p1) p2)
		   (list (list p0 p2) p1)
		   (list (list p1 p2) p0))))))

(defun remove-2simplices (2sxids)
  (dolist (2sxid 2sxids)
    (remove-2simplex 2sxid)))

;; Get the points in the three spatial 2-simplex that neighbor the spatial 2-simplex
;; with the given points
;; We read-in (p0 p1 p2) and return ((foo p1 p2) (foo p0 p2) (foo p0 p1))

(defun 2sx-neighbors-points (points)
  (mapcar (lambda (x) (append (filter (lambda (y) (/= y (car x)))
				      (gethash (cadr x) *1SIMPLEX->POINTS*))
			      (cadr x)))
	  (triples-decompose points)))

(defun line-2sx-points (points)
  (mapcar (lambda (x) (cons x points)) (gethash points *1SIMPLEX->POINTS*)))

;; Introspection
(defun show-2simplex->id-store ()
  (maphash #'(lambda (2sx 2sxid) (format t "~A [~A]~%" 2sx 2sxid)) *2SIMPLEX->ID*))

(defun show-id->2simplex-store ()
  (maphash #'(lambda (2sxid 2sx) (format t "[~A] ~A~%" 2sxid 2sx)) *ID->2SIMPLEX*))

;;
;; Code for 3 simplices
;;

;; 3simplex
;; (type tmlo tmhi (p0 p1 p2 p3) (n0 n1 n2 n3) (t0 t1 t2 t3))
;; type = 1,2,3 tm[lo|hi] - [lo|hi] spatial time slice
;; pj = index of one of the points in the simplex
;;      order is not specified except that low points come before high
;; nj = id of the 3sx that borders on the face that does not have pj
;; tj = id of the 2sx that does not have pj

;; There are several versions of make-3simplex
;; Versions 3, 4, and 5 are currently in use
;; The versions behave as follows:
;; internal : Just performs side-effects. Only called directly when loading
;;          : from file.
;; 1 : DEPRECATED superseded by v5
;; 2 : DEPRECATED superseded by v5
;; 3 : (type tmlo tmhitmp p0tmp p1tmp p2tmp p3tmp) does a very particular
;;   : wrapping during init
;; 4 : DEPRECATED superseded by internal
;; 5 : (simplex-data) = (type, tmlo, tmhi, (points))
;;   : returned as first part of move data

;; The v5 family returns the created id. The other families always return nil.

;; CA: the manner of passing arguments is awful

(defun make-3simplex (type tmlo tmhi p0 p1 p2 p3)
  (error "make-3simplex is deprecated. Does not work properly with HL hash tables.")
  (let ((t0type nil) (t1type nil) (t2type nil) (t3type nil) (sx3id (next-3simplex-id)))
    (ecase type
      (1 (setf t0type 2SXSPATIAL) (setf t1type 2SX12) (setf t2type 2SX12)(setf t3type 2SX12))
      (2 (setf t0type 2SX12) (setf t1type 2SX12) (setf t2type 2SX21)(setf t3type 2SX21))
      (3 (setf t0type 2SX21) (setf t1type 2SX21) (setf t2type 2SX21)(setf t3type 2SXSPATIAL)))
    (setf (gethash sx3id *ID->3SIMPLEX*)
	  (list type tmlo tmhi
		(list p0 p1 p2 p3)
		(list 0 0 0 0)
		(list (make-2simplex t0type tmlo tmhi p1 p2 p3)
		      (make-2simplex t1type tmlo tmhi p0 p2 p3)
		      (make-2simplex t2type tmlo tmhi p0 p1 p3)
		      (make-2simplex t3type tmlo tmhi p0 p1 p2))))
    sx3id))

;; same as above except the points are "packed" into a list; useful during the 
;; moves, when the points of the simplex are computed via append / unions / 
;; intersections in the form of a list
(defun make-3simplex-v2 (type tmlo tmhi pts)
  (error "make-3simplex-v2 is deprecated. Does not work properly with HL hash tables.")
  (let ((t0type nil) 
	(t1type nil) 
	(t2type nil) 
	(t3type nil) 
	(sx3id (next-3simplex-id)))
    (ecase type
      (1 (setf t0type 2SXSPATIAL t1type 2SX12 t2type 2SX12 t3type 2SX12))
      (2 (setf t0type 2SX12 t1type 2SX12 t2type 2SX21 t3type 2SX21))
      (3 (setf t0type 2SX21 t1type 2SX21 t2type 2SX21 t3type 2SXSPATIAL)))
    (setf (gethash sx3id *ID->3SIMPLEX*)
	  (list type tmlo tmhi
		(copy-list pts)
		(list 0 0 0 0)
		(list (make-2simplex t0type tmlo tmhi 
				     (nth 1 pts) (nth 2 pts) (nth 3 pts))
		      (make-2simplex t1type tmlo tmhi 
				     (nth 0 pts) (nth 2 pts) (nth 3 pts))
		      (make-2simplex t2type tmlo tmhi 
				     (nth 0 pts) (nth 1 pts) (nth 3 pts))
		      (make-2simplex t3type tmlo tmhi 
				     (nth 0 pts) (nth 1 pts) (nth 2 pts)))))
    sx3id))

;; The side effects necessary to create a three-simplex
;; Updates various hash tables.
(defun make-3simplex-internal (type tmlo tmhi p0 p1 p2 p3 n0 n1 n2 n3 sx3id)
  (let ((t0type nil) (t1type nil) (t2type nil) (t3type nil))
    (ecase type
      (1 (setf t0type 2SXSPATIAL t1type 2SX12 t2type 2SX12 t3type 2SX12))
      (2 (setf t0type 2SX12 t1type 2SX12 t2type 2SX21 t3type 2SX21))
      (3 (setf t0type 2SX21 t1type 2SX21 t2type 2SX21 t3type 2SXSPATIAL)))
    (when (= type 3SX22)               
      (gethash-update (list p0 p1) *1SIMPLEX->N3SX22UD* '(1 0))
      (gethash-update (list p2 p3) *1SIMPLEX->N3SX22UD* '(0 1))
      )
    (setf (gethash sx3id *ID->3SIMPLEX*)
	  (list type tmlo tmhi
		(list p0 p1 p2 p3)
		(list n0 n1 n2 n3)
		(list (if (= type 1) 
			  (make-2simplex t0type tmhi tmhi p1 p2 p3)
			  (make-2simplex t0type tmlo tmhi p1 p2 p3))
		      (make-2simplex t1type tmlo tmhi p0 p2 p3)
		      (make-2simplex t2type tmlo tmhi p0 p1 p3)
		      (if (= type 3)
			  (make-2simplex t3type tmlo tmlo p0 p1 p2)
			  (make-2simplex t3type tmlo tmhi p0 p1 p2)))))))


;; this version is used only during initialization. If periodic b.c. are 
;; specified, it adjusts the points on the final time slice, since the t=T 
;; slice is identified with t=0 slice.
(defun make-3simplex-v3 (type tmlo tmhitmp p0tmp p1tmp p2tmp p3tmp)
  (let ((sx3id (next-3simplex-id))
	(p0 p0tmp) 
	(p1 p1tmp) 
	(p2 p2tmp) 
	(p3 p3tmp) 
	(tmhi tmhitmp))
    (when (and (string= BCTYPE "PERIODIC") (= NUM-T tmhitmp))
      (setf tmhi 0)
      (cond ((= 1 type)
	     (decf p1 (* N0-PER-SLICE NUM-T)) 
	     (decf p2 (* N0-PER-SLICE NUM-T)) 
	     (decf p3 (* N0-PER-SLICE NUM-T)))
	    ((= 2 type)
	     (decf p2 (* N0-PER-SLICE NUM-T)) 
	     (decf p3 (* N0-PER-SLICE NUM-T)))
	    ((= 3 type)
	     (decf p3 (* N0-PER-SLICE NUM-T)))))
    (make-3simplex-internal type tmlo tmhi p0 p1 p2 p3 0 0 0 0 sx3id)))

;; this version is used for loading the simplex data from file
;; Deprecated! Use make-3simplex-internal
(defun make-3simplex-v4 (type tmlo tmhi p0 p1 p2 p3 n0 n1 n2 n3 sx3id)
  (error "make-3simplex-v4 is deprecated. Does not work properly with HL hash tables.")
  (let ((t0type nil) (t1type nil) (t2type nil) (t3type nil))
    (ecase type
      (1 (setf t0type 2SXSPATIAL t1type 2SX12 t2type 2SX12 t3type 2SX12))
      (2 (setf t0type 2SX12 t1type 2SX12 t2type 2SX21 t3type 2SX21))
      (3 (setf t0type 2SX21 t1type 2SX21 t2type 2SX21 t3type 2SXSPATIAL)))
    (mapc (lambda (x) (gethash-push x *1SIMPLEX->3SIMPLICES* sx3id)) (pairs (list p0 p1 p2 p3)))
    (setf (gethash sx3id *ID->3SIMPLEX*)
	  (list type tmlo tmhi
		(list p0 p1 p2 p3)
		(list n0 n1 n2 n3)
		(list (if (= type 1) 
			  (make-2simplex t0type tmhi tmhi p1 p2 p3)
			  (make-2simplex t0type tmlo tmhi p1 p2 p3))
		      (make-2simplex t1type tmlo tmhi p0 p2 p3)
		      (make-2simplex t2type tmlo tmhi p0 p1 p3)
		      (if (= type 3)
			  (make-2simplex t3type tmlo tmlo p0 p1 p2)
			  (make-2simplex t3type tmlo tmhi p0 p1 p2)))))))

;; a replacement for v2
(defun make-3simplex-v5 (simplex-data)
  (let* ((type (first simplex-data))
	 (tmlo (second simplex-data))
	 (tmhi (third simplex-data))
	 (pts (fourth simplex-data))
	 (p0 (nth 0 pts))
	 (p1 (nth 1 pts))
	 (p2 (nth 2 pts))
	 (p3 (nth 3 pts))
	 (sx3id (next-3simplex-id)))
    (make-3simplex-internal type tmlo tmhi p0 p1 p2 p3 0 0 0 0 sx3id)
    sx3id))

;; simplex-data-list = ((typ tmlo tmhi (p0 p1 p2 p3)) (typ tmlo tmhi (p0 p1 p2 p3))...)
;; the ids of the simplices are returned
(defun make-3simplices-in-bulk (simplex-data-list)
  (let ((3sxids nil))
    (dolist (simplex-data simplex-data-list)
      (push (make-3simplex-v5 simplex-data) 3sxids))
    3sxids))

;; Macro for accessing fields
(defmacro 3sx-type (sx) `(nth 0 ,sx))
(defmacro 3sx-tmlo (sx) `(nth 1 ,sx))
(defmacro 3sx-tmhi (sx) `(nth 2 ,sx))
(defmacro 3sx-points (sx) `(nth 3 ,sx))
(defmacro 3sx-sx3ids (sx) `(nth 4 ,sx))
(defmacro 3sx-sx2ids (sx) `(nth 5 ,sx))
(defmacro 3sx-lopts (sx) `(subseq (3sx-points ,sx) 0 (3sx-type ,sx)))
(defmacro 3sx-hipts (sx) `(subseq (3sx-points ,sx) (3sx-type ,sx)))

(defmacro get-3simplex (sxid)
  `(gethash ,sxid *ID->3SIMPLEX*))

(defmacro remove-3simplex (sx3id)
  `(let ((3sx (gethash ,sx3id *ID->3SIMPLEX*))) 
     (when (= (3sx-type 3sx) 3SX22)
       (gethash-update (3sx-lopts 3sx) *1SIMPLEX->N3SX22UD* '(-1 0))
       (gethash-update (3sx-hipts 3sx) *1SIMPLEX->N3SX22UD* '(0 -1)))
     (remhash ,sx3id *ID->3SIMPLEX*)
     (recycle-3simplex-id ,sx3id)))

(defmacro remove-3simplices (3sxids)
  `(dolist (3sxid ,3sxids)
     (remove-3simplex 3sxid)))

(defun show-id->3simplex-store ()
  (maphash #'(lambda (3sxid 3sx) (format t "[~A] ~A~%" 3sxid 3sx)) 
	   *ID->3SIMPLEX*))

;; CA this needs a more specific name
(defmacro nth-point (sx n)
  `(nth ,n (3sx-points ,sx)))

;;
;; Operations on three simplices
;;

;; Takes two simplex ids and tells them about each other if they share a 
;; single triangle.
(defun connect-3simplices (sx1id sx2id)
  (let ((sx1 nil) (sx2 nil))
    (when (and (setf sx1 (get-3simplex sx1id)) (setf sx2 (get-3simplex sx2id)))
      (let ((2sxlinkid (intersection (3sx-sx2ids sx1) (3sx-sx2ids sx2))))
	(assert (< (length 2sxlinkid) 2))
	(when (= 1 (length 2sxlinkid)) ; Can ASSERT not > 1
	  (let ((pos1 (position (first 2sxlinkid) (3sx-sx2ids sx1)))
		(pos2 (position (first 2sxlinkid) (3sx-sx2ids sx2))))
	    (setf (nth pos1 (3sx-sx3ids sx1)) sx2id 
		  (nth pos2 (3sx-sx3ids sx2)) sx1id)))))))

;; Do these simplices know that they're connected to each other?
(defmacro 3simplices-connected? (sxid1 sxid2)
  `(let ((sx1 nil) (sx2 nil))
     (and (setf sx1 (get-3simplex ,sxid1)) (setf sx2 (get-3simplex ,sxid2))
	  (find ,sxid1 (3sx-sx3ids sx2)) (find ,sxid2 (3sx-sx3ids sx1)))))

;; if the 3-simplices are connected, returns the id of the linking 2-simplex 
;; else returns 0
(defmacro link-id (sxid1 sxid2)
  `(let ((sx1 nil) (sx2 nil) (link nil))
     (if (and (setf sx1 (get-3simplex ,sxid1)) 
	      (setf sx2 (get-3simplex ,sxid2))
	      (setf link (intersection (3sx-sx2ids sx2) (3sx-sx2ids sx1))))
	 (first link)
	 0)))

(defun neighbors-of-type (sx type)
  (let ((nbors nil)
	(nsx nil)
	(nids (3sx-sx3ids sx)))
    (for (n 0 3)
      (when (and (setf nsx (get-3simplex (nth n nids))) (= type (3sx-type nsx)))
	(pushnew (nth n nids) nbors)))
    nbors))


;;
;; Operations on lists of simplices
;; CA These should never be used on long lists
;;

;; An O(n^2) operation that connects any possible pair in the list.
;; CA You're Doing It Wrong^{TM}
(defun connect-3simplices-within-list (sx1ids)
  (for (n 0 (1- (length sx1ids)))
       (for (m (1+ n) (1- (length sx1ids)))
	    (connect-3simplices (nth n sx1ids) (nth m sx1ids)))))

;; As above, but O(mn)
(defun connect-3simplices-across-lists (sx1ids sx2ids)
  (dolist (sx1id sx1ids)
    (dolist (sx2id sx2ids)
      (connect-3simplices sx1id sx2id))))

;;
;; Linear searches for simplices
;; CA: These should never be used in important pathways
;;

;; Gets a list of the ids of the three-simplices that connect these
;; time slices. Linear search of all simplices.
;; CA: You're Doing It Wrong
(defun get-simplices-in-sandwich (tlo thi)
  (let ((sxids '()))
    (maphash #'(lambda (id sx)
		 (when (and (= (3sx-tmlo sx) (mod tlo NUM-T)) 
			    (= (3sx-tmhi sx) (mod thi NUM-T)))
		   (push id sxids)))
	     *ID->3SIMPLEX*)
    sxids))

;; Like above but filters for type
;; CA stop repeating yourself
(defun get-simplices-in-sandwich-of-type (tlo thi typ)
  (let ((sxids '()))
    (maphash #'(lambda (id sx)
		 (when (and (= (3sx-tmlo sx) (mod tlo NUM-T)) 
			    (= (3sx-tmhi sx) (mod thi NUM-T))
			    (= (3sx-type sx) typ))
		   (push id sxids)))
	     *ID->3SIMPLEX*)
    sxids))

;; CA ...
(defun get-2simplices-in-sandwich-of-type (tlo thi typ)
  (let ((sxids '()))
    (maphash #'(lambda (id sx)
		 (when (and (= (2sx-tmlo sx) (mod tlo NUM-T)) 
			    (= (2sx-tmhi sx) (mod thi NUM-T))
			    (= (2sx-type sx) typ))
		   (push id sxids)))
	     *ID->2SIMPLEX*)
    sxids))

;; Returns (1ids 2ids 3ids) where 1ids is a list of (1,3) ids in the sandwich
;; etc.
;; CA ...
(defun get-simplices-in-sandwich-ordered-by-type (tlo thi)
  (let ((1ids nil) (2ids nil) (3ids nil))
    (maphash #'(lambda (id sx)
		 (when (and (= (3sx-tmlo sx) (mod tlo NUM-T)) 
			    (= (3sx-tmhi sx) (mod thi NUM-T)))
		   (ecase (3sx-type sx)
		     (1 (push id 1ids))
		     (2 (push id 2ids))
		     (3 (push id 3ids)))))
	     *ID->3SIMPLEX*)
    (values 1ids 2ids 3ids)))

(defun get-simplices-of-type (typ)
  (let ((sxids '()))
    (maphash #'(lambda (id sx)
		 (if (= (3sx-type sx) typ)
		     (push id sxids)))
	     *ID->3SIMPLEX*)
    sxids))

(defun count-simplices-of-type (typ)
  (let ((count 0))
    (maphash #'(lambda (id sx)
		 (declare (ignore id))
		 (if (= (3sx-type sx) typ)
		     (incf count)))
	     *ID->3SIMPLEX*)
    count))

(defun count-simplices-in-sandwich (tlo thi)
  (let ((count 0))
    (maphash #'(lambda (id sx)
		 (declare (ignore id))
		 (when (and (= (3sx-tmlo sx) (mod tlo NUM-T)) 
			    (= (3sx-tmhi sx) (mod thi NUM-T)))
		     (incf count)))
	     *ID->3SIMPLEX*)
    count))

(defun count-simplices-of-all-types ()
  (let ((1count 0) (2count 0) (3count 0))
    (maphash #'(lambda (id sx)
		 (declare (ignore id))
		 (ecase (3sx-type sx)
		   (1 (incf 1count))
		   (2 (incf 2count))
		   (3 (incf 3count))))
	     *ID->3SIMPLEX*)
    (list 1count 2count 3count (+ 1count 2count 3count))))

;;
;; Sandwich level bulk operations
;; CA these should never be used in important pathways
;;

;; in a given sandwich
;; (1,3) can be connected to a (2,2) and cannot be connected to a (3,1)
;; a (2,2) can be connected to a (1,3) and a (3,1)
;; a (3,2) can be connected to a (2,2) and cannot be connected to a (1,3)
(defun connect-simplices-in-sandwich (tlo thi)
  (multiple-value-bind (1ids 2ids 3ids) 
      (get-simplices-in-sandwich-ordered-by-type tlo thi)
    (connect-3simplices-across-lists 1ids 2ids)
    (connect-3simplices-across-lists 2ids 3ids)
    (connect-3simplices-within-list 1ids)
    (connect-3simplices-within-list 2ids)
    (connect-3simplices-within-list 3ids)))

(defun connect-simplices-in-adjacent-sandwiches (tl tm th)
  (connect-3simplices-across-lists 
   (get-simplices-in-sandwich-of-type tl tm 1)
   (get-simplices-in-sandwich-of-type tm th 3)))

;;
;; Correctness checks
;; These are not currently in use
;; CA: we should make a formal set of conditionally compiled checks
;;

(defun check-13-and-31 (tlo thi)
  (let ((13ids (get-simplices-in-sandwich-of-type tlo thi 1))
	(31ids (get-simplices-in-sandwich-of-type tlo thi 3))
	(problem-ids '()))
    (dolist (s13 13ids)
      (dolist (d13 13ids)
	(when (and (/= s13 d13) (set-equal? 
				 (subseq (3sx-points (get-3simplex s13)) 1)
				 (subseq (3sx-points (get-3simplex d13)) 1)))
	  (push (list s13 d13) problem-ids))))
    (dolist (s31 31ids)
      (dolist (d31 31ids)
	(when (and (/= s31 d31) (set-equal? 
				 (subseq (3sx-points (get-3simplex s31)) 0 3)
				 (subseq (3sx-points (get-3simplex d31)) 0 3)))
	  (push (list s31 d31) problem-ids))))
    problem-ids))

(defun check-all-slices-for-problem-simplices ()
  (for (ts 0 (- NUM-T 1))
       (format t "slice ~A has ~A problem simplices~%" ts 
	       (check-13-and-31 ts (1+ ts)))
       (finish-output)))

(defun check-all-slices-for-simplices-with-missing-neighbors ()
  (let ((problem-ids '()))
    (maphash #'(lambda (id sx)
		 (for (n 0 3)
		   (if (= 0 (nth n (3sx-sx3ids sx)))
		       (push id problem-ids))))
	     *ID->3SIMPLEX*)
    problem-ids))

;;
;; For writing to and reading from files
;;

(defun save-spacetime-to-file (outfile)
  "Writes the current spacetime to a file."
  (format outfile "~A ~A ~A ~A ~A ~A ~A ~A ~A ~A ~A ~A ~A ~A ~A ~A ~A ~A ~A~%" 
	  BCTYPE STOPOLOGY NUM-T N-INIT *LAST-USED-POINT* *LAST-USED-3SXID* 
	  N0 N1-SL N1-TL N2-SL N2-TL N3-TL-31 N3-TL-22 *eps* *k0* *k3* *alpha* *lambda* *mu*)
  (maphash #'(lambda (k v)
	       (format outfile "~A ~A ~A ~A ~A ~A ~A ~A ~A ~A ~A ~A~%" 
		       (3sx-type v) (3sx-tmlo v) (3sx-tmhi v)
		       (nth-point v 0) (nth-point v 1) 
		       (nth-point v 2) (nth-point v 3)
		       (nth 0 (3sx-sx3ids v)) (nth 1 (3sx-sx3ids v)) 
		       (nth 2 (3sx-sx3ids v)) (nth 3 (3sx-sx3ids v)) k))
	   *ID->3SIMPLEX*))

(defun parse-parameters-line (line)
  (with-input-from-string (s line)
    (let ((data (loop
		   :for num := (read s nil nil)
		   :while num
		   :collect num)))
      (setf BCTYPE (nth 0 data) STOPOLOGY (nth 1 data) NUM-T (nth 2 data) 
	    N-INIT (nth 3 data) *LAST-USED-POINT* (nth 4 data) 
	    *LAST-USED-3SXID* (nth 5 data) N0 (nth 6 data) N1-SL (nth 7 data) 
	    N1-TL (nth 8 data) N2-SL (nth 9 data) N2-TL (nth 10 data)
	    N3-TL-31 (nth 11 data) N3-TL-22 (nth 12 data) *eps* (nth 13 data))
      (set-k0-k3-alpha-lambda-mu (nth 14 data) (nth 15 data) (nth 16 data) (nth 17 data) (nth 18 data))
      (setf *SVOLS* (make-array num-t :initial-element 0)))))

(defun make-3simplex-from-linedata (data)
  (make-3simplex-internal (nth 0 data) (nth 1 data) (nth 2 data) (nth 3 data) 
			  (nth 4 data) (nth 5 data) (nth 6 data) (nth 7 data)
			  (nth 8 data) (nth 9 data) (nth 10 data) (nth 11 data)))

(defun parse-simplex-data-line (line func)
  (with-input-from-string (s line)
    (let ((data (loop
		   :for num := (read s nil nil)
		   :while num
		   :collect num)))
      (funcall func data))))

(defun load-spacetime-from-file (infile &optional (init-func (lambda () ())) (line-func 'make-3simplex-from-linedata))
  (parse-parameters-line (read-line infile nil))
  (funcall init-func)
  (loop for line = (read-line infile nil)
     while line do (parse-simplex-data-line line line-func)))

;;
;; Functions for spatial 2-simplices
;; Spatial 2-simplex is a spatial triangle; this information is needed for 
;; computing the spectral and hausdorff dimensions of the spatial slice, among 
;; other things.
;; CA I don't see any of this used in production
;; The only call to it is from generate-data-v3 in montecarlo.lisp
;; DO NOT USE THESE without making sure they update the right hash tables
;;

(defun make-s2simplex (triangle-time triangle-points)
  (let ((stid (next-s2simplex-id)))
    (setf (gethash stid *ID->SPATIAL-2SIMPLEX*) 
	  (list triangle-time (copy-list triangle-points) (list 0 0 0)))
    stid))

;; this version is used for loading the simplex data from file
(defun make-s2simplex-v2 (s2sxid s2sxtm p0 p1 p2 n0 n1 n2)
  (setf (gethash s2sxid *ID->SPATIAL-2SIMPLEX*)
	(list s2sxtm (list p0 p1 p2) (list n0 n1 n2))))

(defmacro get-s2simplex (stid) `(gethash ,stid *ID->SPATIAL-2SIMPLEX*))
(defmacro s2sx-time (s2sx) `(first ,s2sx))
(defmacro s2sx-points (s2sx) `(second ,s2sx))
(defmacro s2sx-sx2ids (s2sx) `(third ,s2sx))


(defun connect-s2simplices (st1id st2id)
  (let ((st1 nil) (st2 nil))
    (when (and (setf st1 (get-s2simplex st1id)) 
	       (setf st2 (get-s2simplex st2id)))
      (let* ((points1 (s2sx-points st1))
	     (points2 (s2sx-points st2))
	     (line (intersection points1 points2)))
	(when (= 2 (length line))
	  (let ((pos1 (position (first (set-difference points1 line)) points1))
		(pos2 (position (first (set-difference points2 line)) points2)))
	    (setf (nth pos1 (s2sx-sx2ids st1)) st2id 
		  (nth pos2 (s2sx-sx2ids st2)) st1id)))))))

(defun connect-s2simplices-within-list (sx1ids)
  (for (n 0 (1- (length sx1ids)))
       (for (m (1+ n) (1- (length sx1ids)))
	    (connect-s2simplices (nth n sx1ids) (nth m sx1ids)))))

(defun count-s2simplices-in-slice (ts)
  "counts the number of spaeclike 2-simplices (i.e. spatial triangles) in spatial slice ts"
  (let ((count 0))
    (maphash #'(lambda (id sx)
		 (declare (ignore id))
		 (when (= (s2sx-time sx) (mod ts NUM-T)) 
		   (incf count)))
	     *ID->SPATIAL-2SIMPLEX*)
    count))

(defun parse-s2simplex-parameters-line (line)
  "parses line for parameters but does nothing with them. The parameters are
stored in the 2-simplex data file for purposes of identification only"
  (with-input-from-string (s line)))

(defun parse-s2simplex-data-line (line)
  "parses line for spatial 2-simplex data"
  (with-input-from-string (s line)
    (let ((data (loop
		   :for num := (read s nil nil)
		   :while num
		   :collect num)))
      (make-s2simplex-v2 (nth 0 data) (nth 1 data) (nth 2 data) (nth 3 data) 
			 (nth 4 data) (nth 5 data) (nth 6 data) (nth 7 data)))))

(defun load-s2simplex-data-from-file (infile)
  "loads spatial 2-simplex data from infile"
  (parse-s2simplex-parameters-line (read-line infile nil))
  (clrhash *ID->SPATIAL-2SIMPLEX*)
  (loop for line = (read-line infile nil)
     while line do (parse-s2simplex-data-line line)))

(defun save-s2simplex-data-to-file (outfile)
  "saves the spatial 2-simplex data to outfile"
  (format outfile "~A ~A ~A ~A ~A ~A ~A ~A~%" BCTYPE STOPOLOGY NUM-T N-INIT 
	  N2-SL (hash-table-count *ID->SPATIAL-2SIMPLEX*) *k0* *k3*)
  (maphash #'(lambda (k v)
	       (let ((pts (s2sx-points v))
		     (nbors (s2sx-sx2ids v)))
		 (format outfile "~A ~A ~A ~A ~A ~A ~A ~A~%"
			 k (s2sx-time v) (nth 0 pts) (nth 1 pts) (nth 2 pts) 
			 (nth 0 nbors) (nth 1 nbors) (nth 2 nbors))))
	   *ID->SPATIAL-2SIMPLEX*))

(defun 3sx2p1->s2sx2p1 ()
  "3sx2p1->s2sx2p1 generates the spatial 2-simplex information for each spatial slice 
from the 3-simplex data for the entire spacetime."
  (clrhash *ID->SPATIAL-2SIMPLEX*)
  (setf *LAST-USED-S2SXID* 0)
  (for (ts 0 (1- NUM-T))
       (let ((31simplices (get-simplices-in-sandwich-of-type ts (1+ ts) 3)) ;; list of ids
	     (spatial-triangles '()))
	 (dolist (31simplex 31simplices)
	   (push (make-s2simplex ts (3sx-lopts (get-3simplex 31simplex))) 
		 spatial-triangles))
	 (connect-s2simplices-within-list spatial-triangles))))

(defun generate-s2sx2p1-files (3sx2p1files)
  "3sx2p1files is a list of .3sx2p1 files. For each file in this list, this
function generates a .s2sx2p1 file. The prefix for the .3sx2p1 and the .s2sx2p1
file are identical"
  (loop for line = (read-line 3sx2p1files nil)
       while line do
       (let ((outfilename (change-file-suffix line "s2sx2p1")))
	 (with-open-file (indatafile line :direction :input)
	   (load-spacetime-from-file indatafile))
	 (3sx2p1->s2sx2p1)
	 (clrhash *ID->3SIMPLEX*)
	 (sb-ext:gc :gen 1000)
	 (with-open-file (outdatafile outfilename :direction :output)
	   (save-s2simplex-data-to-file outdatafile))
	 (format t "finished 3sx2p1->s2sx2p1 for ~A at ~A~%"
		 outfilename (cdt-now-str)))))
	     ;; try-a->b methods returns the following list, IFF the move can be 
;; successfully made (new3sxdata nbors old3sxids old2sxids fvector)
;; Otherwise it returns nil

;; The parameters are as follows
;; new3sxdata : The data necessary to construct the new simplices
;;            : ((type, tmlo, tmhi, (points)), ...)
;; nbors : The 3simplex ids that neighbor the subcomplex that was modified
;; old3sxids : The 3simplices in the modified subcomplex
;; old2sxids : The triangles internal to the modified subcomplex
;; fvector : the delta vector for our modification
;; point id : for the point, if any, created or deleted

(defmacro movedata-new3sxdata (movedata)
  `(nth 0 ,movedata))

(defmacro movedata-nbors (movedata)
  `(nth 1 ,movedata))

(defmacro movedata-old3sxids (movedata)
  `(nth 2 ,movedata))

(defmacro movedata-old2sxids (movedata)
  `(nth 3 ,movedata))

(defmacro movedata-fvector (movedata)
  `(nth 4 ,movedata))

(defmacro movedata-pointid (movedata)
  `(nth 5 ,movedata))

;; 2011-04-12 
;; sxdata should be changed to 
;; new3sxids 3sxnbors old3sxids old2sxids fvect news2sxids s2sxnbors olds2sxids
;; the last 3 elements will be null for 2->3 and 3->2 moves because new 
;; spatial triangles come into existence and old ones are deleted for the
;; 2->6, 6->2 and the 4->4 moves

(defun 2plus1move (mtype sxdata)
  ;; *Went here before
;;  (when (= mtype 26MTYPE)
;;      (setf (gethash (movedata-pointid sxdata) *0SIMPLEX->NS2SX*) 
;;	    (list (2sx-tmlo (gethash (car (movedata-old2sxids sxdata)) *ID->2SIMPLEX*)) '(0))))
;;  (when (= mtype 62MTYPE)
;;      (remhash (movedata-pointid sxdata) *0SIMPLEX->NS2SX*))
  (let ((new3sxids (make-3simplices-in-bulk (first sxdata))))
    (connect-3simplices-within-list new3sxids)
    (connect-3simplices-across-lists new3sxids (second sxdata))
    ;; * Next two lines moved from previous position
    (remove-3simplices (third sxdata))
    (remove-2simplices (fourth sxdata))
    (update-f-vector (fifth sxdata))))

;;
;; The code that suppports 2->6 moves
;;

;; Looks at the chosen simplex and identifies all the useful subcomplexes
;; In our case, there's only one relevant subcomplex that you can possibly
;; be part of

(defun 2->6-subcomplex (sxid)
  (let ((subcmplx nil)
	(sx nil))
    (when (setf sx (get-3simplex sxid))
      (cond ((= 1 (3sx-type sx))
	     (when (/= 0 (nth 0 (3sx-sx3ids sx)))
	       (push (list sxid (nth 0 (3sx-sx3ids sx))) subcmplx)))
	    ((= 3 (3sx-type sx))
	     (when (/= 0 (nth 3 (3sx-sx3ids sx)))
	       (push (list (nth 3 (3sx-sx3ids sx)) sxid) subcmplx)))))
    subcmplx))

;; We'll try the operation on each of the supplied subcomplexes, hoping 
;; that it works on one. As soon as it works, we'll break and return
;; the results.
			  
(defun try-2->6 (sxid)
  (dolist (curr (2->6-subcomplex sxid))
    (let* ((sx13 (get-3simplex (first curr)))
	   (sx31 (get-3simplex (second curr)))
	   (old-internal-triangle (nth 0 (3sx-sx2ids sx13)))
	   (nbors (set-difference (union (3sx-sx3ids sx13) (3sx-sx3ids sx31)) 
				  curr))
	   (lopt (nth-point sx13 0))
	   (hipt (nth-point sx31 3))
	   (bigtr (2sx-points (get-2simplex old-internal-triangle)))
	   (13tmlo (3sx-tmlo sx13)) (13tmhi (3sx-tmhi sx13))
	   (31tmlo (3sx-tmlo sx31)) (31tmhi (3sx-tmhi sx31))
	   (newpt (next-pt))
	   (newsxdata 
	    (list (list 1 13tmlo 13tmhi 
			(list lopt newpt (first bigtr) (second bigtr)))
		  (list 1 13tmlo 13tmhi 
			(list lopt newpt (second bigtr) (third bigtr)))
		  (list 1 13tmlo 13tmhi 
			(list lopt newpt (third bigtr) (first bigtr)))
		  (list 3 31tmlo 31tmhi 
			(list (first bigtr) (second bigtr) newpt hipt))
		  (list 3 31tmlo 31tmhi 
			(list (second bigtr) (third bigtr) newpt hipt))
		  (list 3 31tmlo 31tmhi 
			(list (third bigtr) (first bigtr) newpt hipt)))))
      (return-from try-2->6 (list newsxdata nbors curr 
				  (list old-internal-triangle) DF26 newpt)))))

;;
;; Code that supports a 6->2 move
;;

(defun 6->2-subcomplex (sxid)
  "returns a list of the form ((13id1 13id2 13id3 31id3 31id2 31id1)...)"
  (let ((sx nil)
	(subcmplx nil))
    (when (setf sx (get-3simplex sxid))
      (cond ((= 3SX13 (3sx-type sx))
	     (when (= (elt *SVOLS* (3sx-tmhi sx)) 4)
	       (return-from 6->2-subcomplex))
	     (let ((31id (nth 0 (3sx-sx3ids sx)))
		   (31sx nil))
	       (when (setf 31sx (get-3simplex 31id))
		 (let ((13nbors (neighbors-of-type sx 1)))
		   (unless (< (length 13nbors) 2)
		     (do-tuples/c (currid nextid) 13nbors
		       (let ((curr nil) (next nil))
			 (when (and (setf curr (get-3simplex currid)) (setf next (get-3simplex nextid))
				    (3simplices-connected? currid nextid)
				    (3simplices-connected? (nth 0 (3sx-sx3ids curr))
							   (nth 0 (3sx-sx3ids next)))
				    (3simplices-connected? (nth 0 (3sx-sx3ids curr)) 31id)
				    (3simplices-connected? (nth 0 (3sx-sx3ids next)) 31id))
			   (pushnew (list sxid currid nextid 
					  (nth 0 (3sx-sx3ids next)) 
					  (nth 0 (3sx-sx3ids curr))
					  31id)
				    subcmplx :test #'set-equal?)))))))))
	    ((= 3SX31 (3sx-type sx))
	     (when (= (elt *SVOLS* (3sx-tmlo sx)) 4)
	       (return-from 6->2-subcomplex))
	     (let ((13id (nth 3 (3sx-sx3ids sx)))
		   (13sx nil))
	       (when (setf 13sx (get-3simplex 13id))
		 (let ((31nbors (neighbors-of-type sx 3)))
		   (unless (< (length 31nbors) 2)
		     (do-tuples/c (currid nextid) 31nbors
		       (let ((curr nil) (next nil))
			 (when (and (setf curr (get-3simplex currid)) (setf next (get-3simplex nextid))
				    (3simplices-connected? currid nextid)
				    (3simplices-connected? (nth 3 (3sx-sx3ids curr))
							   (nth 3 (3sx-sx3ids next)))
				    (3simplices-connected? (nth 3 (3sx-sx3ids curr)) 13id)
				    (3simplices-connected? (nth 3 (3sx-sx3ids next)) 13id))
			   (pushnew (list 13id (nth 3 (3sx-sx3ids next)) (nth 3 (3sx-sx3ids curr))
					  currid nextid sxid) 
				    subcmplx :test #'set-equal?)))))))))))
    subcmplx))
	    
(defun try-6->2 (sxid)
  (let ((subcmplx (6->2-subcomplex sxid))
	(old-internal-triangles nil)
	(new-internal-triangle nil)
	(old-pt nil)
	(nbors nil)
	(newsxdata nil)
	(lopt nil)
	(hipt nil)
	(13tmlo nil) (13tmhi nil) (31tmlo nil) (31tmhi nil))
    (unless (null subcmplx)
      (dolist (curr subcmplx)
	(setf old-internal-triangles (list (link-id (first curr) (sixth curr))
					   (link-id (first curr) (second curr))
					   (link-id (first curr) (third curr))
					   (link-id (second curr) (third curr))
					   (link-id (second curr) (fifth curr))
					   (link-id (third curr) (fourth curr))
					   (link-id (fourth curr) (fifth curr))
					   (link-id (fourth curr) (sixth curr))
					   (link-id (fifth curr) (sixth curr))))
	(setf nbors (set-difference (unions (3sx-sx3ids (get-3simplex (first curr)))
					    (3sx-sx3ids (get-3simplex (second curr)))
					    (3sx-sx3ids (get-3simplex (third curr)))
					    (3sx-sx3ids (get-3simplex (fourth curr)))
					    (3sx-sx3ids (get-3simplex (fifth curr)))
					    (3sx-sx3ids (get-3simplex (sixth curr))))
				    (list 0 (first curr) (second curr) (third curr) 
					  (fourth curr) (fifth curr) (sixth curr))))
	(setf 13tmlo (3sx-tmlo (get-3simplex (first curr))))
	(setf 13tmhi (3sx-tmhi (get-3simplex (first curr))))
	(setf 31tmlo (3sx-tmlo (get-3simplex (sixth curr))))
	(setf 31tmhi (3sx-tmhi (get-3simplex (sixth curr))))
	(setf lopt (nth-point (get-3simplex (first curr)) 0))
	(setf hipt (nth-point (get-3simplex (sixth curr)) 3))
	(setf old-pt (car (intersections 
			   (2sx-points (get-2simplex (first old-internal-triangles)))
			   (2sx-points (get-2simplex (fifth old-internal-triangles)))
			   (2sx-points (get-2simplex (sixth old-internal-triangles))))))
	(setf new-internal-triangle 
	      (list 0 13tmlo 13tmhi
		    (set-difference 
		     (unions (2sx-points (get-2simplex (first old-internal-triangles)))
			     (2sx-points (get-2simplex (fifth old-internal-triangles)))
			     (2sx-points (get-2simplex (sixth old-internal-triangles))))
		     (list old-pt))))
	(setf newsxdata (list (list 1 13tmlo 13tmhi (cons lopt (2sx-points new-internal-triangle))) 
			      (list 3 31tmlo 31tmhi (append (2sx-points new-internal-triangle) (list hipt)))))
	(return-from try-6->2 (list newsxdata nbors curr old-internal-triangles DF62 old-pt))))))
  
;;
;; Code that supports a 4->4 move
;;

(defun 4->4-subcomplex (sxid)
  "returns a list of the form ((13id1 13id2 31id2 31id1)...)"
  (let ((sx nil)
	(subcmplx nil))
    (when (setf sx (get-3simplex sxid))
      (cond ((= 1 (3sx-type sx))
	     (let ((31id (nth 0 (3sx-sx3ids sx))) (31sx nil))
	       (when (setf 31sx (get-3simplex 31id))
		 (let ((13nbors (neighbors-of-type sx 1)))
		   (dolist (13nbor 13nbors)
		     (when (3simplices-connected? (nth 0 (3sx-sx3ids (get-3simplex 13nbor))) 31id)
		       (pushnew (list sxid 13nbor (nth 0 (3sx-sx3ids (get-3simplex 13nbor))) 31id)
				subcmplx :test #'set-equal?)))))))
	    ((= 3 (3sx-type sx))
	     (let ((13id (nth 3 (3sx-sx3ids sx))) (13sx nil))
	       (when (setf 13sx (get-3simplex 13id))
		 (let ((31nbors (neighbors-of-type sx 3)))
		   (dolist (31nbor 31nbors)
		     (when (3simplices-connected? (nth 3 (3sx-sx3ids (get-3simplex 31nbor))) 13id)
		       (pushnew (list 13id (nth 3 (3sx-sx3ids (get-3simplex 31nbor))) 31nbor sxid)
				subcmplx :test #'set-equal?)))))))))
    subcmplx))

(defun try-4->4 (sxid)
  (let ((subcmplx (4->4-subcomplex sxid))
	(old-internal-triangles nil)
	(new-internal-sl-triangle-1 nil) (new-internal-sl-triangle-2 nil)
	(new-internal-tl-triangle-1 nil) (new-internal-tl-triangle-2 nil)
	(shared nil) (unshared nil) (nbors nil) (newsxdata nil) (lopt nil) (hipt nil)
	(13tmlo nil) (13tmhi nil) (31tmlo nil) (31tmhi nil))
    (unless (null subcmplx)
      (dolist (curr subcmplx)
	(setf old-internal-triangles (list (link-id (first curr) (fourth curr)) ;; spacelike
					   (link-id (second curr) (third curr)) ;; spacelike
					   (link-id (first curr) (second curr)) ;; timelike
					   (link-id (third curr) (fourth curr)))) ;; timelike
	(setf 13tmlo (3sx-tmlo (get-3simplex (first curr))))
	(setf 13tmhi (3sx-tmhi (get-3simplex (first curr))))
	(setf 31tmlo (3sx-tmlo (get-3simplex (fourth curr))))
	(setf 31tmhi (3sx-tmhi (get-3simplex (fourth curr))))
	(setf lopt (nth-point (get-3simplex (first curr)) 0))
	(setf hipt (nth-point (get-3simplex (fourth curr)) 3))
	(setf nbors (set-difference (unions (3sx-sx3ids (get-3simplex (first curr)))
					    (3sx-sx3ids (get-3simplex (second curr)))
					    (3sx-sx3ids (get-3simplex (third curr)))
					    (3sx-sx3ids (get-3simplex (fourth curr))))
				    (list 0 (first curr) (second curr) (third curr) (fourth curr)))) 
	(setf shared (intersection (2sx-points (get-2simplex (first old-internal-triangles)))
				   (2sx-points (get-2simplex (second old-internal-triangles)))))
	(setf unshared (set-exclusive-or (2sx-points (get-2simplex (first old-internal-triangles)))
					 (2sx-points (get-2simplex (second old-internal-triangles)))))
	(setf new-internal-sl-triangle-1 (list 0 13tmlo 13tmhi (append unshared (butlast shared))))
	(setf new-internal-sl-triangle-2 (list 0 13tmlo 13tmhi (append unshared (last shared))))
	(setf new-internal-tl-triangle-1 (list 1 13tmlo 13tmhi (cons lopt unshared)))
	(setf new-internal-tl-triangle-2 (list 2 31tmlo 31tmhi (append unshared (cons hipt nil))))
	
	(unless (gethash unshared *1SIMPLEX?*)
	  (setf newsxdata (list (list 1 13tmlo 13tmhi (cons lopt (2sx-points new-internal-sl-triangle-1)))
				(list 1 13tmlo 13tmhi (cons lopt (2sx-points new-internal-sl-triangle-2)))
				(list 3 31tmlo 31tmhi (append (2sx-points new-internal-sl-triangle-1) 
							      (list hipt)))
				(list 3 31tmlo 31tmhi (append (2sx-points new-internal-sl-triangle-2) 
							      (list hipt)))))
	  (return-from try-4->4 (list newsxdata nbors curr old-internal-triangles DF44 nil)))))))

;;
;; Code that supports a 2->3 move
;;

(defun 2->3-subcomplex (sxid)
  "returns a list of the form ((1or3 13or31id 22id)...) where the first number 1 or 3 tells us about the
type of the simplex participating in the move"
  (let ((sx nil)
	(subcmplx nil))
    (when (setf sx (get-3simplex sxid))
      (cond ((or (= 1 (3sx-type sx)) (= 3 (3sx-type sx)))
	     (let ((22nbors (neighbors-of-type sx 2)))
	       (dolist (22nbor 22nbors)
		 (pushnew (list (3sx-type sx) sxid 22nbor) subcmplx :test #'set-equal?))))
	    ((= 2 (3sx-type sx))
	     (let ((13nbors (append (neighbors-of-type sx 1))))
	       (dolist (13nbor 13nbors)
		 (pushnew (list 1 13nbor sxid) subcmplx :test #'set-equal?)))
	     (let ((31nbors (append (neighbors-of-type sx 3))))
	       (dolist (31nbor 31nbors)
		 (pushnew (list 3 31nbor sxid) subcmplx :test #'set-equal?))))))
    subcmplx))

;; (1 | 2 3 4) (+) (1 5 | 3 4) --> (5 | 2 3 4) (+) (1 5 | 2 3) (+) (1 5 | 2 4)
(defun 2->3-move-internal-12 (13id 22id)
  "the 2,3 move performed on a (1,3) simplex attached to a (2,2) simplex"
  (let ((13sx nil) (22sx nil))
    (when (and (setf 13sx (get-3simplex 13id)) (setf 22sx (get-3simplex 22id)))
      (let* ((old2 (link-id 13id 22id))
	     (pts234 (3sx-hipts 13sx))    ;; hi points of the 1,3 simplex
	     (pts34 (3sx-hipts 22sx))  ;; hi points of the 2,2 simplex
	     (pts15 (3sx-lopts 22sx));; lo points of the 2,2 simplex
	     (pt5 (set-difference pts15 (3sx-lopts 13sx)))
	     (pt2 (set-difference pts234 pts34))
	     (new-internal-tlt-1 (list 1 (3sx-tmlo 13sx) (3sx-tmhi 13sx) 
				       (append pt5 pt2 (butlast pts34))))
	     (new-internal-tlt-2 (list 1 (3sx-tmlo 13sx) (3sx-tmhi 13sx)
				       (append pt5 pt2 (last pts34))))
	     (new-internal-tlt-3 (list 2 (3sx-tmlo 13sx) (3sx-tmhi 13sx) (append pts15 pt2)))
	     (new-13-pts (append pt5 pts234))
	     (new-22-pts-1 (append pts15 pt2 (butlast pts34)))
	     (new-22-pts-2 (append pts15 pt2 (last pts34)))
	     (nbors (set-difference (union (3sx-sx3ids 13sx) (3sx-sx3ids 22sx)) (list 0 13id 22id)))
	     (newsxdata nil))
	(unless (gethash (list (car pt2) (car pt5)) *1SIMPLEX?*)
	  (setf newsxdata (list (list 1 (3sx-tmlo 13sx) (3sx-tmhi 13sx) new-13-pts) 
				(list 2 (3sx-tmlo 22sx) (3sx-tmhi 22sx) new-22-pts-1) 
				(list 2 (3sx-tmlo 22sx) (3sx-tmhi 22sx) new-22-pts-2)))
	  (return-from 2->3-move-internal-12 (list newsxdata nbors (list 13id 22id) (list old2) DF23 nil)))))))

;; (2 3 4 | 1) (+) (3 4 | 1 5) --> (2 3 4 | 5) (+) (2 3 | 1 5) (+) (2 4 | 1 5)
(defun 2->3-move-internal-32 (31id 22id)
  "the 2,3 move performed on a (3,1) simplex attached to a (2,2) simplex"
    (let ((31sx nil) (22sx nil))
      (when (and (setf 31sx (get-3simplex 31id)) (setf 22sx (get-3simplex 22id)))
	(let* ((old2 (link-id 31id 22id))
	       (pts234 (3sx-lopts 31sx))    ;; lo points of the 3,1 simplex
	       (pts34 (3sx-lopts 22sx))  ;; lo points of the 2,2 simplex
	       (pts15 (3sx-hipts 22sx));; hi points of the 2,2 simplex
	       (pt5 (set-difference pts15 (3sx-hipts 31sx)))
	       (pt2 (set-difference pts234 pts34))
	       (new-internal-tlt-1 (list 2 (3sx-tmlo 31sx) (3sx-tmhi 31sx) 
					 (append pt2 (butlast pts34) pt5)))
	       (new-internal-tlt-2 (list 2 (3sx-tmlo 31sx) (3sx-tmhi 31sx)
					 (append pt2 (last pts34) pt5)))
	       (new-internal-tlt-3 (list 1 (3sx-tmlo 31sx) (3sx-tmhi 31sx) (append pt2 pts15)))
	       (new-31-pts (append pts234 pt5))
	       (new-22-pts-1 (append pt2 (butlast pts34) pts15))
	       (new-22-pts-2 (append pt2 (last pts34) pts15))
	       (nbors (set-difference (union (3sx-sx3ids 31sx) (3sx-sx3ids 22sx)) (list 0 31id 22id)))
	       (newsxdata nil))
	  (unless (gethash (list (car pt2) (car pt5)) *1SIMPLEX?*)
	    (setf newsxdata (list (list 3 (3sx-tmlo 31sx) (3sx-tmhi 31sx) new-31-pts) 
				  (list 2 (3sx-tmlo 22sx) (3sx-tmhi 22sx) new-22-pts-1) 
				  (list 2 (3sx-tmlo 22sx) (3sx-tmhi 22sx) new-22-pts-2)))
	    (return-from 2->3-move-internal-32 (list newsxdata nbors (list 31id 22id) (list old2) DF23 nil)))))))

(defun try-2->3 (sxid)
  (let ((subcmplx (2->3-subcomplex sxid))
	(movedata nil))
    (unless (null subcmplx)
      (dolist (curr subcmplx)
	(cond ((= 1 (first curr))
	       (setf movedata (2->3-move-internal-12 (second curr) (third curr)))
	       (when movedata
		 (return-from try-2->3 movedata)))
	      ((= 3 (first curr))
	       (setf movedata (2->3-move-internal-32 (second curr) (third curr)))
	       (when movedata
		 (return-from try-2->3 movedata))))))))

;;
;; Code that supports a 3->2 move
;;

(defun 3->2-subcomplex (sxid)
  "returns a list of the form ((1or3 13or31id 22id1 22id2)...) where the first number 1 or 3 tells us 
about the type of the simplex participating in the move"
  (let ((sx nil)
	(subcmplx nil))
    (when (setf sx (get-3simplex sxid))
      (cond ((or (= 1 (3sx-type sx)) (= 3 (3sx-type sx)))
	     (let ((22nbors (neighbors-of-type sx 2)))
	       (dolist (22nbor 22nbors)
		 (let ((22sx (get-3simplex 22nbor)))
		   (when 22sx
		     (let ((22nborsof22nbor (neighbors-of-type 22sx 2)))
		       (dolist (22nborof22nbor 22nborsof22nbor)
			 (when (3simplices-connected? 22nborof22nbor sxid)
			   (pushnew (list (3sx-type sx) sxid 22nbor 22nborof22nbor) subcmplx 
				    :test #'set-equal?)))))))))
	    ((= 2 (3sx-type sx))
	     (let ((22nbors (neighbors-of-type sx 2)))
	       (dolist (22nbor 22nbors)
		 (let ((22sx (get-3simplex 22nbor)))
		   (when 22sx
		     (let ((13nborsof22nbor (neighbors-of-type 22sx 1)))
		       (dolist (13nborof22nbor 13nborsof22nbor)
			 (when (3simplices-connected? 13nborof22nbor sxid)
			   (pushnew (list 1 13nborof22nbor sxid 22nbor) subcmplx :test #'set-equal?))))
		     (let ((31nborsof22nbor (neighbors-of-type 22sx 3)))
		       (dolist (31nborof22nbor 31nborsof22nbor)
			 (when (3simplices-connected? 31nborof22nbor sxid)
			   (pushnew (list 3 31nborof22nbor sxid 22nbor) subcmplx :test #'set-equal?)))))))))))
    subcmplx))
;; (5 | 2 3 4) (+) (1 5 | 2 3) (+) (1 5 | 2 4) --> (1 | 2 3 4) (+) (1 5 | 3 4)
(defun 3->2-move-internal-122 (13id 22id1 22id2)
  "the (3,2) move performed on a (1,3) simplex attached to two (2,2) simplices"
  (let ((13sx nil) (22sx1 nil) (22sx2 nil))
    (when (and (setf 13sx (get-3simplex 13id)) 
	       (setf 22sx1 (get-3simplex 22id1))
	       (setf 22sx2 (get-3simplex 22id2)))
      (let* ((old2s (list (link-id 13id 22id1) (link-id 13id 22id2) (link-id 22id1 22id2)))
	     (pts234 (3sx-hipts 13sx))
	     (pts15 (3sx-lopts 22sx1))
	     (pt1 (set-difference pts15 (3sx-lopts 13sx)))
	     (pts34 (set-exclusive-or (3sx-hipts 22sx1) (3sx-hipts 22sx2)))
	     (new-internal-triangle (list 1 (3sx-tmlo 13sx) (3sx-tmhi 13sx) (append pt1 pts34)))
	     (new-13-pts (append pt1 pts234))
	     (new-22-pts (append pts15 pts34))
	     (nbors (set-difference (unions (3sx-sx3ids 13sx) (3sx-sx3ids 22sx1) (3sx-sx3ids 22sx2)) 
				    (list 0 13id 22id1 22id2)))
	     (newsxdata nil))
	(unless (gethash new-internal-triangle *2SIMPLEX->ID*)
	  (setf newsxdata (list (list 1 (3sx-tmlo 13sx) (3sx-tmhi 13sx) new-13-pts) 
				(list 2 (3sx-tmlo 22sx1) (3sx-tmhi 22sx2) new-22-pts)))
	  (return-from 3->2-move-internal-122 (list newsxdata nbors (list 13id 22id1 22id2) old2s DF32 nil)))))))
	
(defun 3->2-move-internal-322 (31id 22id1 22id2)
  "the (3,2) move performed on a (3,1) simplex attached to two (2,2) simplices"
  (let ((31sx nil) (22sx1 nil) (22sx2 nil))
    (when (and (setf 31sx (get-3simplex 31id)) 
	       (setf 22sx1 (get-3simplex 22id1))
	       (setf 22sx2 (get-3simplex 22id2)))
      (let* ((old2s (list (link-id 31id 22id1) (link-id 31id 22id2) (link-id 22id1 22id2)))
	     (pts234 (3sx-lopts 31sx)) ;; lo points of 3,1
	     (pts15 (3sx-hipts 22sx1))   ;; hi points of 2,2
	     (pt1 (set-difference pts15 (3sx-hipts 31sx))) ;; 2,2 hi - 3,1 hi
	     (pts34 (set-exclusive-or (3sx-lopts 22sx1) (3sx-lopts 22sx2)))
	     (new-internal-triangle (list 2 (3sx-tmlo 31sx) (3sx-tmhi 31sx) (append pts34 pt1)))
	     (new-31-pts (append pts234 pt1))
	     (new-22-pts (append pts34 pts15))
	     (nbors (set-difference (unions (3sx-sx3ids 31sx) (3sx-sx3ids 22sx1) (3sx-sx3ids 22sx2)) 
				    (list 0 31id 22id1 22id2)))
	     (newsxdata nil))
	(unless (gethash new-internal-triangle *2SIMPLEX->ID*)
	  (setf newsxdata (list (list 3 (3sx-tmlo 31sx) (3sx-tmhi 31sx) new-31-pts) 
				(list 2 (3sx-tmlo 22sx1) (3sx-tmhi 22sx2) new-22-pts))) 
	  (return-from 3->2-move-internal-322 (list newsxdata nbors (list 31id 22id1 22id2) old2s DF32 nil)))))))

(defun try-3->2 (sxid)
  (let ((subcmplx (3->2-subcomplex sxid))
	(movedata nil))
    (unless (null subcmplx)
      (dolist (curr subcmplx)
	(cond ((= 1 (first curr))
	       (setf movedata (3->2-move-internal-122 (second curr) (third curr) (fourth curr)))
	       (when movedata
		 (return-from try-3->2 movedata)))
	      ((= 3 (first curr))
	       (setf movedata (3->2-move-internal-322 (second curr) (third curr) (fourth curr)))
	       (when movedata 
		 (return-from try-3->2 movedata))))))))
;...........................................................................................................
; cdt-2plus1-initialization.lisp
;...........................................................................................................

(defun initialize-S2-triangulation ()
  (for (n 0 (1- (/ NUM-T 2)))
       (dolist (fourpts S2-1/2-31)
	 (make-3simplex-v3 3 (* 2 n) (1+ (* 2 n))                       ;-----o------------ t = 1
			   (+ (* 2 n N0-PER-SLICE) (first fourpts))     ;    / \
			   (+ (* 2 n N0-PER-SLICE) (second fourpts))    ;   /   \
			   (+ (* 2 n N0-PER-SLICE) (third fourpts))     ;  /     \
			   (+ (* 2 n N0-PER-SLICE) (fourth fourpts))))  ;-o-------o-------- t = 0
       (dolist (fourpts S2-1/2-22)
	 (make-3simplex-v3 2 (* 2 n) (1+ (* 2 n))                       ;-----o-----o------ t = 1
			   (+ (* 2 n N0-PER-SLICE) (first fourpts))     ;    /     / 
			   (+ (* 2 n N0-PER-SLICE) (second fourpts))    ;   /     /     
			   (+ (* 2 n N0-PER-SLICE) (third fourpts))     ;  /     /  
			   (+ (* 2 n N0-PER-SLICE) (fourth fourpts))))  ;-o-----o---------- t = 0

       (dolist (fourpts S2-1/2-13)
	 (make-3simplex-v3 1 (* 2 n) (1+ (* 2 n))                       ;-o-------o-------- t = 1
			   (+ (* 2 n N0-PER-SLICE) (first fourpts))     ;  \     /
			   (+ (* 2 n N0-PER-SLICE) (second fourpts))    ;   \   /
			   (+ (* 2 n N0-PER-SLICE) (third fourpts))     ;    \ /
			   (+ (* 2 n N0-PER-SLICE) (fourth fourpts))))  ;-----o------------ t = 0
       (dolist (fourpts S2-1/2-31)
	 (make-3simplex-v3 1 (1+ (* 2 n)) (2+ (* 2 n))                      ;-o-------o-------- t = 2
			   (+ (* 2 n N0-PER-SLICE) (fourth fourpts))        ;  \     /
			   (+ (* 2 (+ n 1) N0-PER-SLICE) (first fourpts))   ;   \   /
			   (+ (* 2 (+ n 1) N0-PER-SLICE) (second fourpts))  ;    \ /
			   (+ (* 2 (+ n 1) N0-PER-SLICE) (third fourpts)))) ;-----o------------ t = 1
       (dolist (fourpts S2-1/2-22)
	 (make-3simplex-v3 2 (1+ (* 2 n)) (2+ (* 2 n))                      ;-----o-----o------ t = 2
			   (+ (* 2 n N0-PER-SLICE) (third fourpts))         ;      \     \  
			   (+ (* 2 n N0-PER-SLICE) (fourth fourpts))        ;       \     \
			   (+ (* 2 (+ n 1) N0-PER-SLICE) (first fourpts))   ;        \     \
			   (+ (* 2 (+ n 1) N0-PER-SLICE) (second fourpts))));---------o-----o-- t = 1
       (dolist (fourpts S2-1/2-13)
	 (make-3simplex-v3 3 (1+ (* 2 n)) (2+ (* 2 n))                      ;-----o------------ t = 2
			   (+ (* 2 n N0-PER-SLICE) (second fourpts))        ;    / \
			   (+ (* 2 n N0-PER-SLICE) (third fourpts))         ;   /   \
			   (+ (* 2 n N0-PER-SLICE) (fourth fourpts))        ;  /     \
			   (+ (* 2 (+ n 1) N0-PER-SLICE) (first fourpts)))));-o-------o-------- t = 1
  
  (when (string= BCTYPE "PERIODIC")

    (for (ts 0 (- NUM-T 1))
	 (connect-simplices-in-sandwich ts (1+ ts) )
	 (connect-simplices-in-adjacent-sandwiches ts (+ ts 1) (+ ts 2)))
  
    (set-last-used-pt (* NUM-T N0-PER-SLICE))
    
    (set-f-vector (* NUM-T N0-PER-SLICE)                               ; N0
		  (* NUM-T N1-SL-PER-SLICE)                            ; N1-SL
		  (* NUM-T N1-TL-PER-SLICE)                            ; N1-TL
		  (* NUM-T N2-SL-PER-SLICE)                            ; N2-SL
		  (* NUM-T N2-TL-PER-SLICE)                            ; N2-TL
		  (* NUM-T (+ N3-TL-13-PER-SLICE N3-TL-31-PER-SLICE))  ; N3-TL-13 + N3-TL-31
		  (* NUM-T N3-TL-22-PER-SLICE)))                       ; N3-TL-22

  (when (string= BCTYPE "OPEN")
    (error "open boundary conditions not yet implemented")))

(defun initialize-T2-triangulation (num-time-intervals boundary-conditions)
  (setf NUM-T num-time-intervals)
  (setf BCTYPE (string-upcase boundary-conditions))
  (error "IxT2 and S1xT2 spacetimes not yet implemented"))

(defun set-T-slices-with-V-volume (&key
				   num-time-slices
				   target-volume
				   spatial-topology
				   boundary-conditions)
  (setf STOPOLOGY (string-upcase spatial-topology))
  (setf TARG-VOL-IDEAL target-volume)
  (setf NUM-T num-time-slices)
  (setf BCTYPE (string-upcase boundary-conditions)))

(defun initialize (&optional (version nil))
  (when version
    (error "Initialize does not support multiple versions."))

  (reset-global-state)

  (when (string= STOPOLOGY "S2")
    (initialize-S2-triangulation))

  (when (string= STOPOLOGY "T2")
    (initialize-T2-triangulation))

;;  (maphash #'(lambda (id sx)
;;	       (declare (ignore sx))
;;	       (update-subcomplexes-for-3simplex id))
;;	   *ID->3SIMPLEX*)

  (format t "initial count = ~A ~%" (count-simplices-of-all-types))

  (loop named tv
     do
       (dolist (id23 (get-simplices-of-type 2))
	 (let ((movedata nil))
	   (when (setf movedata (try-2->3 id23))
	     (2plus1move 23MTYPE movedata)))
	 (if (> (N3) TARG-VOL-IDEAL)
	     (return-from tv)))
       (dolist (id26 (get-simplices-of-type 3))
	 (let ((movedata nil))
	   (when (setf movedata (try-2->6 id26))
	     (2plus1move 26MTYPE movedata)))
	 (if (> (N3) TARG-VOL-IDEAL)
	     (return-from tv)))

       (dolist (id23 (get-simplices-of-type 2))
	 (let ((movedata nil))
	   (when (setf movedata (try-2->3 id23))
	     (2plus1move 23MTYPE movedata)))
	 (if (> (N3) TARG-VOL-IDEAL)
	     (return-from tv))))


  (format t "final count = ~A~%" (count-simplices-of-all-types))

  (setf N-INIT (N3)))
; cdt-2plus1-montecarlo.lisp

(defun try-move (sxid mtype)
  "Try the move of the given numerical type on the given simplex."
  (ecase mtype
    (0 (try-2->6 sxid))
    (1 (try-2->3 sxid))
    (2 (try-4->4 sxid))
    (3 (try-3->2 sxid))
    (4 (try-6->2 sxid))))

;; Deprecated
(defun random-move (nsweeps)
  (loop :for sweepnum :from 1 :to nsweeps
     do
     (let* ((id (random *LAST-USED-3SXID*))
	    (mtype (select-move))
	    (sx (get-3simplex id))
	    (movedata nil))

       (incf CURRENT-MOVE-NUMBER)
       (when (and sx (setf movedata (try-move id mtype)))
	 (2plus1move movedata))
       (when (= 0 (mod sweepnum 1000))
	 (format t "finished ~A of ~A sweeps with count ~A~%" sweepnum nsweeps 
		 (count-simplices-of-all-types))
	 (finish-output)))))
;; Calculate the action of a given spacetime
(let* ((2alpha+1 (+ (* 2 *alpha*) 1))
	 (4alpha+1 (+ (* 4 *alpha*) 1))
	 (4alpha+2 (+ (* 4 *alpha*) 2))
	 (3alpha+1 (+ (* 3 *alpha*) 1))
	 (arcsin-1 (asin (/ (* *-i* (wrsqrt (* 8 2alpha+1))) 4alpha+1)))
	 (arccos-1 (acos (/ *-i* (wrsqrt (* 3 4alpha+1)))))
	 (arccos-2 (acos (/ -1 4alpha+1)))
	 (arccos-3 (acos (/ 2alpha+1 4alpha+1))))
  (defun action-depr (num1-sl num1-tl num3-31 num3-22)
    (let ((out (- (* *k* (+ (- (* *2pi/i* num1-sl) 
			       (* *2/i* num3-22 arcsin-1) 
			       (* *3/i* num3-31 arccos-1))
			    (* (wrsqrt *alpha*) (- (* 2 pi num1-tl) 
						   (* 4 num3-22 arccos-2) 
						   (* 3 num3-31 arccos-3)))))
		  (* (/ *litL* 12) (+ (* num3-22 (wrsqrt 4alpha+2)) 
				      (* num3-31 (wrsqrt 3alpha+1)))))))
      out)))


;; Calculate the action of a given spacetime
(let* ((2alpha+1 (+ (* 2 *alpha*) 1))
       (4alpha+1 (+ (* 4 *alpha*) 1))
       (4alpha+2 (+ (* 4 *alpha*) 2))
       (3alpha+1 (+ (* 3 *alpha*) 1))
       (arcsin-1 (asin (/ (* *-i* (wrsqrt (* 8 2alpha+1))) 4alpha+1)))
       (arccos-1 (acos (/ *-i* (wrsqrt (* 3 4alpha+1)))))
       (arccos-2 (acos (/ -1 4alpha+1)))
       (arccos-3 (acos (/ 2alpha+1 4alpha+1)))
       (A (* *2pi/i* *k*))
       (B (* (wrsqrt *alpha*) 2 pi *k*))
       (C (- (+ (* *3/i* arccos-1 *k*) (* (wrsqrt *alpha*) 3 arccos-3 *k*) (* (/ *litL* 12) (wrsqrt 3alpha+1)))))
       (D (- (+ (* *2/i* arcsin-1 *k*) (* (wrsqrt *alpha*) 4 arccos-2 *k*) (* (/ *litL* 12) (wrsqrt 4alpha+2)))))
       (side-effect (format t "REGGE COUPLING :: ~A ~A ~A ~A ~%" A B C D)))
  (defun action (num1-sl num1-tl num3-31 num3-22)
    (let ((out (+ (* A num1-sl) (* B num1-tl) (* C num3-31) (* D num3-22))))
;;      (assert (= out (action-depr num1-sl num1-tl num3-31 num3-22)))
;;      (format t "~A ~A ~%" (+ (* *2/i* arcsin-1) (* (wrsqrt *alpha*) 4 arccos-2))
;;	      (/ (wrsqrt 4alpha+2) 12))
;;      (format t "REGGE COUPLING :: ~A ~A ~A ~A ~%" A B C D)
      (format t "~A vs. ~A ~%" out (funcall *regge-action* num1-sl num1-tl num3-31 num3-22))
      out)))

;; Given a spatial edge and the index of the time-slice that contains it
;; Return the counts of (up down) 2,2 simplices that contain that edge
(defun N22UD (edge)
  ;; Note: This is dangerous. If the edge isn't listed in the hash table
  ;; we assume that it is simply surrounded by 13 simplices and is not
  ;; part of any 22s. Hopefully true, but it could supress errors.
  (or (gethash edge *1SIMPLEX->N3SX22UD*) '(0 0)))

(defun N22UD-triag (triag)
  (reduce (lambda (x y) (mapcar #'+ x y)) (mapcar (lambda (x) (or (gethash x *1SIMPLEX->N3SX22UD*) '(0 0))) (pairs triag))))

;; Given a vertex. Return the count of spatial 2simplices that contain it
(defun NS2SX (vertex)
  (car (gethash vertex *0SIMPLEX->NS2SX*)))

(let* ((aa (* -2 (/ *k* 4) (- 1 *lambda*) (sqrt *eta*) *a*))
       (bb (* 2 (/ *k* 4) *mu* (sqrt *eta*) (/ 1 *a*)))
       (A (+ (* bb (/ -8 3) (expt pi 2)) (* aa 2 (expt (- (* 3 pi) (* 6 *theta31*)) 2))))
       (B (* bb (/ 2 9) (expt pi 2)))
       (C (* aa -2 (expt *theta22* 2))))
  (defun horava-action-move-26 (movedata)
    (let* ((s2sx (gethash (car (movedata-old2sxids movedata)) *ID->2SIMPLEX*))
;;	   (err (assert (= (2sx-type s2sx) 2SXSPATIAL)))
	   (points (2sx-points s2sx))
;;	   (err (assert (= time (2sx-tmhi s2sx))))
	   (n22 (mapcar (lambda (x) (N22UD x)) (pairs points)))
	   (n22u (mapcar 'car n22))
	   (n22d (mapcar 'cadr n22)))
;;	   (err (assert (> (sum n22u) 2)))
;;	   (err (assert (> (sum n22d) 2))))
;;      (format t "~A ~%" (mapcar 'NS2SX points))
      (+ A
	 (* B
	    (sum (mapcar 'NS2SX points)))
	 (* C
	    (+ (sum (mapcar (lambda (y) (apply #'* y)) (pairs n22u)))
	       (sum (mapcar (lambda (z) (apply #'* z)) (pairs n22d)))))))))

(let* ((aa (* 2 (/ *k* 4) (- 1 *lambda*) (sqrt *eta*) *a*))
       (bb (* 2 (/ *k* 4) *mu* (sqrt *eta*) (/ 1 *a*)))
       (A (+ (* bb (/ 10 3) (expt pi 2)) (* aa 2 (expt (- (* 3 pi) (* 6 *theta31*)) 2))))
       (B (* bb (/ -2 9) (expt pi 2)))
       (C (* aa -2 (expt *theta22* 2))))
  (defun horava-action-move-62 (movedata)
    (let* ((new3sx (car (movedata-new3sxdata movedata)))
;;	   (err (assert (/= (3sx-type new3sx) 3SX22)))
	   (points (if (= (3sx-type new3sx) 3SX13) (3sx-hipts new3sx) (3sx-lopts new3sx)))
	   (n22 (mapcar (lambda (x) (N22UD x)) (pairs points)))
	   (n22u (mapcar 'car n22))
	   (n22d (mapcar 'cadr n22)))
;;	   (err (assert (> (sum n22u) 2)))
;;	   (err (assert (> (sum n22d) 2))))
;;      (format t "~A ~%" (mapcar 'NS2SX points))
      (+ A
	 (* B
	    (sum (mapcar 'NS2SX points)))
	 (* C
	    (+ (sum (mapcar (lambda (x) (apply #'* x)) (pairs n22u)))
	       (sum (mapcar (lambda (x) (apply #'* x)) (pairs n22d)))))))))

(let* ((aa (* -2 (/ *k* 4) (- 1 *lambda*)))
       (A (* aa (+ (expt *theta22* 2) (* 2 *theta22* (- (* 3 pi) (* 6 *theta31*))))))
       (B (* aa -2 (expt *theta22* 2))))
  (defun horava-action-move-32 (movedata)
    (let* ((new3sxdata (movedata-new3sxdata movedata))
	   (new3sx (if (= (3sx-type (car new3sxdata)) 3SX22) (cadr new3sxdata) (car new3sxdata)))
;;	   (err (assert (/= (3sx-type new3sx) 3SX22)))
	   (points (if (= (3sx-type new3sx) 3SX13) (3sx-hipts new3sx) (3sx-lopts new3sx)))
	   (n22 (mapcar (lambda (x) (N22UD x)) (pairs points)))
	   (n22u (mapcar 'car n22)))
      (+ A
	 (* B (sum n22u))))))

(let* ((aa (* -2 (/ *k* 4) (- 1 *lambda*)))
       (A (* aa (+ (expt *theta22* 2) (* -2 *theta22* (- (* 3 pi) (* 6 *theta31*))))))
       (B (* aa 2 (expt *theta22* 2))))
  (defun horava-action-move-23 (movedata)
    (let* ((new3sx (car 
		    (filter (lambda (x) (/= (3sx-type x) 3SX22)) (movedata-new3sxdata movedata))))
	   (points (if (= (3sx-type new3sx) 3SX13) (3sx-hipts new3sx) (3sx-lopts new3sx)))
	   (n22 (mapcar (lambda (x) (N22UD x)) (pairs points)))
	   (n22u (mapcar 'car n22)))
      (+ A
	 (* B (sum n22u))))))

(let* ((A (* 2 (/ *k* 4) *mu* (sqrt *eta*) (/ 1 *a*) (/ 2 9) (expt pi 2)))
       (B (* -2 (/ *k* 4) (- 1 *lambda*) (sqrt *eta*) *a* 2 (expt *theta22* 2))))
  (defun horava-action-move-44 (movedata)
    (let* ((olds2sxs (filter (lambda (x) (= (2sx-type x) 2SXSPATIAL)) (mapcar (lambda (x) (gethash x *ID->2SIMPLEX*)) (movedata-old2sxids movedata))))
	   (oldpoints (mapcar (lambda (x) (2sx-points x)) olds2sxs))
	   (pts35 (apply 'intersections oldpoints))
	   (pts24 (differences (apply 'unions oldpoints) pts35))
	   (pt3 (car pts35))
	   (pt5 (cadr pts35))
	   (pt2 (car pts24))
	   (pt4 (cadr pts24))
;;	   (time (2sx-tmlo (car olds2sxs))) ;; Triangle is spatial so tmlo = tmhi
	   (nud25 (N22UD (list pt2 pt5)))
	   (nud45 (N22UD (list pt4 pt5)))
	   (nud34 (N22UD (list pt3 pt4)))
	   (nud23 (N22UD (list pt2 pt3))))
;;      (format t "3: ~A 5: ~A 2: ~A 4: ~A ~%" (NS2SX pt3) (NS2SX pt5) (NS2SX pt2) (NS2SX pt4))
      (+ (* A (+ 2 (-(NS2SX pt3)) (-(NS2SX pt5)) (NS2SX pt2) (NS2SX pt4)))
	 (* B (+ (* (car nud25) (car nud45)) (* (car nud34) (car nud23)) 
		 (- (* (car nud25) (car nud23))) (- (* (car nud45) (car nud34)))
		 (* (cadr nud25) (cadr nud45)) (* (cadr nud34) (cadr nud23))
		 (- (* (cadr nud25) (cadr nud23))) (- (* (cadr nud45) (cadr nud34)))))))))

(let* ((aa (* -2 (/ *k* 4) (- 1 *lambda*) (sqrt *eta*) *a*))
       (bb (* 2 (/ *k* 4) *mu* (sqrt *eta*) (/ 1 *a*)))
       (A (+ (* bb (/ 2 3) (expt pi 2)) (* aa 2 (expt (- (* 3 pi) (* 6 *theta31*)) 2))))
       (B (* bb 4 (expt pi 2)))
       (C (* aa -2 (expt *theta22* 2))))
  (defun horava-action-move-26-new (movedata)
    (let* ((s2sx (gethash (car (movedata-old2sxids movedata)) *ID->2SIMPLEX*))
;;	   (err (assert (= (2sx-type s2sx) 2SXSPATIAL)))
	   (points (2sx-points s2sx))
;;	   (err (assert (= time (2sx-tmhi s2sx))))
	   (n22 (mapcar (lambda (x) (N22UD x)) (pairs points)))
	   (n22u (mapcar 'car n22))
	   (n22d (mapcar 'cadr n22)))
;;	   (err (assert (> (sum n22u) 2)))
;;	   (err (assert (> (sum n22d) 2))))
;;      (format t "~A ~A ~A ~A ~A ~%" aa bb A B C)
;;      (format t "~A ~%" (mapcar 'NS2SX points))
      (+ A
	 (* B
	    (sum (mapcar (lambda (x) (- (/ 1 (+ x 1)) (/ 1 x))) (mapcar 'NS2SX points))))
	 (* C
	    (+ (sum (mapcar (lambda (y) (apply #'* y)) (pairs n22u)))
	       (sum (mapcar (lambda (z) (apply #'* z)) (pairs n22d)))))))))

(let* ((aa (* 2 (/ *k* 4) (- 1 *lambda*) (sqrt *eta*) *a*))
       (bb (* 2 (/ *k* 4) *mu* (sqrt *eta*) (/ 1 *a*)))
       (A (+ (* bb (/ -2 3) (expt pi 2)) (* aa 2 (expt (- (* 3 pi) (* 6 *theta31*)) 2))))
       (B (* bb 4 (expt pi 2)))
       (C (* aa -2 (expt *theta22* 2))))
  (defun horava-action-move-62-new (movedata)
    (let* ((new3sx (car (movedata-new3sxdata movedata)))
;;	   (err (assert (/= (3sx-type new3sx) 3SX22)))
	   (points (if (= (3sx-type new3sx) 3SX13) (3sx-hipts new3sx) (3sx-lopts new3sx)))
	   (n22 (mapcar (lambda (x) (N22UD x)) (pairs points)))
	   (n22u (mapcar 'car n22))
	   (n22d (mapcar 'cadr n22)))
;;	   (err (assert (> (sum n22u) 2)))
;;	   (err (assert (> (sum n22d) 2))))
;;      (format t "~A ~%" (mapcar 'NS2SX points))
      (+ A
	 (* B
	    (sum (mapcar (lambda (x) (- (/ 1 (- x 1)) (/ 1 x))) (mapcar 'NS2SX points))))
	 (* C
	    (+ (sum (mapcar (lambda (x) (apply #'* x)) (pairs n22u)))
	       (sum (mapcar (lambda (x) (apply #'* x)) (pairs n22d)))))))))

(let* ((A (* 2 (/ *k* 4) *mu* (sqrt *eta*) (/ 1 *a*) 4 (expt pi 2)))
       (B (* -2 (/ *k* 4) (- 1 *lambda*) (sqrt *eta*) *a* 2 (expt *theta22* 2))))
  (flet ((pos (x) (- (/ 1 (+ x 1)) (/ 1 x)))
	 (neg (x) (- (/ 1 (- x 1)) (/ 1 x))))
    (defun horava-action-move-44-new (movedata)
      (let* ((olds2sxs (filter (lambda (x) (= (2sx-type x) 2SXSPATIAL)) (mapcar (lambda (x) (gethash x *ID->2SIMPLEX*)) (movedata-old2sxids movedata))))
	     (oldpoints (mapcar (lambda (x) (2sx-points x)) olds2sxs))
	     (pts35 (apply 'intersections oldpoints))
	     (pts24 (differences (apply 'unions oldpoints) pts35))
	     (pt3 (car pts35))
	     (pt5 (cadr pts35))
	     (pt2 (car pts24))
	     (pt4 (cadr pts24))
;;	     (time (2sx-tmlo (car olds2sxs))) ;; Triangle is spatial so tmlo = tmhi
	     (nud25 (N22UD (list pt2 pt5)))
	     (nud45 (N22UD (list pt4 pt5)))
	     (nud34 (N22UD (list pt3 pt4)))
	     (nud23 (N22UD (list pt2 pt3))))
;;	(format t "3: ~A 5: ~A 2: ~A 4: ~A ~%" (NS2SX pt3) (NS2SX pt5) (NS2SX pt2) (NS2SX pt4))
	(+ (* A (+ (neg (NS2SX pt3)) (neg (NS2SX pt5)) (pos (NS2SX pt2)) (pos (NS2SX pt4))))
	   (* B (+ (* (car nud25) (car nud45)) (* (car nud34) (car nud23)) 
		   (- (* (car nud25) (car nud23))) (- (* (car nud45) (car nud34)))
		   (* (cadr nud25) (cadr nud45)) (* (cadr nud34) (cadr nud23))
		   (- (* (cadr nud25) (cadr nud23))) (- (* (cadr nud45) (cadr nud34))))))))))

;; Calculates the R^2 contribution from a set of vertices with the given
;; sizes of entourage

;; Version 1 uses the volume share average prescription
(defun build-r2-pluggable ()
  (let ((A (* 2 (/ *k* 4) *mu* (/ (sqrt *eta*) *a*)))
	(B (* 4 (expt pi 2)))
	(C (* (/ -4 3) (expt pi 2)))
	(D (* (/ 1 9) (expt pi 2))))
    (defun r2-term-v1 (ns2sx-list)
      (* A (sum (mapcar (lambda (x)
			  (+ (/ B x) C (* D x)))
			ns2sx-list))))))

;; Calculates the K^2 contribution from a set of triangles with the given
;; up and down counts
;; This is the complete contribution, including all leading coefficients

;; Version 1 uses the prescription from Patrick's note of July 28, 2011
(defun build-k2-pluggable ()
  (let ((A (* -2 (/ *k* 4) (- 1 *lambda*) *a*))
	(B (- (* 3 pi) (* 6 *THETA31*)))
	(C (- *THETA22*))
	(D (* 2 (sqrt 2) (sqrt (- (* 3 *eta*) 1))))
	(E (sqrt (- (* 2 *eta*) 1))))
    (defun k2-term-v1 (n22ud-list)
      ;;    (format t "~A ~%" n22ud-list)
      (* A (sum (mapcar (lambda (x) 
			  (sum (mapcar (lambda (y)
					 (/ (expt (+ B (* C y)) 2) (+ D (* E y))))
				       x)))
			n22ud-list))))))

(defun horava-action-move-26-pluggable (movedata)
  (let* ((s2sx (gethash (car (movedata-old2sxids movedata)) *ID->2SIMPLEX*))
;;	   (err (assert (= (2sx-type s2sx) 2SXSPATIAL)))
	 (points (2sx-points s2sx))
	 (out 0))
;;	   (err (assert (= time (2sx-tmhi s2sx))))
    (setf out (+ (zero-or *mu* 
			  (let ((ns (mapcar 'NS2SX points)))
			    (- (funcall *r2-pluggable* (cons 3 (mapcar (lambda (x y) (+ x y)) ns '(1 1 1))))
			       (funcall *r2-pluggable* ns))))
		 (zero-or (- 1 *lambda*)
			  (let ((n22 (mapcar (lambda (x) (N22UD x)) (pairs points))))
			    (- (funcall *k2-pluggable* n22)
			       (funcall *k2-pluggable* (list (reduce (lambda (x y) (mapcar #'+ x y)) n22))))))))
;;    (format t "Check: ~A vs. ~A ~%" out (horava-action-move-26-new movedata))
    out))

(defun horava-action-move-62-pluggable (movedata)
  (let* ((new3sx (car (movedata-new3sxdata movedata)))
;;	   (err (assert (/= (3sx-type new3sx) 3SX22)))
	 (points (if (= (3sx-type new3sx) 3SX13) (3sx-hipts new3sx) (3sx-lopts new3sx)))
	 (out 0))
;;      (format t "~A ~%" (mapcar 'NS2SX points))
    (setf out (+ (zero-or *mu* 
			  (let ((ns (mapcar 'NS2SX points)))
			    (- (funcall *r2-pluggable* (mapcar (lambda (x y) (+ x y)) ns '(-1 -1 -1)))
			       (funcall *r2-pluggable* (cons 3 ns)))))
		 (zero-or (- 1 *lambda*) 
			  (let ((n22 (mapcar (lambda (x) (N22UD x)) (pairs points))))
			    (- (funcall *k2-pluggable* (list (reduce (lambda (x y) (mapcar #'+ x y)) n22)))
			       (funcall *k2-pluggable* n22))))))
;;    (format t "Check: ~A vs. ~A ~%" out (horava-action-move-62-new movedata))
    out))

;; Before a 44 move, the two spatial triangles share edge 35
(defun horava-action-move-44-pluggable (movedata)
  (let* ((3sxup (filter (lambda (x) (= (3sx-type x) 3SX31)) (movedata-new3sxdata movedata)))
	 (newpoints (mapcar (lambda (x) (3sx-lopts x)) 3sxup))
	 (pts24 (apply 'intersections newpoints))
	 (pts35 (differences (apply 'unions newpoints) pts24))
	 (pt3 (car pts35))
	 (pt5 (cadr pts35))
	 (pt2 (car pts24))
	 (pt4 (cadr pts24))
	 (out 0))
;;	     (time (2sx-tmlo (car olds2sxs))) ;; Triangle is spatial so tmlo = tmhi
;;	(format t "3: ~A 5: ~A 2: ~A 4: ~A ~%" (NS2SX pt3) (NS2SX pt5) (NS2SX pt2) (NS2SX pt4))
    (setf out (+ (zero-or *mu*
			  (let ((ns (mapcar 'NS2SX (list pt2 pt3 pt4 pt5))))
			    (- (funcall *r2-pluggable* (mapcar (lambda (x y) (+ x y)) ns '(1 -1 1 -1)))
			       (funcall *r2-pluggable* ns))))
		 ;;TODO FIX
		 (zero-or (- 1 *lambda*) 
			  (let ((nud25 (N22UD (list pt2 pt5)))
				(nud45 (N22UD (list pt4 pt5)))
				(nud34 (N22UD (list pt3 pt4)))
				(nud23 (N22UD (list pt2 pt3))))
			    (- (funcall *k2-pluggable* (list (mapcar (lambda (x y) (+ x y)) nud23 nud34)
						 (mapcar (lambda (x y) (+ x y)) nud25 nud45)))
			       (funcall *k2-pluggable* (list (mapcar (lambda (x y) (+ x y)) nud23 nud25)
						 (mapcar (lambda (x y) (+ x y)) nud34 nud45))))))))
;;    (format t "Check: ~A vs. ~A ~%" out (horava-action-move-44-new movedata))
    out))

;; For variable names, we assume that the vertices in the two simplices 
;; that you perform a 2->3 move on are:
;; 3 simplex type 2: (2 3 4 5)
;; 3 simplex type 3: (1 2 3 4)
;; where 1, 2, 3 are on the same time-slice
(let ((UP 0) (DOWN 1))
  (defun horava-action-move-23-pluggable (movedata)
    (zero-or 
     (- 1 *lambda*)
     (let* ((old3sx (mapcar (lambda (x) (gethash x *ID->3SIMPLEX*)) (movedata-old3sxids movedata)))
	    (old3sx2 (if (= (3sx-type (car old3sx)) 3SX22) (car old3sx) (cadr old3sx)))
	    (old3sx3 (if (/= (3sx-type (car old3sx)) 3SX22) (car old3sx) (cadr old3sx)))
	    (direction (if (= (3sx-type old3sx3) 3SX31) UP DOWN))
	    (pt1 (car (differences (3sx-points old3sx3) (3sx-points old3sx2))))
	    (pts23 (intersections (if (= direction UP) (3sx-lopts old3sx3) (3sx-hipts old3sx3)) (3sx-points old3sx2)))
	    (pts45 (if (= direction UP) (3sx-hipts old3sx2) (3sx-lopts old3sx2)))
	    (pts123 (cons pt1 pts23))
	    (pts123-n22ud (N22UD-triag pts123))
	    (neighbors123 (2sx-neighbors-points pts123))
	    (neighbors123-n22ud (mapcar 'N22UD-triag neighbors123))
	    (neighbors45 (line-2sx-points pts45))
;;	    (side (format t "23 MOVE : ~A ~%" direction))
	    (neighbors45-n22ud (mapcar 'N22UD-triag neighbors45))
	    (incr (if (= direction UP) '(1 0) '(0 1)))
	    (decr (if (= direction UP) '(-1 0) '(0 -1))))
       (- (funcall *k2-pluggable* (append (list (mapcar (lambda (x y) (+ x y)) pts123-n22ud incr))
			      (mapcar (lambda (p q) (mapcar (lambda (x y) (+ x y)) p q)) neighbors123-n22ud (list decr incr incr))
			      (mapcar (lambda (p q) (mapcar (lambda (x y) (+ x y)) p q)) neighbors45-n22ud (list incr incr))))
	  (funcall *k2-pluggable* (append (list pts123-n22ud) neighbors123-n22ud neighbors45-n22ud)))))))

(let ((UP 0) (DOWN 1))
  (defun horava-action-move-32-pluggable (movedata)
    (zero-or 
     (- 1 *lambda*)
     (let* ((new3sxdata (movedata-new3sxdata movedata))
	    (new3sx2 (if (= (3sx-type (car new3sxdata)) 3SX22) (car new3sxdata) (cadr new3sxdata)))
	    (new3sx3 (car 
		     (filter (lambda (x) (/= (3sx-type x) 3SX22)) (movedata-new3sxdata movedata))))
	    (direction (if (= (3sx-type new3sx3) 3SX31) UP DOWN))
	    (pt1 (car (differences (3sx-points new3sx3) (3sx-points new3sx2))))
	    (pts23 (intersection (if (= direction UP) (3sx-lopts new3sx3) (3sx-hipts new3sx3)) (3sx-points new3sx2)))
	    (pts45 (if (= direction UP) (3sx-hipts new3sx2) (3sx-lopts new3sx2)))
	    (pts123 (cons pt1 pts23))
	    (pts123-n22ud (N22UD-triag pts123))
	    (neighbors123 (2sx-neighbors-points pts123))
	    (neighbors123-n22ud (mapcar 'N22UD-triag neighbors123))
	    (neighbors45 (line-2sx-points pts45))
;;	    (side (format t "32 MOVE : ~A ~%" direction))
	    (neighbors45-n22ud (mapcar 'N22UD-triag neighbors45))
	    (incr (if (= direction UP) '(1 0) '(0 1)))
	    (decr (if (= direction UP) '(-1 0) '(0 -1))))
       (- (funcall *k2-pluggable* (append (list (mapcar (lambda (x y) (+ x y)) pts123-n22ud decr))
			      (mapcar (lambda (p q) (mapcar (lambda (x y) (+ x y)) p q)) neighbors123-n22ud (list incr decr decr))
			      (mapcar (lambda (p q) (mapcar (lambda (x y) (+ x y)) p q)) neighbors45-n22ud (list decr decr))))
	  (funcall *k2-pluggable* (append (list pts123-n22ud) neighbors123-n22ud neighbors45-n22ud)))))))

(defparameter horava-action-oldr2 (list 'horava-action-move-26
					'horava-action-move-23
					'horava-action-move-44
					'horava-action-move-32
					'horava-action-move-62))

(defparameter horava-action-newr2 (list 'horava-action-move-26-new
					'horava-action-move-23
					'horava-action-move-44-new
					'horava-action-move-32
					'horava-action-move-62-new))

(defparameter horava-action-pluggable (list 'horava-action-move-26-pluggable
					'horava-action-move-23-pluggable
					'horava-action-move-44-pluggable
					'horava-action-move-32-pluggable
					'horava-action-move-62-pluggable))

(defparameter horava-action-move horava-action-pluggable)

(defun accept-move? (mtype movedata)
  "Decide whether to accept a legal move based on its MC weight."
  (let ((cur-vol (N3))
	(final-vol 0)
	(delta-action 0.0)
	(delta-damping 0.0)
	(HL-delta-action 0.0)
	(tot-neg-delta-action 0.0))
;;    (format t "MOVE :: ~A ~%" mtype)
    (cond ((= 26MTYPE mtype) ;; 2->6 move
	   (setf final-vol (+ cur-vol 4))
	   (setf delta-damping (- (damping (+ (N3) 4)) (damping (N3))))
	   (setf delta-action (- (funcall *regge-action* (+ N1-SL 3) (+ N1-TL 2) (+ N3-TL-31 4) (+ N3-TL-22 0)) 
				    (funcall *regge-action* N1-SL N1-TL N3-TL-31 N3-TL-22)))
	   (setf HL-delta-action (apply (nth 0 horava-action-move) (list movedata))))
	  ((= 23MTYPE mtype) ;; 2->3 move
	   (setf final-vol (+ cur-vol 1))
	   (setf delta-damping (- (damping (+ (N3) 1)) (damping (N3))))
	   (setf delta-action (- (funcall *regge-action* (+ N1-SL 0) (+ N1-TL 1) (+ N3-TL-31 0) (+ N3-TL-22 1)) 
				    (funcall *regge-action* N1-SL N1-TL N3-TL-31 N3-TL-22)))
	   (setf HL-delta-action (apply (nth 1 horava-action-move) (list movedata))))
	  ((= 44MTYPE mtype) ;; 4->4 move
	   (setf final-vol cur-vol)
	   (setf delta-damping (- (damping (+ (N3) 0)) (damping (N3))))
	   (setf delta-action (- (funcall *regge-action* (+ N1-SL 0) (+ N1-TL 0) (+ N3-TL-31 0) (+ N3-TL-22 0)) 
				    (funcall *regge-action* N1-SL N1-TL N3-TL-31 N3-TL-22)))
	   (setf HL-delta-action (apply (nth 2 horava-action-move) (list movedata))))
	  ((= 32MTYPE mtype) ;; 3->2 move
	   (setf final-vol (+ cur-vol -1))
	   (setf delta-damping (- (damping (+ (N3) -1)) (damping (N3))))
	   (setf delta-action (- (funcall *regge-action* (+ N1-SL 0) (+ N1-TL -1) (+ N3-TL-31 0) (+ N3-TL-22 -1)) 
				    (funcall *regge-action* N1-SL N1-TL N3-TL-31 N3-TL-22)))
	   (setf HL-delta-action (apply (nth 3 horava-action-move) (list movedata))))
	  ((= 62MTYPE mtype) ;; 6->2 move
	   (setf final-vol (+ cur-vol -4))
	   (setf delta-damping (- (damping (+ (N3) -4)) (damping (N3))))
	   (setf delta-action (- (funcall *regge-action* (+ N1-SL -3) (+ N1-TL -2) (+ N3-TL-31 -4) (+ N3-TL-22 0)) 
				 (funcall *regge-action* N1-SL N1-TL N3-TL-31 N3-TL-22)))
	   (setf HL-delta-action (apply (nth 4 horava-action-move) (list movedata)))))
    (setf tot-neg-delta-action (- (realpart (* *i* delta-action)) HL-delta-action delta-damping))
;;    (format t "Cur vol: ~A, Final vol: ~A, Target Vol: ~A ~%" cur-vol final-vol N-INIT)
;;    (format t "NEG-DELTA-ACTION :: ~$ ~% ~%" tot-neg-delta-action) ;;COMP
;;    (format t "Neg change in Regge action :: ~A ~%" (realpart (* *i* delta-action)))
;;    (format t "HL change :: ~A ~%" HL-delta-action)
    (or (> tot-neg-delta-action 0) (< (random 1.0) (* (exp tot-neg-delta-action))))))

(defun sweep ()
  "Attempt to make N-INIT moves."
  (let ((num-attempted 0))
    (while (< num-attempted N-INIT)
      (incf THRASHING)
      (let* ((sxid (random *LAST-USED-3SXID*))
	     (mtype (select-move))
	     (movedata (try-move sxid mtype))
	     (timeout 0))
;;	(format t "Randomly selected sxid: ~A from ~A~%" sxid *LAST-USED-3SXID*) ;;COMP
	(while (null movedata)
	  (incf timeout)
	  (incf THRASHING)
	  (setf sxid (random *LAST-USED-3SXID*))
	  (if *FIX-SELECTIONS*
	      (when (= timeout 500)
		(format t "Stuck on ~A ~%" mtype))
	      (setf mtype (select-move)))
	  (setf movedata (try-move sxid mtype)))
;;	  (format t "Randomly selected sxid: ~A from ~A~%" sxid *LAST-USED-3SXID*)) ;;COMP
	(incf num-attempted) ;; number-of-attempted-moves-counter for this sweep
	(incf (nth mtype ATTEMPTED-MOVES)) ;; number of moves of mtype that have been attempted
	(when (accept-move? mtype movedata)
;;	  (format t "Accepted move.~%") ;;COMP
	  (incf (nth mtype SUCCESSFUL-MOVES))
;;	  (format t "Move: ~A~%" mtype)
	  (2plus1move mtype movedata))))))

;; We have a family of generate-* functions.
;; These are the top-level functions that perform (sweep)
;; Each collects particular data for a particular use-case

;; data-console: Periodically prints simplex counts and accept ratios to console.
;;               Good for tuning k0 and k3 parameters.
;; data : Periodically re-writes a single file using save-spacetime-to-file and 
;;        a second file with the current progress
;; data-v2 : Like data, but writes a new file each time so no progress file needed.
;; data-v3 : Like data-v2, but at each checkpoint it also writes out a file using
;;         : save-s2simplex-data-to-file
;; movie-data : Periodically appends to a file a new list of count-simplices-in-sandwich
;;            : writes out a second file with progress data
;; movie-data-console : SUPERSEDED Like movie-data but puts everything on the console
;; movie-data-console-v2 : Like movie-data-console but counts spatial triangles per slice (FAST!)
;;                         and saves a spacetime at the end.
;; spacetime-and-movie-data : a combination of data and movie-data

;; Following is to be used when tuning the k0 and k3 parameters. Not during data runs.
(defun generate-data-console (&optional (start-sweep 1))
  "Runs the program, printing the status out to the console every 10 sweeps."
  (for (ns start-sweep (+ start-sweep NUM-SWEEPS -1))
       (sweep)
       (when (= 0 (mod ns 10))
	 (format t "start = ~A end = ~A current = ~A count = ~A ~A\%~%"
		 start-sweep (+ start-sweep NUM-SWEEPS -1) ns 
		 ;(count-simplices-of-all-types) (percent-tv)));SUCCESSFUL-MOVES));(accept-ratios)))
		 (count-simplices-of-all-types) (accept-ratios)))
       (finish-output)))

;; generate-data should be called after setting the values for eps, k0, k3,
;; NUM-SWEEPS and calling one of the initialize-xx-slices.
(defun generate-data (&optional (start-sweep 1))
  (let ((datafilestr (concatenate 'string (generate-filename start-sweep) 3SXEXT))
	(progfilestr (concatenate 'string (generate-filename start-sweep) PRGEXT))
	(end-sweep (+ start-sweep NUM-SWEEPS -1)))
    (for (ns start-sweep end-sweep)
	 (sweep)
	 (when (= 0 (mod ns SAVE-EVERY-N-SWEEPS))
	   (with-open-file (datafile datafilestr 
				     :direction :output
				     :if-exists :supersede)
	     (save-spacetime-to-file datafile))
	   (with-open-file (progfile progfilestr
				     :direction :output
				     :if-exists :supersede)
	     (format progfile "~A/~A/~A ~A~%"
		     start-sweep ns end-sweep (count-simplices-of-all-types)))))))

;; generate-data-v2 is similar to generate-data except it creates a fresh data file every 
;; SAVE-EVERY-N-SWEEPS. since a fresh datafile is created, there is no need to maintain a seprate progress
;; file.
(defun generate-data-v2 (&optional (start-sweep 1))
  (setf SIM-START-TIME (cdt-now-str))
  (let ((end-sweep (+ start-sweep NUM-SWEEPS -1)))
    (when (= 1 start-sweep) ;; save the initial spacetime contents if this is a brand new run
      (with-open-file (datafile (concatenate 'string (generate-filename-v2 start-sweep 0) 3SXEXT)
				:direction :output
				:if-exists :supersede)
	(save-spacetime-to-file datafile)))
    (for (ns start-sweep end-sweep)
	 (sweep)
	 (when (= 0 (mod ns SAVE-EVERY-N-SWEEPS))
	   (with-open-file (datafile (concatenate 'string (generate-filename-v2 start-sweep ns) 3SXEXT)
				     :direction :output
				     :if-exists :supersede)
	     (save-spacetime-to-file datafile))))))

;; generate-data-v3 is similar to generate-data-v2 except it also creates an additional data file every 
;; SAVE-EVERY-N-SWEEPS that contains the spatial 2-simplex information for each spatial slice.
(defun generate-data-v3 (&optional (start-sweep 1))
  (setf SIM-START-TIME (cdt-now-str))
  (let ((end-sweep (+ start-sweep NUM-SWEEPS -1)))
    (when (= 1 start-sweep) ;; save the initial spacetime contents if this is a brand new run
      (with-open-file (datafile (concatenate 'string (generate-filename-v2 start-sweep 0) 3SXEXT)
				:direction :output
				:if-exists :supersede)
	(save-spacetime-to-file datafile)))
    (for (ns start-sweep end-sweep)
	 (sweep)
	 (when (= 0 (mod ns SAVE-EVERY-N-SWEEPS))
	   (let ((filename (generate-filename-v2 start-sweep ns)))
	     (with-open-file (datafile (concatenate 'string filename 3SXEXT)
				       :direction :output
				       :if-exists :supersede)
	       (save-spacetime-to-file datafile))
	     (3sx2p1->s2sx2p1)
	     (with-open-file (datafile (concatenate 'string filename S2SXEXT)
				       :direction :output
				       :if-exists :supersede)
	       (save-s2simplex-data-to-file datafile)))))))

;; generate-movie-data saves number of simplices every SAVE-EVERY-N-SWEEPS
(defun generate-movie-data (&optional (start-sweep 1))
  (setf SAVE-EVERY-N-SWEEPS 10)
  (let ((moviefilestr (concatenate 'string (generate-filename start-sweep) MOVEXT))
	(trackfilestr (concatenate 'string (generate-filename start-sweep) PRGEXT))
	(end-sweep (+ start-sweep NUM-SWEEPS -1)))

    ;; open and close the file for :append to work
    (with-open-file (moviefile moviefilestr 
			       :direction :output
			       :if-exists :supersede)
      ;; record the initial data only if start-sweep = 1
      (when (= start-sweep 1)

	(for (ts 0 (1- NUM-T))
	     (format moviefile "~A " (count-simplices-in-sandwich ts (1+ ts))))
	(format moviefile "~%")))
    
    (for (ns start-sweep end-sweep)
	 (sweep)
	 (when (= 0 (mod ns SAVE-EVERY-N-SWEEPS))
	   (with-open-file (moviefile moviefilestr 
				      :direction :output
				      :if-exists :append)
	     (for (ts 0 (1- NUM-T))
		  (format moviefile "~A " (count-simplices-in-sandwich ts (1+ ts))))
	     (format moviefile "~%"))
	   (with-open-file (trackfile trackfilestr
				      :direction :output
				      :if-exists :supersede)
	     (format trackfile "~A/~A/~A ~A~%"
		     start-sweep ns end-sweep (count-simplices-of-all-types)))))))

(defun generate-movie-data-console (&optional (start-sweep 1))
  (when (= 1 start-sweep)
    (reset-move-counts))
  (let ((datafilestr (concatenate 'string (generate-filename start-sweep) 3SXEXT))
	(end-sweep (+ start-sweep NUM-SWEEPS -1)))
       (for (ns start-sweep end-sweep)
	   (sweep)
;;	   (when (and (> ns 5) (< (N3) (* .8 *TARG-VOL*)))
;;	     (return-from generate-movie-data-console))
	   (when (= 0 (mod ns 10))
	     (format t "~A/~A " ns end-sweep)
	     (for (ts 0 (1- NUM-T))
		  (format t "~A " (count-simplices-in-sandwich ts (1+ ts))))
	     (format t "~% ~A ~A~%~%" (f-vector) (accept-ratios))))
       (with-open-file (datafile datafilestr 
				 :direction :output
				 :if-exists :supersede)
	 (save-spacetime-to-file datafile))))

(defun generate-movie-data-console-v2 (&optional (start-sweep 1))
  (when (= 1 start-sweep)
    (reset-move-counts))
  (let ((datafilestr (concatenate 'string (generate-filename start-sweep) 3SXEXT))
	(end-sweep (+ start-sweep NUM-SWEEPS -1)))
       (for (ns start-sweep end-sweep)
	    (sweep)
;;	    (update-pattempted)
;;	    (set-k0-k3-alpha-lambda-mu *k0* *k3* *alpha* *lambda* *mu*)
;;	    (format t "PATTEMPTED :: ~A ~%" PATTEMPTED)
;;	   (when (and (> ns 5) (< (N3) (* .8 *TARG-VOL*)))
;;	     (return-from generate-movie-data-console))
	   (when (= 0 (mod ns 10))
	     (format t "~A/~A " ns end-sweep)
	     (format t "Spatial volumes: ~A ~% " *SVOLS*)
	     (format t "~% ~A ~A |~A| ~%~%" (f-vector) (accept-ratios) THRASHING))
	   (when (= 0 (mod ns SAVE-EVERY-N-SWEEPS))
	     (with-open-file (datafile datafilestr 
				       :direction :output
				       :if-exists :supersede)
	       (save-spacetime-to-file datafile))))
       (with-open-file (datafile datafilestr 
				 :direction :output
				 :if-exists :supersede)
	 (save-spacetime-to-file datafile))))

(defun generate-spacetime-and-movie-data (&optional (start-sweep 1))
  (let* ((end-sweep (+ start-sweep NUM-SWEEPS -1))
	 (datafilestr (concatenate 'string (generate-filename start-sweep) 3SXEXT))
	 (trackfilestr (concatenate 'string (generate-filename start-sweep) PRGEXT))
	 (moviefilestr (concatenate 'string (generate-filename start-sweep) MOVEXT)))
    
    ;; open and close the file, for :append to work properly
    (with-open-file (moviefile moviefilestr 
			       :direction :output
			       :if-exists :supersede)
      ;; record the initial data only if start-sweep = 1
      (when (= 1 start-sweep)
	(for (ts 0 (1- NUM-T))
	     (format moviefile "~A " (count-simplices-in-sandwich ts (1+ ts))))
	(format moviefile "~%")))
    
    (for (ns start-sweep end-sweep)
	 (sweep)
	 (when (= 0 (mod ns SAVE-EVERY-N-SWEEPS))
	   (with-open-file (datafile datafilestr 
				     :direction :output
				     :if-exists :supersede)
	     (save-spacetime-to-file datafile))
	   
	   (with-open-file (trackfile trackfilestr
				      :direction :output
				      :if-exists :supersede)
	     (format trackfile "~A/~A/~A ~A~%" start-sweep ns end-sweep (count-simplices-of-all-types)))
	   
	   (with-open-file (moviefile moviefilestr 
				      :direction :output
				      :if-exists :append)
	     (for (ts 0 (1- NUM-T))
		  (format moviefile "~A " (count-simplices-in-sandwich ts (1+ ts))))
	     (format moviefile "~%"))))))

#|
(defun calculate-order-parameter (&optional (start-sweep 1))
  (let* ((end-sweep (+ start-sweep NUM-SWEEPS -1))
	 (order-parameter 0.0)
	 (datafilestr (format nil 
			      "~A-~A-op-T~A_V~A_eps~A_kz~A_kt~A_sweeps~Ato~A.op" 
			      *topology* *boundary-conditions*
			      NUM-T N-INIT *eps* *k0* *k3* start-sweep end-sweep))
	 (trackfilestr (format nil 
			       "~A-~A-op-T~A_V~A_eps~A_kz~A_kt~A_sweeps~Ato~A.progress" 
			       *topology* *boundary-conditions*
			       NUM-T N-INIT *eps* *k0* *k3* start-sweep end-sweep)))
    (do ((ns start-sweep (1+ ns)) 
	 (tot 0.0 (incf tot (/ N3-TL-22 (N3)))))
	((> ns end-sweep) (setf order-parameter (/ tot NUM-SWEEPS)))
      (sweep)
      (when (= 0 (mod ns SAVE-EVERY-N-SWEEPS))
	(with-open-file (trackfile trackfilestr
				   :direction :output
				   :if-exists :supersede)
	  (format trackfile "start = ~A end = ~A current = ~A~%"
		  start-sweep end-sweep ns))))
    (with-open-file (datafile datafilestr
			      :direction :output
			      :if-exists :supersede)
      (format datafile "T=~A V=~A eps=~A k0=~A k3=~A start=~A end=~A op=~A~%" 
	      NUM-T N-INIT *eps* *k0* *k3* start-sweep end-sweep order-parameter))))

(defun calculate-volume-volume-correlator (&optional (start-sweep 1))
  (let* ((end-sweep (+ start-sweep NUM-SWEEPS -1))
	 (vvparams (make-array (+ NUM-T 1) :initial-element 0.0))
	 (dfilestr (format nil 
			   "~A-~A-vv-T~A_V~A_eps~A_kz~A_kt~A_sweeps~Ato~A.vv" 
			   *topology* *boundary-conditions*
			   NUM-T N-INIT *eps* *k0* *k3* start-sweep end-sweep))
	 (tfilestr (format nil 
			   "~A-~A-vv-T~A_V~A_eps~A_kz~A_kt~A_sweeps~Ato~A.prog" 
			   *topology* *boundary-conditions*
			   NUM-T N-INIT *eps* *k0* *k3* start-sweep end-sweep)))
    (do ((ns start-sweep (1+ ns)))
	((> ns end-sweep)
	 (do ((j 0 (1+ j))) ((> j NUM-T))
	   (setf (aref vvparams j) 
		 (/ (aref vvparams j) (* NUM-T NUM-T (/ NUM-SWEEPS 100))))))
      (sweep)
      (when (= 0 (mod ns SAVE-EVERY-N-SWEEPS))
	(do ((col 0 (incf col))) ((> col NUM-T))
	  (do ((ts 1/2 (1+ ts))) ((> ts NUM-T))
	    (incf (svref vvparams col)
		  (* (count-simplices-at-time ts)
		     (count-simplices-at-time-pbc (+ ts
						     (- col (/ NUM-T 2))))))))
	(with-open-file (tfile tfilestr 
			       :direction :output 
			       :if-exists :supersede)
	  (format tfile "start = ~A end = ~A current = ~A~%"
		  start-sweep end-sweep ns))))
    (with-open-file (dfile dfilestr
			   :direction :output
			   :if-exists :supersede)
      (format dfile "T=~A V=~A eps=~A k0=~A k3=~A start=~A end=~A vvp=~A~%" 
	      NUM-T N-INIT *eps* *k0* *k3* start-sweep end-sweep vvparams))))

(defun compute-spatial-slice-hausdorff-dimension ()
  "compute the hausdorff dimension of all the spatial slices")

(defun compute-thin-sandwich-hausdorff-dimension ()
  "a thin sandwich consists of two adjacent spatial slices")

(defun compute-spacetime-hausdorff-dimension ()
  "the hausdorff dimension of the entire spacetime")

(defun compute-spatial-slice-spectral-dimension ()
  "compute the spectral dimension of all the spatial slices")

(defun compute-thin-sandwich-spectral-dimension ()
  "a thin sandwich consists of two adjacent spatial slices")

(defun compute-spacetime-spectral-dimension ()
  "the spectral dimension of the entire spacetime")

|#