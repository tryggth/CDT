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
