(load "../utilities.lisp")
(load "globals.lisp")
(load "simplex.lisp")
(load "moves.lisp")
(load "initialization.lisp")
(load "montecarlo.lisp")
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

;;------------------------------------------------------------------------------
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
;;-----------------------------------------------------------------------------
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
(defparameter SAVE-EVERY-N-SWEEPS 10 "save every 10 sweeps by default")
(defparameter NUM-T 666666 
  "number of time slices --- set to a non-zero value so (mod ts NUM-T) works")
(defparameter N-INIT 0 
  "initial volume of spacetime; we try to keep the volume close to this number")
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
(defparameter *k0* 0.0)
(defparameter *k3* 0.0)
(defparameter *eps* 0.02)
(defparameter *a* 1.0)
(defparameter *alpha* -1.0)
(defparameter *i* #C(0.0 1.0)) ;; complex number i
(defparameter *-i* #C(0.0 -1.0)) ;; complex number -i
(defparameter *2/i* (/ 2 *i*))
(defparameter *2pi/i* (* *2/i* pi))
(defparameter *3/i* (/ 3 *i*))
(defparameter ROOT2 (sqrt 2.0))
(defparameter KAPPA (/ (acos (/ 1 3)) pi))
(defparameter 6ROOT2 (* 6.0 ROOT2))
(defparameter 3KAPPAMINUS1 (- (* 3 KAPPA) 1))
;;; wrsqrt is the "wick rotated" sqrt function. Basically wrsqrt(x) = -i*sqrt(-x) when x < 0 and not 
;;; i*sqrt(-x). So wrsqrt(-1) = -i
(defmacro wrsqrt (val)
  `(if (< ,val 0)
       (* -1 ,*i* (sqrt (* -1 ,val)))
       (sqrt ,val)))

(defvar action nil)

;; STOPOLOGY-BCTYPE-NUMT-NINIT-k0-k3-eps-alpha-startsweep-endsweep-hostname-currenttime
(defun generate-filename (&optional (start-sweep 1) (end-sweep (+ start-sweep NUM-SWEEPS -1)))
  (format nil "~A-~A-T~3,'0d-V~6,'0d-~A-~A-~A-~A-~9,'0d-~9,'0d-on-~A-started~A" 
	  STOPOLOGY BCTYPE NUM-T N-INIT *k0* *k3* *eps* *alpha* start-sweep end-sweep (hostname) (cdt-now-str)))

;; STOPOLOGY-BCTYPE-NUMT-NINIT-k0-k3-eps-alpha-startsweep-currsweep-endsweep-hostname-starttime-currenttime
(defun generate-filename-v2 (&optional (ssweep 1) (csweep 0) (esweep (+ ssweep NUM-SWEEPS -1)))
  (format nil "~A-~A-T~3,'0d-V~6,'0d-~A-~A-~A-~A-~9,'0d-~9,'0d-~9,'0d-on-~A-start~A-curr~A" 
	  STOPOLOGY BCTYPE NUM-T N-INIT *k0* *k3* *eps* *alpha* ssweep csweep esweep 
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
(declaim (optimize (speed 3)
		   (compilation-speed 0)
		   (debug 0)
		   (safety 0)))
;; sl 1-simplex is (tslice (p0 p1))
;; tl 1-simplex is (type tmlo (p0 p1)) where type = 1
;; sl 2-simplex is (tslice (p0 p1 p2))
;; tl 2-simplex is (type tmlo (p0 p1 p2)) where type = 1,2
;; tl 3-simplex is (type tmlo tmhi (p0 p1 p2 p3) (n0 n1 n2 n3))
;; where type = 1,2,3,4 
;; nj = id of the 4sx that does not have pj

(defun make-2simplices (sxtype tmlo tmhi p0 p1 p2 p3)
  "makes all 2-simplices with the specified data iff they don't already exist"
  (cond ((= 1 sxtype) ; (p0 | p1 p2 p3)
	 (setf (gethash `(1 ,tmlo (,p0 ,p1 ,p2)) *TL2SIMPLEX->ID*) 0)
	 (setf (gethash `(1 ,tmlo (,p0 ,p1 ,p3)) *TL2SIMPLEX->ID*) 0)
	 (setf (gethash `(1 ,tmlo (,p0 ,p2 ,p3)) *TL2SIMPLEX->ID*) 0)
	 (setf (gethash `(,tmhi (,p1 ,p2 ,p3)) *SL2SIMPLEX->ID*) 0))
	((= 2 sxtype) ; (p0 p1 | p2 p3)
	 (setf (gethash `(2 ,tmlo (,p0 ,p1 ,p2)) *TL2SIMPLEX->ID*) 0)
	 (setf (gethash `(2 ,tmlo (,p0 ,p1 ,p3)) *TL2SIMPLEX->ID*) 0)
	 (setf (gethash `(1 ,tmlo (,p0 ,p2 ,p3)) *TL2SIMPLEX->ID*) 0)
	 (setf (gethash `(1 ,tmlo (,p1 ,p2 ,p3)) *TL2SIMPLEX->ID*) 0))
	((= 3 sxtype) ; (p0 p1 p2 | p3)
	 (setf (gethash `(,tmlo (,p0 ,p1 ,p2)) *SL2SIMPLEX->ID*) 0)
	 (setf (gethash `(2 ,tmlo (,p0 ,p1 ,p3)) *TL2SIMPLEX->ID*) 0)
	 (setf (gethash `(2 ,tmlo (,p0 ,p2 ,p3)) *TL2SIMPLEX->ID*) 0)
	 (setf (gethash `(2 ,tmlo (,p1 ,p2 ,p3)) *TL2SIMPLEX->ID*) 0))))

(defun make-1simplices (sxtype tmlo tmhi p0 p1 p2 p3)
  "makes all 1-simplices with the specified data iff they don't already exist"
  (cond ((= 1 sxtype) ;; p0 | p1 p2 p3
	 (setf (gethash `(1 ,tmlo (,p0 ,p1)) *TL1SIMPLEX->ID*) 0)
	 (setf (gethash `(1 ,tmlo (,p0 ,p2)) *TL1SIMPLEX->ID*) 0)
	 (setf (gethash `(1 ,tmlo (,p0 ,p3)) *TL1SIMPLEX->ID*) 0)
	 (setf (gethash `(,tmhi (,p1 ,p2)) *SL1SIMPLEX->ID*) 0)
	 (setf (gethash `(,tmhi (,p1 ,p3)) *SL1SIMPLEX->ID*) 0)
	 (setf (gethash `(,tmhi (,p2 ,p3)) *SL1SIMPLEX->ID*) 0))
	((= 2 sxtype) ;; p0 p1 | p2 p3
	 (setf (gethash `(,tmlo (,p0 ,p1)) *SL1SIMPLEX->ID*) 0)
	 (setf (gethash `(1 ,tmlo (,p0 ,p2)) *TL1SIMPLEX->ID*) 0)
	 (setf (gethash `(1 ,tmlo (,p0 ,p3)) *TL1SIMPLEX->ID*) 0)
	 (setf (gethash `(1 ,tmlo (,p1 ,p2)) *TL1SIMPLEX->ID*) 0)
	 (setf (gethash `(1 ,tmlo (,p1 ,p3)) *TL1SIMPLEX->ID*) 0)
	 (setf (gethash `(,tmhi (,p2 ,p3)) *SL1SIMPLEX->ID*) 0))
	((= 3 sxtype) ;; p0 p1 p2 | p3
	 (setf (gethash `(,tmlo (,p0 ,p1)) *SL1SIMPLEX->ID*) 0)
	 (setf (gethash `(,tmlo (,p0 ,p2)) *SL1SIMPLEX->ID*) 0)
	 (setf (gethash `(1 ,tmlo (,p0 ,p3)) *TL1SIMPLEX->ID*) 0)
	 (setf (gethash `(,tmlo (,p1 ,p2)) *SL1SIMPLEX->ID*) 0)
	 (setf (gethash `(1 ,tmlo (,p1 ,p3)) *TL1SIMPLEX->ID*) 0)
	 (setf (gethash `(1 ,tmlo (,p2 ,p3)) *TL1SIMPLEX->ID*) 0))))

(defun remove-tl2simplex (tl2sx)
  (remhash tl2sx *TL2SIMPLEX->ID*))
(defun remove-sl2simplex (sl2sx)
  (remhash sl2sx *SL2SIMPLEX->ID*))
(defun remove-tl2simplices (tl2sxs)
  (dolist (tl2sx tl2sxs)
    (remove-tl2simplex tl2sx)))
(defun remove-sl2simplices (sl2sxs)
  (dolist (sl2sx sl2sxs)
    (remove-sl2simplex sl2sx)))

(defun remove-tl1simplex (tl1sx)
  (remhash tl1sx *TL1SIMPLEX->ID*))
(defun remove-sl1simplex (sl1sx)
  (remhash sl1sx *SL1SIMPLEX->ID*))
(defun remove-tl1simplices (tl1sxs)
  (dolist (tl1sx tl1sxs)
    (remove-tl1simplex tl1sx)))
(defun remove-sl1simplices (sl1sxs)
  (dolist (sl1sx sl1sxs)
    (remove-sl1simplex sl1sx)))

(defun make-3simplex (sxtype tmlo tmhi p0 p1 p2 p3)
  (let ((sx3id (next-3simplex-id)))
    (setf (gethash sx3id *ID->3SIMPLEX*)
	  (list sxtype tmlo tmhi (list p0 p1 p2 p3) (list 0 0 0 0)))
    (make-2simplices sxtype tmlo tmhi p0 p1 p2 p3)
    (make-1simplices sxtype tmlo tmhi p0 p1 p2 p3)
    sx3id))

;; same as above except the points are "packed" into a list
(defun make-3simplex-v2 (sxtype tmlo tmhi pts)
    (let ((p0 (nth 0 pts))
	(p1 (nth 1 pts))
	(p2 (nth 2 pts))
	(p3 (nth 3 pts)))
    (make-3simplex sxtype tmlo tmhi p0 p1 p2 p3)))

;; this version is used only during initialization. If periodic b.c. are 
;; specified, it adjusts the points on the final time slice, since the t=T 
;; slice is identified with t=0 slice.
(defun make-3simplex-v3 (sxtype tmlo tmhitmp p0tmp p1tmp p2tmp p3tmp)
  (let ((p0 p0tmp) (p1 p1tmp) (p2 p2tmp) (p3 p3tmp) (tmhi tmhitmp))
    (when (and (string= BCTYPE "PERIODIC") (= NUM-T tmhitmp))
      (setf tmhi 0)
      (cond ((= 1 sxtype)
	     (decf p1 (* N0-PER-SLICE NUM-T)) 
	     (decf p2 (* N0-PER-SLICE NUM-T)) 
	     (decf p3 (* N0-PER-SLICE NUM-T)))
	    ((= 2 sxtype)
	     (decf p2 (* N0-PER-SLICE NUM-T)) 
	     (decf p3 (* N0-PER-SLICE NUM-T)))
	    ((= 3 sxtype)
	     (decf p3 (* N0-PER-SLICE NUM-T)))))
    (make-3simplex sxtype tmlo tmhi p0 p1 p2 p3)))

;; this version is used for loading the simplex data from file
(defun make-3simplex-v4 (sxtype tmlo tmhi p0 p1 p2 p3 n0 n1 n2 n3 sx3id)
  (setf (gethash sx3id *ID->3SIMPLEX*)
	(list sxtype tmlo tmhi (list p0 p1 p2 p3) (list n0 n1 n2 n3)))
  (make-2simplices sxtype tmlo tmhi p0 p1 p2 p3)
  (make-1simplices sxtype tmlo tmhi p0 p1 p2 p3))

;; all the simplex data, not just the points, is packed into a list
(defun make-3simplex-v5 (simplex-data)
  (let ((sxtype (first simplex-data))
	(tmlo (second simplex-data))
	(tmhi (third simplex-data))
	(pts (fourth simplex-data)))
    (make-3simplex-v2 sxtype tmlo tmhi pts)))

;; simplex-data-list = ((typ tmlo tmhi (p0 p1 p2 p3))...)
;; the ids of the simplices are returned
(defun make-3simplices-in-bulk (simplex-data-list)
  (let ((3sxids nil))
    (dolist (simplex-data simplex-data-list)
      (push (make-3simplex-v5 simplex-data) 3sxids))
    3sxids))

(defmacro 3sx-type (sx) `(first ,sx))
(defmacro 3sx-tmlo (sx) `(second ,sx))
(defmacro 3sx-tmhi (sx) `(third ,sx))
(defmacro 3sx-points (sx) `(fourth ,sx))
(defmacro 3sx-sx3ids (sx) `(fifth ,sx))
(defmacro 3sx-lopts (sx) `(subseq (3sx-points ,sx) 0 (3sx-type ,sx)))
(defmacro 3sx-hipts (sx) `(subseq (3sx-points ,sx) (3sx-type ,sx)))
(defmacro nth-point (sx n) `(nth ,n (3sx-points ,sx)))
(defmacro nth-neighbor (sx n) `(nth ,n (3sx-sx3ids ,sx)))

(defmacro get-3simplex (sxid)
  `(gethash ,sxid *ID->3SIMPLEX*))
(defmacro remove-3simplex (3sxid)
  `(progn
     (remhash ,3sxid *ID->3SIMPLEX*)
     (recycle-3simplex-id ,3sxid)))
(defmacro remove-3simplices (3sxids)
  `(dolist (3sxid ,3sxids)
     (remhash 3sxid *ID->3SIMPLEX*)
     (recycle-3simplex-id 3sxid)))


(defun show-id->3simplex-store ()
  (maphash #'(lambda (3sxid 3sx) 
	       (cond ((= 1 (3sx-type 3sx))
		      (format t "[~A] (~A ~A ~A (~A|~A ~A ~A) (~A ~A ~A ~A))~%"
			      3sxid (3sx-type 3sx) (3sx-tmlo 3sx) (3sx-tmhi 3sx)
			      (nth-point 3sx 0) (nth-point 3sx 1) 
			      (nth-point 3sx 2) (nth-point 3sx 3)
			      (nth-neighbor 3sx 0) (nth-neighbor 3sx 1) 
			      (nth-neighbor 3sx 2) (nth-neighbor 3sx 3)))
		     ((= 2 (3sx-type 3sx))
		      (format t "[~A] (~A ~A ~A (~A ~A|~A ~A) (~A ~A ~A ~A))~%"
			      3sxid (3sx-type 3sx) (3sx-tmlo 3sx) (3sx-tmhi 3sx)
			      (nth-point 3sx 0) (nth-point 3sx 1) 
			      (nth-point 3sx 2) (nth-point 3sx 3)
			      (nth-neighbor 3sx 0) (nth-neighbor 3sx 1) 
			      (nth-neighbor 3sx 2) (nth-neighbor 3sx 3)))
		     ((= 3 (3sx-type 3sx))
		      (format t "[~A] (~A ~A ~A (~A ~A ~A|~A) (~A ~A ~A ~A))~%"
			      3sxid (3sx-type 3sx) (3sx-tmlo 3sx) (3sx-tmhi 3sx)
			      (nth-point 3sx 0) (nth-point 3sx 1) 
			      (nth-point 3sx 2) (nth-point 3sx 3)
			      (nth-neighbor 3sx 0) (nth-neighbor 3sx 1) 
			      (nth-neighbor 3sx 2) (nth-neighbor 3sx 3)))))
	   *ID->3SIMPLEX*))

(defun show-tl2simplex-store ()
  (let ((count 1))
    (maphash #'(lambda (tl2sx id)
		 (format t "[~A] ~A ~A~%" count tl2sx id)
		 (incf count))
	     *TL2SIMPLEX->ID*)))
(defun show-sl2simplex-store ()
  (let ((count 1))
    (maphash #'(lambda (sl2sx id)
		 (format t "[~A] ~A ~A~%" count sl2sx id)
		 (incf count))
	     *SL2SIMPLEX->ID*)))

(defun show-tl1simplex-store ()
  (let ((count 1))
    (maphash #'(lambda (tl1sx id)
		 (format t "[~A] ~A ~A~%" count tl1sx id)
		 (incf count))
	     *TL1SIMPLEX->ID*)))
(defun show-sl1simplex-store ()
  (let ((count 1))
    (maphash #'(lambda (sl1sx id)
		 (format t "[~A] ~A ~A~%" count sl1sx id)
		 (incf count))
	     *SL1SIMPLEX->ID*)))

(defun connect-3simplices (sx1id sx2id)
  (let ((sx1 nil) (sx2 nil))
    (when (and (setf sx1 (get-3simplex sx1id))
	       (setf sx2 (get-3simplex sx2id)))
      (let ((tri (intersection (3sx-points sx1) (3sx-points sx2))))
	(when (= 3 (length tri))
	  (let* ((pts1 (3sx-points sx1))
		 (pts2 (3sx-points sx2))
		 (pos1 (position (first (set-difference pts1 tri)) pts1))
		 (pos2 (position (first (set-difference pts2 tri)) pts2)))
	    (setf (nth pos1 (3sx-sx3ids sx1)) sx2id 
		  (nth pos2 (3sx-sx3ids sx2)) sx1id)))))))

(defun connect-3simplices-within-list (sx1ids)
  (for (n 0 (1- (length sx1ids)))
       (for (m (1+ n) (1- (length sx1ids)))
	    (connect-3simplices (nth n sx1ids) (nth m sx1ids)))))

(defun connect-3simplices-across-lists (sx1ids sx2ids)
  (dolist (sx1id sx1ids)
    (dolist (sx2id sx2ids)
      (connect-3simplices sx1id sx2id))))

(defmacro 3simplices-connected? (sxid1 sxid2)
  `(let ((sx1 nil) (sx2 nil))
     (and (setf sx1 (get-3simplex ,sxid1)) 
	  (setf sx2 (get-3simplex ,sxid2))
	  (find ,sxid1 (3sx-sx3ids sx2)) 
	  (find ,sxid2 (3sx-sx3ids sx1)))))

(defun get-simplices-in-sandwich (tlo thi)
  (let ((sxids '()))
    (maphash #'(lambda (id sx)
		 (when (and (= (3sx-tmlo sx) (mod tlo NUM-T)) 
			    (= (3sx-tmhi sx) (mod thi NUM-T)))
		   (push id sxids)))
	     *ID->3SIMPLEX*)
    sxids))

(defun get-simplices-in-sandwich-of-type (tlo thi typ)
  (let ((sxids '()))
    (maphash #'(lambda (id sx)
		 (when (and (= (3sx-tmlo sx) (mod tlo NUM-T)) 
			    (= (3sx-tmhi sx) (mod thi NUM-T))
			    (= (3sx-type sx) typ))
		   (push id sxids)))
	     *ID->3SIMPLEX*)
    sxids))

;; returns (1ids 2ids 3ids) where 1ids is a list of (1,3) ids in the sandwich
;;  etc.
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

(defun neighbors-of-type (sx sxtype)
  (let ((nbors nil)
	(nsx nil)
	(nids (3sx-sx3ids sx)))
    (for (n 0 3)
      (when (and (setf nsx (get-3simplex (nth n nids))) 
		 (= sxtype (3sx-type nsx)))
	(pushnew (nth n nids) nbors)))
    nbors))

(defun save-spacetime-to-file (outfile)
  (format outfile "~A ~A ~A ~A ~A ~A ~A ~A ~A ~A ~A ~A ~A ~A ~A ~A ~A~%" 
	  BCTYPE STOPOLOGY NUM-T N-INIT *LAST-USED-POINT* *LAST-USED-3SXID* 
	  N0 N1-SL N1-TL N2-SL N2-TL N3-TL-31 N3-TL-22 *eps* *k0* *k3* *alpha*)
  (maphash #'(lambda (k v)
	       (format outfile "~A ~A ~A ~A ~A ~A ~A ~A ~A ~A ~A ~A~%" 
		       (3sx-type v) (3sx-tmlo v) (3sx-tmhi v)
		       (nth-point v 0) (nth-point v 1) 
		       (nth-point v 2) (nth-point v 3)
		       (nth-neighbor v 0) (nth-neighbor v 1)
		       (nth-neighbor v 2) (nth-neighbor v 3) k))
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
      (set-k0-k3-alpha (nth 14 data) (nth 15 data) (nth 16 data)))))

(defun parse-simplex-data-line (line)
  (with-input-from-string (s line)
    (let ((data (loop
		   :for num := (read s nil nil)
		   :while num
		   :collect num)))
      (make-3simplex-v4 (nth 0 data) (nth 1 data) (nth 2 data) (nth 3 data) 
			(nth 4 data) (nth 5 data) (nth 6 data) (nth 7 data)
			(nth 8 data) (nth 9 data) (nth 10 data) (nth 11 data)))))

(defun load-spacetime-from-file (infile)
  (parse-parameters-line (read-line infile nil))
  (loop for line = (read-line infile nil)
     while line do (parse-simplex-data-line line)))

;; spatial 2-simplex is a spatial triangle; this information is needed for 
;; computing the spectral and hausdorff dimensions of the spatial slice, among 
;; other things.
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

	     ;; try-a->b methods returns the following list, IFF the move can be successfully 
;; made 
;; (new3sxids nbors old3sxids oldTL2sxs oldSL2sxs oldTL1sxs oldSL1sx fvector)

(defun 2plus1move (sxdata)
  (let ((new3sxids (make-3simplices-in-bulk (first sxdata))))
    (connect-3simplices-within-list new3sxids)
    (connect-3simplices-across-lists new3sxids (second sxdata))
    (remove-3simplices (third sxdata))
    (remove-tl2simplices (fourth sxdata))
    (remove-sl2simplices (fifth sxdata))
    (remove-tl1simplices (sixth sxdata))
    (remove-sl1simplices (seventh sxdata))
    (update-f-vector (eighth sxdata))))

;;---------------------------------------------------------------------------[0]
;; 5
;;-------------- t+1
;; 2 3 4
;;-------------- t
;; 1
;;-------------- t-1
;;
;; input : either a (1,3) or a (3,1) 3-simplex
;;
;; 1(1,3) + 1(3,1) -> 3(1,3) + 3(3,1)
;;
;; (1|2 3 4) + (2 3 4|5) 
;; ->
;; (1|2 3 6) + (2 3 6|5) +
;; (1|2 6 4) + (2 6 4|5) +
;; (1|6 3 4) + (6 3 4|5)
;;---------------------------------------------------------------------------[0]
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
				  
(defun try-2->6 (sxid)
;;  (setf CURRENT-MOVE-IDENTIFIER "try-2->6")
  (dolist (curr (2->6-subcomplex sxid))
    (let* ((sx13 (get-3simplex (first curr)))
	   (sx31 (get-3simplex (second curr)))
	   (nbors (set-difference 
		   (union (3sx-sx3ids sx13) (3sx-sx3ids sx31)) 
		   curr))
	   (lopt (nth-point sx13 0))
	   (hipt (nth-point sx31 3))
	   (bigtr (3sx-lopts sx31))
	   (13tmlo (3sx-tmlo sx13)) 
	   (13tmhi (3sx-tmhi sx13))
	   (31tmlo (3sx-tmlo sx31)) 
	   (31tmhi (3sx-tmhi sx31))
	   (newpt (next-pt))
	   (oldSL2sxs `((,31tmlo ,bigtr)))
	   (newsxdata 
	    `((1 ,13tmlo ,13tmhi (,lopt ,newpt ,@(circular-subseq bigtr 0 2)))
	      (1 ,13tmlo ,13tmhi (,lopt ,newpt ,@(circular-subseq bigtr 1 2)))
	      (1 ,13tmlo ,13tmhi (,lopt ,newpt ,@(circular-subseq bigtr 2 2)))
	      (3 ,31tmlo ,31tmhi (,@(circular-subseq bigtr 0 2) ,newpt ,hipt))
	      (3 ,31tmlo ,31tmhi (,@(circular-subseq bigtr 1 2) ,newpt ,hipt))
	      (3 ,31tmlo ,31tmhi (,@(circular-subseq bigtr 2 2) ,newpt ,hipt)))))
      (return-from try-2->6 
	(list newsxdata nbors curr 
	      nil oldSL2sxs
	      nil nil
	      DF26)))))
;;---------------------------------------------------------------------------[1]
;; 5
;;-------------- t+1
;; 2 3 4
;;-------------- t
;; 1
;;-------------- t-1
;;
;; input : either a (1,3) or a (3,1) 3-simplex
;;
;; 3(1,3) + 3(3,1) -> 1(1,3) + 1(3,1)
;;
;; (1|2 3 6) + (2 3 6|5) +
;; (1|2 6 4) + (2 6 4|5) +
;; (1|6 3 4) + (6 3 4|5)
;; ->
;; (1|2 3 4) + (2 3 4|5) 
;;---------------------------------------------------------------------------[1]
(defun 6->2-subcomplex (sxid)
  "returns a list of the form ((13id1 13id2 13id3 31id3 31id2 31id1)...)"
  (let ((sx nil)
	(subcmplx nil))
    (when (setf sx (get-3simplex sxid))
      (cond ((= 1 (3sx-type sx))
	     (let ((31id (nth 0 (3sx-sx3ids sx)))
		   (31sx nil))
	       (when (setf 31sx (get-3simplex 31id))
		 (let ((13nbors (neighbors-of-type sx 1)))
		   (unless (< (length 13nbors) 2)
		     (do-tuples/c (currid nextid) 13nbors
		       (let ((curr nil) (next nil))
			 (when (and (setf curr (get-3simplex currid)) 
				    (setf next (get-3simplex nextid))
				    (3simplices-connected? currid nextid)
				    (3simplices-connected? 
				     (nth 0 (3sx-sx3ids curr))
				     (nth 0 (3sx-sx3ids next)))
				    (3simplices-connected? 
				     (nth 0 (3sx-sx3ids curr)) 31id)
				    (3simplices-connected? 
				     (nth 0 (3sx-sx3ids next)) 31id))
			   (pushnew (list sxid currid nextid 
					  (nth 0 (3sx-sx3ids next)) 
					  (nth 0 (3sx-sx3ids curr))
					  31id)
				    subcmplx :test #'set-equal?)))))))))
	    ((= 3 (3sx-type sx))
	     (let ((13id (nth 3 (3sx-sx3ids sx)))
		   (13sx nil))
	       (when (setf 13sx (get-3simplex 13id))
		 (let ((31nbors (neighbors-of-type sx 3)))
		   (unless (< (length 31nbors) 2)
		     (do-tuples/c (currid nextid) 31nbors
		       (let ((curr nil) (next nil))
			 (when (and (setf curr (get-3simplex currid)) 
				    (setf next (get-3simplex nextid))
				    (3simplices-connected? currid nextid)
				    (3simplices-connected? 
				     (nth 3 (3sx-sx3ids curr))
				     (nth 3 (3sx-sx3ids next)))
				    (3simplices-connected? 
				     (nth 3 (3sx-sx3ids curr)) 13id)
				    (3simplices-connected? 
				     (nth 3 (3sx-sx3ids next)) 13id))
			   (pushnew (list 13id 
					  (nth 3 (3sx-sx3ids next)) 
					  (nth 3 (3sx-sx3ids curr))
					  currid nextid sxid) 
				    subcmplx :test #'set-equal?)))))))))))
    subcmplx))
	    
(defun try-6->2 (sxid)
  ;;  (setf CURRENT-MOVE-IDENTIFIER "try-8->2")
  (dolist (subcx (6->2-subcomplex sxid))
    (let* ((1id1 (first subcx)) 
	   (1id2 (second subcx)) 
	   (1id3 (third subcx)) 
	   (3id3 (fourth subcx)) 
	   (3id2 (fifth subcx)) 
	   (3id1 (sixth subcx)) 
	   (1sx1 (get-3simplex 1id1)) 
	   (1sx2 (get-3simplex 1id2)) 
	   (1sx3 (get-3simplex 1id3)) 
	   (3sx1 (get-3simplex 3id1)) 
	   (3sx2 (get-3simplex 3id2))
	   (3sx3 (get-3simplex 3id3)) 
	   (13tmlo (3sx-tmlo 1sx1)) 
	   (13tmhi (3sx-tmhi 1sx1))
	   (31tmlo (3sx-tmlo 3sx1)) 
	   (31tmhi (3sx-tmhi 3sx1))
	   (nbors (set-difference 
		   (unions (3sx-sx3ids 1sx1) (3sx-sx3ids 1sx2) 
			   (3sx-sx3ids 1sx3) (3sx-sx3ids 3sx1) 
			   (3sx-sx3ids 3sx2) (3sx-sx3ids 3sx3))
		   (list 0 1id1 1id2 1id3 3id1 3id2 3id3)))
	   (pt1 (3sx-lopts 1sx1)) ; (1)
	   (pt5 (3sx-hipts 3sx1)) ; (5)
	   (pt6 (intersections (3sx-hipts 1sx1) (3sx-hipts 1sx2)
			       (3sx-hipts 1sx3))) ; (6)
	   (pts234 (set-difference (unions (3sx-hipts 1sx1) (3sx-hipts 1sx2)
					   (3sx-hipts 1sx3))
				   pt6)) ; (2 3 4)
	   (oldSL2sxs `((,13tmhi (,@pt6 ,@(circular-subseq pts234 0 2)))
			(,13tmhi (,@pt6 ,@(circular-subseq pts234 1 2)))
			(,13tmhi (,@pt6 ,@(circular-subseq pts234 2 2)))))
	   (oldTL2sxs `((1 ,13tmlo (,@pt1 ,(first pts234) ,@pt6))
			(1 ,13tmlo (,@pt1 ,(second pts234) ,@pt6))
			(1 ,13tmlo (,@pt1 ,(third pts234) ,@pt6))
			(2 ,31tmlo (,(first pts234) ,@pt6 ,@pt5))
			(2 ,31tmlo (,(second pts234) ,@pt6 ,@pt5))
			(2 ,31tmlo (,(third pts234) ,@pt6 ,@pt5))))
	   (oldSL1sxs `((,13tmhi (,(first pts234) ,@pt6))
			(,13tmhi (,(second pts234) ,@pt6))
			(,13tmhi (,(third pts234) ,@pt6))))
	   (oldTL1sxs `((1 ,13tmlo (,@pt1 ,@pt6))
			(1 ,31tmlo (,@pt6 ,@pt5))))
	   (newsxdata nil))
      (unless (gethash `(,13tmhi ,pts234) *SL2SIMPLEX->ID*)
	(setf newsxdata 
	      `((1 ,13tmlo ,13tmhi (,@pt1 ,@pts234))
		(3 ,31tmlo ,31tmhi (,@pts234 ,@pt5))))
	(return-from try-6->2 
	  (list newsxdata nbors subcx
		oldTL2sxs oldSL2sxs
		oldTL1sxs oldSL1sxs
		DF62))))))

;;---------------------------------------------------------------------------[2]
;; 6
;;-------------- t+1
;; 2 3 4 5
;;-------------- t
;; 1
;;-------------- t-1
;;
;; input : either a (1,3) or a (3,1) 3-simplex
;;
;; 2(1,3) + 2(3,1) -> 2(1,3) + 2(3,1)
;;
;; (1|2 3 4) + (1|3 4 5) +
;; (2 3 4|6) + (3 4 5|6)
;; ->
;; (1|2 4 5) + (1|2 3 5) +
;; (2 4 5|6) + (2 3 5|6)
;;---------------------------------------------------------------------------[2]
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
		     (when (3simplices-connected? 
			    (nth 0 (3sx-sx3ids (get-3simplex 13nbor))) 31id)
		       (pushnew (list sxid 13nbor 
				      (nth 0 (3sx-sx3ids (get-3simplex 13nbor)))
				      31id)
				subcmplx :test #'set-equal?)))))))
	    ((= 3 (3sx-type sx))
	     (let ((13id (nth 3 (3sx-sx3ids sx))) (13sx nil))
	       (when (setf 13sx (get-3simplex 13id))
		 (let ((31nbors (neighbors-of-type sx 3)))
		   (dolist (31nbor 31nbors)
		     (when (3simplices-connected? 
			    (nth 3 (3sx-sx3ids (get-3simplex 31nbor))) 13id)
		       (pushnew (list 13id 
				      (nth 3 (3sx-sx3ids (get-3simplex 31nbor)))
				      31nbor sxid)
				subcmplx :test #'set-equal?)))))))))
    subcmplx))

(defun try-4->4 (sxid)
;;  (setf CURRENT-MOVE-IDENTIFIER "try-4->4")
  (dolist (subcx (4->4-subcomplex sxid))
    (let* ((1id1 (first subcx)) 
	   (1id2 (second subcx)) 
	   (3id2 (third subcx)) 
	   (3id1 (fourth subcx)) 
	   (1sx1 (get-3simplex 1id1)) 
	   (1sx2 (get-3simplex 1id2)) 
	   (3sx1 (get-3simplex 3id1)) 
	   (3sx2 (get-3simplex 3id2))
	   (13tmlo (3sx-tmlo 1sx1)) 
	   (13tmhi (3sx-tmhi 1sx1))
	   (31tmlo (3sx-tmlo 3sx1)) 
	   (31tmhi (3sx-tmhi 3sx1))
	   (nbors (set-difference 
		   (unions (3sx-sx3ids 1sx1) (3sx-sx3ids 1sx2) 
			   (3sx-sx3ids 3sx1) (3sx-sx3ids 3sx2))
		   (list 0 1id1 1id2 3id1 3id2)))
	   (pt1 (3sx-lopts 1sx1)) ; (1)
	   (pt6 (3sx-hipts 3sx1)) ; (6)
	   (pts25 (set-exclusive-or (3sx-hipts 1sx1) (3sx-hipts 1sx2)));(2 5)
	   (pts34 (intersection (3sx-lopts 3sx1) (3sx-lopts 3sx2)));(3 4)
	   (oldTL2sxs `((1 ,13tmlo (,@pt1 ,@pts34))
			(2 ,31tmlo (,@pts34 ,@pt6))))
	   (oldSL2sxs `((,13tmhi (,(first pts25) ,@pts34))
			(,31tmlo (,@pts34 ,(second pts25)))))
	   (oldSL1sxs `((,13tmhi ,pts34)))
	   (newsxdata nil))
      (unless (gethash `(,31tmlo ,pts25) *SL1SIMPLEX->ID*)
	(setf newsxdata 
	      `((1 ,13tmlo ,13tmhi (,@pt1 ,@pts25 ,(first pts34)))
		(1 ,13tmlo ,13tmhi (,@pt1 ,@pts25 ,(second pts34)))
		(3 ,31tmlo ,31tmhi (,@pts25 ,(first pts34) ,@pt6))
		(3 ,31tmlo ,31tmhi (,@pts25 ,(second pts34) ,@pt6))))
	(return-from try-4->4 
	  (list newsxdata nbors subcx 
		oldTL2sxs oldSL2sxs
		nil oldSL1sxs
		DF44))))))
;;---------------------------------------------------------------------------[3]
;; 2 3 4
;;-------------- t+1
;; 1 5
;;-------------- t
;;
;; input : either a (1,3) or a (3,1) or a (2,2) 3-simplex
;;
;; 1(1,3) + 1(2,2) -> 1(1,3) + 2(2,2)
;;
;; (1|2 3 4) + (1 5|3 4) -> (5|2 3 4) + (1 5|2 3) + (1 5|2 4)
;; or
;; (2 3 4|1) + (3 4|1 5) -> (2 3 4|5) + (2 3|1 5) + (2 4|1 5)
;;---------------------------------------------------------------------------[3]
(defun 2->3-subcomplex (sxid)
  "returns a list of the form ((1or3 13or31id 22id)...) where the first number 1 or 3 tells us about the
type of the simplex participating in the move"
  (let ((sx nil)
	(subcmplx nil))
    (when (setf sx (get-3simplex sxid))
      (cond ((or (= 1 (3sx-type sx)) (= 3 (3sx-type sx)))
	     (let ((22nbors (neighbors-of-type sx 2)))
	       (dolist (22nbor 22nbors)
		 (pushnew (list (3sx-type sx) sxid 22nbor) 
			  subcmplx :test #'set-equal?))))
	    ((= 2 (3sx-type sx))
	     (let ((13nbors (append (neighbors-of-type sx 1))))
	       (dolist (13nbor 13nbors)
		 (pushnew (list 1 13nbor sxid) subcmplx :test #'set-equal?)))
	     (let ((31nbors (append (neighbors-of-type sx 3))))
	       (dolist (31nbor 31nbors)
		 (pushnew (list 3 31nbor sxid) 
			  subcmplx :test #'set-equal?))))))
    subcmplx))

(defun 2->3-move-internal-12 (13id 22id)
;;  (setf CURRENT-MOVE-IDENTIFIER "2->3-move-internal-12")
  (let ((13sx nil) (22sx nil))
    (when (and (setf 13sx (get-3simplex 13id)) 
	       (setf 22sx (get-3simplex 22id)))
      (let* ((pts234 (3sx-hipts 13sx));(2 3 4)
	     (pts34 (3sx-hipts 22sx));(3 4)
	     (pts15 (3sx-lopts 22sx));(1 5)
	     (pt1 (3sx-lopts 13sx));(1)
	     (pt5 (set-difference pts15 pt1));(5)
	     (pt2 (set-difference pts234 pts34));(2)
	     (nbors (set-difference 
		     (union (3sx-sx3ids 13sx) (3sx-sx3ids 22sx)) 
		     (list 0 13id 22id)))
	     (tmlo (3sx-tmlo 22sx))
	     (tmhi (3sx-tmhi 22sx))
	     (oldTL2sxs `((1 ,tmlo (,@pt1 ,@pts34))))
	     (newsxdata nil))
	(unless (gethash `(1 ,tmlo (,@pt5 ,@pt2)) *TL1SIMPLEX->ID*)
	  (setf newsxdata 
		`((1 ,tmlo ,tmhi (,@pt5 ,@pts234))
		  (2 ,tmlo ,tmhi (,@pt1 ,@pt5 ,@pt2 ,(first pts34)))
		  (2 ,tmlo ,tmhi (,@pt1 ,@pt5 ,@pt2 ,(second pts34)))))
	  (return-from 2->3-move-internal-12 
	    (list newsxdata nbors (list 13id 22id) 
		  oldTL2sxs nil
		  nil nil
		  DF23)))))))

(defun 2->3-move-internal-32 (31id 22id)
;;  (setf CURRENT-MOVE-IDENTIFIER "2->3-move-internal-32")
  (let ((31sx nil) (22sx nil))
    (when (and (setf 31sx (get-3simplex 31id)) 
	       (setf 22sx (get-3simplex 22id)))
      (let* ((pts234 (3sx-lopts 31sx));(2 3 4)
	     (pts34 (3sx-lopts 22sx));(3 4)
	     (pts15 (3sx-hipts 22sx));(1 5)
	     (pt1 (3sx-hipts 31sx));(1)
	     (pt5 (set-difference pts15 pt1));(5)
	     (pt2 (set-difference pts234 pts34));(2)
	     (nbors (set-difference 
		     (union (3sx-sx3ids 31sx) (3sx-sx3ids 22sx)) 
		     (list 0 31id 22id)))
	     (tmlo (3sx-tmlo 22sx))
	     (tmhi (3sx-tmhi 22sx))
	     (oldTL2sxs `((2 ,tmlo (,@pts34 ,@pt1))))
	     (newsxdata nil))
	(unless (gethash `(1 ,tmlo (,@pt5 ,@pt2)) *TL1SIMPLEX->ID*)
	  (setf newsxdata 
		`((3 ,tmlo ,tmhi (,@pts234 ,@pt5))
		  (2 ,tmlo ,tmhi (,@pt2 ,(first pts34) ,@pt1 ,@pt5))
		  (2 ,tmlo ,tmhi (,@pt2 ,(second pts34) ,@pt1 ,@pt5))))
	  (return-from 2->3-move-internal-32 
	    (list newsxdata nbors (list 31id 22id) 
		  oldTL2sxs nil
		  nil nil
		  DF23)))))))

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
;;---------------------------------------------------------------------------[4]
;; 2 3 4
;;-------------- t+1
;; 1 5
;;-------------- t
;;
;; input : either a (1,3) or a (3,1) or a (2,2) 3-simplex
;;
;; 1(1,3) + 2(2,2) -> 1(1,3) + 1(2,2)
;;
;; (5|2 3 4) + (1 5|2 3) + (1 5|2 4) -> (1|2 3 4) + (1 5|3 4)
;; or
;; (2 3 4|5) + (2 3|1 5) + (2 4|1 5) -> (2 3 4|1) + (3 4|1 5)
;;---------------------------------------------------------------------------[4]
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
			   (pushnew (list (3sx-type sx) sxid 22nbor 
					  22nborof22nbor) 
				    subcmplx :test #'set-equal?)))))))))
	    ((= 2 (3sx-type sx))
	     (let ((22nbors (neighbors-of-type sx 2)))
	       (dolist (22nbor 22nbors)
		 (let ((22sx (get-3simplex 22nbor)))
		   (when 22sx
		     (let ((13nborsof22nbor (neighbors-of-type 22sx 1)))
		       (dolist (13nborof22nbor 13nborsof22nbor)
			 (when (3simplices-connected? 13nborof22nbor sxid)
			   (pushnew (list 1 13nborof22nbor sxid 22nbor) 
				    subcmplx :test #'set-equal?))))
		     (let ((31nborsof22nbor (neighbors-of-type 22sx 3)))
		       (dolist (31nborof22nbor 31nborsof22nbor)
			 (when (3simplices-connected? 31nborof22nbor sxid)
			   (pushnew (list 3 31nborof22nbor sxid 22nbor) 
				    subcmplx :test #'set-equal?)))))))))))
    subcmplx))

(defun 3->2-move-internal-122 (13id 22id1 22id2)
  "the (3,2) move performed on a (1,3) simplex attached to two (2,2) simplices"
  (let ((13sx nil) (22sx1 nil) (22sx2 nil))
    (when (and (setf 13sx (get-3simplex 13id)) 
	       (setf 22sx1 (get-3simplex 22id1))
	       (setf 22sx2 (get-3simplex 22id2)))
      (let* ((pts234 (3sx-hipts 13sx));(2 3 4)
	     (pts15 (3sx-lopts 22sx1));(1 5)
	     (pt5 (3sx-lopts 13sx));(5)
	     (pt1 (set-difference pts15 pt5));(1)
	     (pts34 (set-exclusive-or 
		     (3sx-hipts 22sx1) (3sx-hipts 22sx2)));(3 4)
	     (pt2 (set-difference pts234 pts34));(2)
	     (nbors (set-difference 
		     (unions (3sx-sx3ids 13sx) (3sx-sx3ids 22sx1) 
			     (3sx-sx3ids 22sx2)) 
		     (list 0 13id 22id1 22id2)))
	     (tmlo (3sx-tmlo 22sx1))
	     (tmhi (3sx-tmhi 22sx1))
	     (oldTL2sxs `((1 ,tmlo (,@pt5 ,@pt2 ,(first pts34)))
			  (1 ,tmlo (,@pt5 ,@pt2 ,(second pts34)))
			  (2 ,tmlo (,@pts15 ,@pt2))))
	     (oldTL1sxs `((1 ,tmlo (,@pt5 ,@pt2))))
	     (newsxdata nil))
	(unless (gethash `(1 ,tmlo (,@pt1 ,@pts34)) *TL2SIMPLEX->ID*)
	  (setf newsxdata 
		`((1 ,tmlo ,tmhi (,@pt1 ,@pts234))
		  (2 ,tmlo ,tmhi (,@pts15 ,@pts34))))
	  (return-from 3->2-move-internal-122 
	    (list newsxdata nbors (list 13id 22id1 22id2)
		  oldTL2sxs nil
		  oldTL1sxs nil
		  DF32)))))))

(defun 3->2-move-internal-322 (31id 22id1 22id2)
  (let ((31sx nil) (22sx1 nil) (22sx2 nil))
    (when (and (setf 31sx (get-3simplex 31id)) 
	       (setf 22sx1 (get-3simplex 22id1))
	       (setf 22sx2 (get-3simplex 22id2)))
      (let* ((pts234 (3sx-lopts 31sx));(2 3 4)
	     (pts15 (3sx-hipts 22sx1));(1 5)
	     (pt5 (3sx-hipts 31sx));(5)
	     (pt1 (set-difference pts15 pt5));(1)
	     (pts34 (set-exclusive-or 
		     (3sx-lopts 22sx1) (3sx-lopts 22sx2)));(3 4)
	     (pt2 (set-difference pts234 pts34))
	     (nbors (set-difference 
		     (unions (3sx-sx3ids 31sx) (3sx-sx3ids 22sx1) 
			     (3sx-sx3ids 22sx2)) 
		     (list 0 31id 22id1 22id2)))
	     (tmlo (3sx-tmlo 22sx1))
	     (tmhi (3sx-tmhi 22sx1))
	     (oldTL2sxs `((2 ,tmlo (,@pt2 ,(first pts34) ,@pt5))
			  (2 ,tmlo (,@pt2 ,(second pts34) ,@pt5))
			  (1 ,tmlo (,@pt2 ,@pts15))))
	     (oldTL1sxs `((1 ,tmlo (,@pt2 ,@pt5))))
	     (newsxdata nil))
	(unless (gethash `(2 ,tmlo (,@pts34 ,@pt1)) *TL2SIMPLEX->ID*)
	  (setf newsxdata 
		`((3 ,tmlo ,tmhi (,@pts234 ,@pt1))
		  (2 ,tmlo ,tmhi (,@pts34 ,@pts15))))
	  (return-from 3->2-move-internal-322 
	    (list newsxdata nbors (list 31id 22id1 22id2)
		  oldTL2sxs nil
		  oldTL1sxs nil
		  DF32)))))))

(defun try-3->2 (sxid)
  (let ((subcmplx (3->2-subcomplex sxid))
	(movedata nil))
    (unless (null subcmplx)
      (dolist (curr subcmplx)
	(cond ((= 1 (first curr))
	       (setf movedata (3->2-move-internal-122 
			       (second curr) (third curr) (fourth curr)))
	       (when movedata
		 (return-from try-3->2 movedata)))
	      ((= 3 (first curr))
	       (setf movedata (3->2-move-internal-322 
			       (second curr) (third curr) (fourth curr)))
	       (when movedata 
		 (return-from try-3->2 movedata))))))))
;...........................................................................................................
; cdt-2plus1-initialization.lisp
;...........................................................................................................
(defun initialize-S2-triangulation (num-time-intervals boundary-conditions)
  (setf NUM-T num-time-intervals)
  (setf BCTYPE (string-upcase boundary-conditions))
  (for (n 0 (1- (/ num-time-intervals 2)))
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

(defun initialize-T-slices-with-V-volume (&key 
					  num-time-slices 
					  target-volume
					  spatial-topology
					  boundary-conditions)
  (setf STOPOLOGY (string-upcase spatial-topology))

  (when (string= STOPOLOGY "S2")
    (initialize-S2-triangulation num-time-slices boundary-conditions))

  (when (string= STOPOLOGY "T2")
    (initialize-T2-triangulation num-time-slices boundary-conditions))

;;  (maphash #'(lambda (id sx)
;;	       (declare (ignore sx))
;;	       (update-subcomplexes-for-3simplex id))
;;	   *ID->3SIMPLEX*)

  (format t "initial count = ~A~%" (count-simplices-of-all-types))
  
  (loop named tv
     do
       (dolist (id23 (get-simplices-of-type 2))
	 (let ((movedata nil))
	   (when (setf movedata (try-2->3 id23))
	     (2plus1move movedata)))
	 (if (> (N3) target-volume)
	     (return-from tv)))

       ;; (4,4) moves to mix things up
       (dolist (id44 (get-simplices-of-type 1))
	 (let ((movedata nil))
	   (when (setf movedata (try-4->4 id44))
	     (2plus1move movedata))))

       (dolist (id26 (get-simplices-of-type 3))
	 (let ((movedata nil))
	   (when (setf movedata (try-2->6 id26))
	     (2plus1move movedata)))
	 (if (> (N3) target-volume)
	     (return-from tv)))

       ;; (4,4) moves to mix things up
       (dolist (id44 (get-simplices-of-type 3))
	 (let ((movedata nil))
	   (when (setf movedata (try-4->4 id44))
	     (2plus1move movedata))))

       (dolist (id23 (get-simplices-of-type 2))
	 (let ((movedata nil))
	   (when (setf movedata (try-2->3 id23))
	     (2plus1move movedata)))
	 (if (> (N3) target-volume)
	     (return-from tv))))

  (format t "final count = ~A~%" (count-simplices-of-all-types))

  (setf N-INIT (N3)))
; cdt-2plus1-montecarlo.lisp

(defun try-move (sxid mtype)
  (ecase mtype
    (0 (try-2->6 sxid))
    (1 (try-2->3 sxid))
    (2 (try-4->4 sxid))
    (3 (try-3->2 sxid))
    (4 (try-6->2 sxid))))

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

(defun accept-move? (mtype)
  (let ((delta-action 0.0)
	(delta-damping 0.0))
    (cond ((= 0 mtype) ;; 2->6 move
	   (setf delta-damping (- (damping (+ (N3) 4)) (damping (N3))))
	   (setf delta-action 
		 (- (action (+ N1-SL 3) (+ N1-TL 2) 
			    (+ N3-TL-31 4) (+ N3-TL-22 0)) 
		    (action N1-SL N1-TL N3-TL-31 N3-TL-22)))) 
	  ((= 1 mtype) ;; 2->3 move
	   (setf delta-damping (- (damping (+ (N3) 1)) (damping (N3))))
	   (setf delta-action (- (action (+ N1-SL 0) (+ N1-TL 1) (+ N3-TL-31 0) (+ N3-TL-22 1)) 
				 (action N1-SL N1-TL N3-TL-31 N3-TL-22))))
	  ((= 2 mtype) ;; 4->4 move
	   (setf delta-damping (- (damping (+ (N3) 0)) (damping (N3))))
	   (setf delta-action (- (action (+ N1-SL 0) (+ N1-TL 0) (+ N3-TL-31 0) (+ N3-TL-22 0)) 
				 (action N1-SL N1-TL N3-TL-31 N3-TL-22))))
	  ((= 3 mtype) ;; 3->2 move
	   (setf delta-damping (- (damping (+ (N3) -1)) (damping (N3))))
	   (setf delta-action (- (action (+ N1-SL 0) (+ N1-TL -1) (+ N3-TL-31 0) (+ N3-TL-22 -1)) 
				 (action N1-SL N1-TL N3-TL-31 N3-TL-22))))
	  ((= 4 mtype) ;; 6->2 move
	   (setf delta-damping (- (damping (+ (N3) -4)) (damping (N3))))
	   (setf delta-action (- (action (+ N1-SL -3) (+ N1-TL -2) (+ N3-TL-31 -4) (+ N3-TL-22 0)) 
				 (action N1-SL N1-TL N3-TL-31 N3-TL-22)))))
    (< (random 1.0) (* (exp (realpart (* *i* delta-action)))
		       (exp (* -1.0 delta-damping))))))

;; a sweep is defined as N-INIT number of attempted moves
(defun sweep ()
  (let ((num-attempted 0))
    (while (< num-attempted N-INIT)
      (let* ((sxid (random *LAST-USED-3SXID*))
	     (mtype (select-move))
	     (movedata (try-move sxid mtype)))
	(while (null movedata) 
	  (setf sxid (random *LAST-USED-3SXID*)
		mtype (select-move) 
		movedata (try-move sxid mtype)))
	(incf num-attempted) ;; number-of-attempted-moves-counter for this sweep
	(incf (nth mtype ATTEMPTED-MOVES)) ;; number of moves of mtype that have been attempted
	(when (accept-move? mtype)
	  (incf (nth mtype SUCCESSFUL-MOVES)) ;; number of moves of mtype that have succeeded
	  (2plus1move movedata))))))

;; following is to be used for tuning the k0 and k3 parameters
(defun generate-data-console (&optional (start-sweep 1))
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

;; generate-data-v3 is similar to generate-data-v2 except it also creates an 
;; additional data file every 
;; SAVE-EVERY-N-SWEEPS that contains the spatial 2-simplex information for 
;; each spatial slice.
(defun generate-data-v3 (&optional (start-sweep 1))
  (setf SIM-START-TIME (cdt-now-str))
  (let ((end-sweep (+ start-sweep NUM-SWEEPS -1)))
    (when (= 1 start-sweep) ;; save the initial spacetime contents if this is a brand new run
      (with-open-file 
	  (datafile 
	   (concatenate 'string (generate-filename-v2 start-sweep 0) 3SXEXT)
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
  (let ((moviefilestr 
	 (concatenate 'string (generate-filename start-sweep) MOVEXT))
	(trackfilestr 
	 (concatenate 'string (generate-filename start-sweep) PRGEXT))
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
  (let ((end-sweep (+ start-sweep NUM-SWEEPS -1)))
       (for (ns start-sweep end-sweep)
	   (sweep)
	   (when (= 0 (mod ns 10))
	     (format t "~A/~A " ns end-sweep)
	     (for (ts 0 (1- NUM-T))
		  (format t "~A " (count-simplices-in-sandwich ts (1+ ts))))
	     (format t "~A ~A~%" (count-simplices-of-all-types) (accept-ratios))))))

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