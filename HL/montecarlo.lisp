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