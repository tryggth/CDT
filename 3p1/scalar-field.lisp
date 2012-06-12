;; scalar-field.lisp
;(load "../utilities.lisp")
;(load "../queue.lisp")
;(load "globals.lisp")
;(load "simplex.lisp")

(defconstant a 1.0)
(defconstant m 0.001)
(defconstant m^2/2 (/ (* m m) 2.0))
(defconstant 1/2a^2 (/ 1.0 (* 2.0 a a)))

(defparameter NUM-ATTEMPTED 0)
(defparameter NUM-ACCEPTED 0)
(defparameter SFALPHA 0.0)
(defparameter 4SX14VOL 0.0)
(defparameter 4SX23VOL 0.0)

(defun compute-alpha-and-4sxvols ()
  "computes the value of alpha 4sxvolumes based on KAPPA-0, DELTA, KAPPA-4"
  (setf SFALPHA 1.0) ;; this needs to be fixed later
  (setf 4SX14VOL (/ (sqrt (+ (* 8.0 SFALPHA) 3.0)) 96.0)) 
  (setf 4SX23VOL (/ (sqrt (+ (* 12.0 SFALPHA) 7.0)) 96.0)))


(defun rndval ()
  "returns a uniform random value in (-1.0,1.0)"
  (/ (- (random 2000000) 1000000) 1000000.0))

(defparameter *VERTEX->NUM4SX* (make-hash-table)
  "the key is a vertex p and the value is (n1 n2 n3 n4) where ni is the number of 
(i,5-i) 4-simplices that p belongs to")

(defparameter *VERTEX->VERTEX* (make-hash-table))
(defparameter *VERTEX->PHI* (make-hash-table))

(defun link-hash (link) (sxhash (sort (copy-list link) #'<)))
(defun link-equal (link1 link2) (set-equal? link1 link2))
(sb-ext:define-hash-table-test link-equal link-hash)

(defparameter *LINK->NUM4SX* (make-hash-table :test 'link-equal)
  "the key is a link (p1 p2) and the value is (n1 n2 n3 n4) where ni is the number of 
(i,5-i) 4-simplices that (p1 p2) belongs to")

(defparameter *VERTEX->COLOR* (make-hash-table)) ; 0=WHITE, 1=GRAY, 2=BLACK
(defparameter *VERTEX->DIST* (make-hash-table)) ; -1 MEANS INFINITY
(defparameter *VERTEX->PRED* (make-hash-table))
(defparameter *DIST->VERTICES* (make-hash-table))

(defun breadth-first-search (s)
  "does a breadth first search starting from source vertex s"
  ; clear the bfs tables
  (clrhash *VERTEX->COLOR*)
  (clrhash *VERTEX->DIST*)
  (clrhash *VERTEX->PRED*)
  (clrhash *DIST->VERTICES*)

  ; populate the bfs tables
  (maphash #'(lambda (k v)
	       (declare (ignore v))
	       (setf (gethash k *VERTEX->COLOR*) 0)
	       (setf (gethash k *VERTEX->DIST*) 'INFINITY)
	       (setf (gethash k *VERTEX->PRED*) nil))
	   *VERTEX->NUM4SX*)

  ; update the bfs tables for the source vertex
  (setf (gethash s *VERTEX->COLOR*) 1) ; color the source vertex gray
  (setf (gethash s *VERTEX->DIST*) 0)
  (setf (gethash s *VERTEX->PRED*) nil)

  (let ((Q (make-instance 'queue))) ; a fifo queue to store the gray vertices
    (enqueue Q s)
    (while (not (queue-empty-p Q))
      (let ((u (dequeue Q)))
	(dolist (v (gethash u *VERTEX->VERTEX*)); for all vs that u is connected to
	  (when (= 0 (gethash v *VERTEX->COLOR*)) ; found a white vertex
	    (setf (gethash v *VERTEX->COLOR*) 1) ; color this vertex gray
	    (setf (gethash v *VERTEX->DIST*) (+ 1 (gethash u *VERTEX->DIST*)))
	    (setf (gethash v *VERTEX->PRED*) u)
	    (enqueue Q v)))
  	(setf (gethash u *VERTEX->COLOR*) 2)))) ; color vertex u black
  (maphash #'(lambda (k v)
	       (pushnew k (gethash v *DIST->VERTICES*)))
	   *VERTEX->DIST*))

(defun compute-phi-phi (x ys)
  "ys=(y1 ... yN); computes (phi(x)*phi(y1) + ... + phi(x)*phi(yN))/N"
  (let ((sum 0.0)
	(phix (gethash x *VERTEX->PHI*)))
    (dolist (y ys)
      (incf sum (* phix (gethash y *VERTEX->PHI*))))
    (/ sum (length ys))))

(defun print-phi-phi-expectation (s outfile)
  "s is the source vertex"
  (maphash #'(lambda (k v)
	       (format outfile "~A,~A~%"
		       k (compute-phi-phi s v)))
	   *DIST->VERTICES*))

(defun generate-data (infilename numthermal numsweeps save-every-n-sweeps)
  "numthermal is the number of sweeps for thermalization and numsweeps is the
total number of sweeps."
  ; clear all existing tables
  (clrhash *ID->4SIMPLEX*)
  (clrhash *TL3SIMPLEX->ID*)
  (clrhash *sL3SIMPLEX->ID*)
  (clrhash *TL2SIMPLEX->ID*)
  (clrhash *sL2SIMPLEX->ID*)
  (clrhash *TL1SIMPLEX->ID*)
  (clrhash *sL1SIMPLEX->ID*)
  (clrhash *VERTEX->NUM4SX*)
  (clrhash *VERTEX->VERTEX*)
  (clrhash *VERTEX->PHI*)
  (clrhash *LINK->NUM4SX*)
  (clrhash *VERTEX->COLOR*)
  (clrhash *VERTEX->DIST*)
  (clrhash *VERTEX->PRED*)
  (clrhash *DIST->VERTICES*)

  ;; load the spacetime data
  (with-open-file (infile infilename :direction :input)
    (load-spacetime-from-file infile))

  ;; compute the value of alpha for this spacetime
  (compute-alpha-and-4sxvols)

  ;; populate scalar field tables
  (populate-scalar-field-tables)

  ;; thermalize
  (dotimes (ns numthermal)
    (sweep))

  (let* ((dotpos (position #\. infilename :from-end t))
	 (prefix (subseq infilename 0 dotpos)))
    ;; save the thermalized setup as our initial data point
    (let ((outfilename (concatenate 'string 
				    prefix 
				    "-sfs-" 
				    (format nil "~A" numthermal)
				    ".sf3p1"))
	  (source-vertex (select-random-site)))
      (breadth-first-search source-vertex)
      (with-open-file (outfile outfilename :direction :output)
	(print-phi-phi-expectation source-vertex outfile)))
    ;; continue sweeping, saving every 10000 sweeps
    (for (ns (1+ numthermal) numsweeps)
	 (sweep)
	 (when (= 0 (mod ns save-every-n-sweeps))
	   (let ((outfilename (concatenate 'string 
					   prefix 
					   "-sfs-" 
					   (format nil "~A" ns)
					   ".sf3p1"))
		 (source-vertex (select-random-site)))
	     (breadth-first-search source-vertex)
	     (with-open-file (outfile outfilename :direction :output)
	       (print-phi-phi-expectation source-vertex outfile)))))))

(defun select-random-site ()
  "selects a random site"
  (let ((marker (random (hash-table-count *VERTEX->NUM4SX*)))
	(counter 0))
    (maphash #'(lambda (k v)
	       (declare (ignore v))
	       (if (= counter marker)
		   (return-from select-random-site k)
		   (incf counter)))
	     *VERTEX->NUM4SX*)))

(defun vertex-volume (i)
  "computes the volume associated with site i. the factor of 0.2 below is because each 
4-simplex contributes only a fifth of its volume to each of its five vertices"
  (let* ((n4sxs (gethash i *VERTEX->NUM4SX*))
         (vol14p41 (* 4SX14VOL 0.2 (+ (first n4sxs) (fourth n4sxs))))
         (vol23p32 (* 4SX23VOL 0.2 (+ (second n4sxs) (third n4sxs)))))
    (+ vol14p41 vol23p32)))

(defun link-volume (ij)
  "computes the volume associated with link ij. the factor of 0.1 below is because each
4-simplex contributes only a tenth of its volume to each of its ten edges"
  (let* ((n4sxs (gethash ij *LINK->NUM4SX*))
         (vol14p41 (* 4SX14VOL 0.1 (+ (first n4sxs) (fourth n4sxs))))
         (vol23p32 (* 4SX23VOL 0.1 (+ (second n4sxs) (third n4sxs)))))
    (+ vol14p41 vol23p32)))

(defun delta-action (site newphi)
  "computes the change in the scalar field action if we change the field value
at site to newphi"
  (let* ((oldphi (gethash site *VERTEX->PHI*))
	 (newphisqr (* newphi newphi))
	 (oldphisqr (* oldphi oldphi))
	 (sitevol (vertex-volume site))
	 (links (mapcar #'(lambda (k) (list site k)) 
			(gethash site *VERTEX->VERTEX*)))
         (linkvols (mapcar #'link-volume links))
	 (phijs (mapcar #'(lambda (k) (gethash k *VERTEX->PHI*))
			(gethash site *VERTEX->VERTEX*)))
	 (term2 (* m^2/2 sitevol (- newphisqr oldphisqr)))
	 (jkterms (mapcar #'(lambda (v phij)
			      (* v (- newphisqr oldphisqr 
				      (* 2 phij (- newphi oldphi)))))
			 linkvols phijs))
	 (term1 (* 1/2a^2 (reduce #'+ jkterms))))
    (+ term1 term2)))

(defun accept-move? (site newphi)
  "if the move is accepted, site->phi = newphi else site->phi = oldphi"
  (let ((delact (delta-action site newphi)))
    (if (<= delact 0.0)
        T
	(if (< (random 1.0) (exp (* -1.0 delact)))
            T
            nil))))

(defun sweep ()
  (setf NUM-ATTEMPTED 0)
  (setf NUM-ACCEPTED 0)
  "perform a sweep number of moves"
  (for (num 1 (hash-table-count *VERTEX->NUM4SX*))
       (incf NUM-ATTEMPTED)
       (let ((site (select-random-site))
	     (newphi (rndval)))
	 (when (accept-move? site newphi)
           (incf NUM-ACCEPTED)
           (setf (gethash site *VERTEX->PHI*) newphi))))
  (* 100.0 (/ NUM-ACCEPTED NUM-ATTEMPTED)))

;--------1---------2---------3---------4---------5---------6---------7---------8

(defun parse-4simplex (4sx)
  (let ((p0 (nth-point 4sx 0))
	(p1 (nth-point 4sx 1))
	(p2 (nth-point 4sx 2))
	(p3 (nth-point 4sx 3))
	(p4 (nth-point 4sx 4))
        (typ (4sx-type 4sx)))
    ;
    (if (gethash p0 *VERTEX->NUM4SX*)
        (incf (nth (1- typ) (gethash p0 *VERTEX->NUM4SX*)))
        (setf (gethash p0 *VERTEX->NUM4SX*) (list 0 0 0 0)))
    (if (gethash p1 *VERTEX->NUM4SX*)
        (incf (nth (1- typ) (gethash p1 *VERTEX->NUM4SX*)))
        (setf (gethash p1 *VERTEX->NUM4SX*) (list 0 0 0 0)))
    (if (gethash p2 *VERTEX->NUM4SX*)
        (incf (nth (1- typ) (gethash p2 *VERTEX->NUM4SX*)))
        (setf (gethash p2 *VERTEX->NUM4SX*) (list 0 0 0 0)))
    (if (gethash p3 *VERTEX->NUM4SX*)
        (incf (nth (1- typ) (gethash p3 *VERTEX->NUM4SX*)))
        (setf (gethash p3 *VERTEX->NUM4SX*) (list 0 0 0 0)))
    (if (gethash p4 *VERTEX->NUM4SX*)
        (incf (nth (1- typ) (gethash p4 *VERTEX->NUM4SX*)))
        (setf (gethash p4 *VERTEX->NUM4SX*) (list 0 0 0 0)))

    ;
    (setf (gethash p0 *VERTEX->PHI*) (rndval))
    (setf (gethash p1 *VERTEX->PHI*) (rndval))
    (setf (gethash p2 *VERTEX->PHI*) (rndval))
    (setf (gethash p3 *VERTEX->PHI*) (rndval))
    (setf (gethash p4 *VERTEX->PHI*) (rndval))
    ;
    (pushnew p1 (gethash p0 *VERTEX->VERTEX*))
    (pushnew p2 (gethash p0 *VERTEX->VERTEX*))
    (pushnew p3 (gethash p0 *VERTEX->VERTEX*))
    (pushnew p4 (gethash p0 *VERTEX->VERTEX*))
    (pushnew p0 (gethash p1 *VERTEX->VERTEX*))
    (pushnew p2 (gethash p1 *VERTEX->VERTEX*))
    (pushnew p3 (gethash p1 *VERTEX->VERTEX*))
    (pushnew p4 (gethash p1 *VERTEX->VERTEX*))
    (pushnew p0 (gethash p2 *VERTEX->VERTEX*))
    (pushnew p1 (gethash p2 *VERTEX->VERTEX*))
    (pushnew p3 (gethash p2 *VERTEX->VERTEX*))
    (pushnew p4 (gethash p2 *VERTEX->VERTEX*))
    (pushnew p0 (gethash p3 *VERTEX->VERTEX*))
    (pushnew p1 (gethash p3 *VERTEX->VERTEX*))
    (pushnew p2 (gethash p3 *VERTEX->VERTEX*))
    (pushnew p4 (gethash p3 *VERTEX->VERTEX*))
    (pushnew p0 (gethash p4 *VERTEX->VERTEX*))
    (pushnew p1 (gethash p4 *VERTEX->VERTEX*))
    (pushnew p2 (gethash p4 *VERTEX->VERTEX*))
    (pushnew p3 (gethash p4 *VERTEX->VERTEX*))
    ;
    (if (gethash `(,p0 ,p1) *LINK->NUM4SX*)
        (incf (nth (1- typ) (gethash `(,p0 ,p1) *LINK->NUM4SX*)))
        (setf (gethash `(,p0 ,p1) *LINK->NUM4SX*) (list 0 0 0 0)))
    (if (gethash `(,p0 ,p2) *LINK->NUM4SX*)
        (incf (nth (1- typ) (gethash `(,p0 ,p2) *LINK->NUM4SX*)))
        (setf (gethash `(,p0 ,p2) *LINK->NUM4SX*) (list 0 0 0 0)))
    (if (gethash `(,p0 ,p3) *LINK->NUM4SX*)
        (incf (nth (1- typ) (gethash `(,p0 ,p3) *LINK->NUM4SX*)))
        (setf (gethash `(,p0 ,p3) *LINK->NUM4SX*) (list 0 0 0 0)))
    (if (gethash `(,p0 ,p4) *LINK->NUM4SX*)
        (incf (nth (1- typ) (gethash `(,p0 ,p4) *LINK->NUM4SX*)))
        (setf (gethash `(,p0 ,p4) *LINK->NUM4SX*) (list 0 0 0 0)))
    (if (gethash `(,p1 ,p2) *LINK->NUM4SX*)
        (incf (nth (1- typ) (gethash `(,p1 ,p2) *LINK->NUM4SX*)))
        (setf (gethash `(,p1 ,p2) *LINK->NUM4SX*) (list 0 0 0 0)))
    (if (gethash `(,p1 ,p3) *LINK->NUM4SX*)
        (incf (nth (1- typ) (gethash `(,p1 ,p3) *LINK->NUM4SX*)))
        (setf (gethash `(,p1 ,p3) *LINK->NUM4SX*) (list 0 0 0 0)))
    (if (gethash `(,p1 ,p4) *LINK->NUM4SX*)
        (incf (nth (1- typ) (gethash `(,p1 ,p4) *LINK->NUM4SX*)))
        (setf (gethash `(,p1 ,p4) *LINK->NUM4SX*) (list 0 0 0 0)))
    (if (gethash `(,p2 ,p3) *LINK->NUM4SX*)
        (incf (nth (1- typ) (gethash `(,p2 ,p3) *LINK->NUM4SX*)))
        (setf (gethash `(,p2 ,p3) *LINK->NUM4SX*) (list 0 0 0 0)))
    (if (gethash `(,p2 ,p4) *LINK->NUM4SX*)
        (incf (nth (1- typ) (gethash `(,p2 ,p4) *LINK->NUM4SX*)))
        (setf (gethash `(,p2 ,p4) *LINK->NUM4SX*) (list 0 0 0 0)))
    (if (gethash `(,p3 ,p4) *LINK->NUM4SX*)
        (incf (nth (1- typ) (gethash `(,p3 ,p4) *LINK->NUM4SX*)))
        (setf (gethash `(,p3 ,p4) *LINK->NUM4SX*) (list 0 0 0 0)))))

(defun populate-scalar-field-tables ()
  "assumes spacetime has been loaded"
  (maphash #'(lambda (k v)
	       (declare (ignore k))
	       (parse-4simplex v))
	   *ID->4SIMPLEX*))
