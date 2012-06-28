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

(defun rndval ()
  "returns a uniform random value in (-1.0,1.0)"
  (/ (- (random 2000000) 1000000) 1000000.0))

;-------------------------------------------------------------------------------
; *VERTEX->VERTEX
;-------------------------------------------------------------------------------
(defparameter *VERTEX->VERTEX* (make-hash-table)
  "the key is a vertex p and the value is a list of vertices that p is 
connected to")

;-------------------------------------------------------------------------------
; *VERTEX->PHI*
;-------------------------------------------------------------------------------
(defparameter *VERTEX->PHI* (make-hash-table)
  "the key is a vertex p and the value is the scalar field value at that
vertex")

;-------------------------------------------------------------------------------
; *VERTEX->NUM3SX*
;-------------------------------------------------------------------------------
(defparameter *VERTEX->NUM3SX* (make-hash-table)
  "the key is a vertex p and the value is the number of 3-simplices that p 
belongs to")

;-------------------------------------------------------------------------------
; *LINK->NUM3SX*
;-------------------------------------------------------------------------------
(defun link-hash (link) (sxhash (sort (copy-list link) #'<)))
(defun link-equal (link1 link2) (set-equal? link1 link2))
(sb-ext:define-hash-table-test link-equal link-hash)
(defparameter *LINK->NUM3SX* (make-hash-table :test 'link-equal)
  "the key is a link (p1 p2) and the value is the number of 3-simplices that 
(p1 p2) belongs to")

;-------------------------------------------------------------------------------
; use in breadth first search
;-------------------------------------------------------------------------------
(defparameter *VERTEX->COLOR* (make-hash-table)) ; 0=WHITE, 1=GRAY, 2=BLACK
(defparameter *VERTEX->DIST* (make-hash-table)) ; -1 MEANS INFINITY
(defparameter *VERTEX->PRED* (make-hash-table))
(defparameter *DIST->VERTICES* (make-hash-table) "key is distance in terms of
link length, value is a list of vertices at that distance, from the selected 
source vertex")

;-------------------------------------------------------------------------------
; used in auto correlation calculations
;-------------------------------------------------------------------------------
(defparameter *VERTEX->MULTI-PHIS* (make-hash-table)
  "this hashtable is used during the autocorrelation computation process; 
the key is a vertex, and the value is a list of phi values, one for each
sweep from lambada = 1 to Lambda. Hence the length of each list is Lambda")
(defparameter *TABLE3* (make-hash-table))
(defparameter *TABLE4* (make-hash-table))
(defparameter *TABLE5* (make-hash-table))
(defparameter *TABLE6* (make-hash-table))

;-------------------------------------------------------------------------------
; the maximum distance in this particular spacetime, from the currently 
; selected source vertex; basically the largest key in *DIST->VERTICES*
;-------------------------------------------------------------------------------
(defparameter *DMAX* 0)

; what is this?
(defparameter *LAMBDA* 0)

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
	   *VERTEX->NUM3SX*)

  ; update the bfs tables for the source vertex
  (setf (gethash s *VERTEX->COLOR*) 1) ; color the source vertex gray
  (setf (gethash s *VERTEX->DIST*) 0)
  (setf (gethash s *VERTEX->PRED*) nil)

  (let ((Q (make-instance 'queue))) ; a fifo queue to store the gray vertices
    (enqueue Q s)
    (while (not (queue-empty-p Q))
      (let ((u (dequeue Q)))
	(dolist (v (gethash u *VERTEX->VERTEX*))
	  (when (= 0 (gethash v *VERTEX->COLOR*)) ; found a white vertex
	    (setf (gethash v *VERTEX->COLOR*) 1) ; color this vertex gray
	    (setf (gethash v *VERTEX->DIST*) (+ 1 (gethash u *VERTEX->DIST*)))
	    (setf (gethash v *VERTEX->PRED*) u)
	    (enqueue Q v)))
  	(setf (gethash u *VERTEX->COLOR*) 2)))) ; color vertex u black
  (maphash #'(lambda (k v)
	       (pushnew k (gethash v *DIST->VERTICES*)))
	   *VERTEX->DIST*))

(defun compute-dmax ()
  "computes the maximum distance between vertices in a spacetime"
  (setf *DMAX* 0)
  (maphash #'(lambda (k v)
	       (declare (ignore v))
	       (when (> k *DMAX*)
		 (setf *DMAX* k)))
	   *DIST->VERTICES*))

(defun compute-phi-phi (x ys)
  "ys=(y1 ... yN); computes (phi(x)*phi(y1) + ... + phi(x)*phi(yN))/N. Here
the x and ys are the vertices, not the scalar field values at those vertices"
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

(defun compute-autocorrel (infilename big-lambda num-thermal)
  "computes the auto correlation time f(s) for the spacetime that infile 
contains. big-lambda is the size of the scalar field configuration collection.
This function will write an output file that has f(s) for s=0 to smax for this
particular spacetime"
  (format t "starting autocorrelation for spacetime ~a LAMBDA= ~A NUMTHERMAL = ~A at ~A~%"
	  infilename big-lambda num-thermal (cdt-now-str))
  ;; clear all existing tables
  (clrhash *ID->3SIMPLEX*)
  (clrhash *TL2SIMPLEX->ID*)
  (clrhash *SL2SIMPLEX->ID*)
  (clrhash *TL1SIMPLEX->ID*)
  (clrhash *SL1SIMPLEX->ID*)
  (clrhash *VERTEX->NUM3SX*)
  (clrhash *VERTEX->VERTEX*)
  (clrhash *VERTEX->PHI*)
  (clrhash *VERTEX->MULTI-PHIS*)
  (clrhash *LINK->NUM3SX*)
  (clrhash *VERTEX->COLOR*)
  (clrhash *VERTEX->DIST*)
  (clrhash *VERTEX->PRED*)
  (clrhash *DIST->VERTICES*)
  
  (setf *LAMBDA* big-lambda)

  ;; load the spacetime data
  (with-open-file (infile infilename :direction :input)
    (load-spacetime-from-file infile))

  ;; compute the value of alpha for this spacetime
  (compute-alpha-and-3sxvols)

  ;; populate scalar field tables
  (populate-scalar-field-tables)

  ;; thermalize
  (dotimes (ns num-thermal)
    (sweep))

  (format t "finished thermalizing sweeps ~A for spacetime ~a at ~A~%"
	  num-thermal infilename (cdt-now-str))

  ;; we are done thermalizing; push the current state of the scalar field
  ;; configuration to *VERTEX->MULTI-PHIS*
  (maphash #'(lambda (k v) (push v (gethash k *VERTEX->MULTI-PHIS*)))
	   *VERTEX->PHI*)
  
  ;; sweep and save for l = 1 to Lambda
  (for (l 1 *LAMBDA*)
       (sweep)

       (maphash #'(lambda (k v) (push v (gethash k *VERTEX->MULTI-PHIS*)))
		*VERTEX->PHI*))

  (format t "finished ~A configuration generation sweeps for spacetime ~A at ~A~%"
	  *LAMBDA* infilename (cdt-now-str))
  
  ;; TABLE1 IS *VERTEX->MULTI-PHIS*

  ;; now we are ready to compute the f(tau,s) for this spacetime
  (let ((s0 (select-random-site)))
    (breadth-first-search s0) ;; every vertex has been assigned a dist w.r.t s0
    (compute-dmax)
    ;; TABLE2 IS *DIST->VERTICES*

    (dotimes (d *DMAX*)
      (dotimes (l *LAMBDA*)
	(compute-phi-phi-v2 s0 l d)))
    (accumulate-each-row-of-table3)
    (compute_f_ac)
    (compute_2autocorreltime)
    (let* ((dotpos (position #\. infilename :from-end t))
	   (prefix (subseq infilename 0 dotpos))
	   (outfilename (concatenate 'string 
				     prefix 
				     "-AC-" 
				     (format nil "LAMBDA-~A-NUMTHERMAL-~A" 
					     *LAMBDA* num-thermal)
				     ".sfac2p1")))
      (with-open-file (outfile outfilename :direction :output)
	(maphash #'(lambda (k v) (format outfile "~A ~A~%" k v))
		 *TABLE6*)))))
    

(defun compute-phi-phi-v2 (v l d)
  (let* ((phi_x (nth l (gethash v *VERTEX->MULTI-PHIS*)))
	 (ys (gethash d *DIST->VERTICES*))
	 (phi_ys (mapcar #'(lambda (y)
			     (nth l (gethash y *VERTEX->MULTI-PHIS*))) 
			 ys))
	 (result (reduce #'+ (mapcar #'(lambda (phi_y) 
					 (* phi_x phi_y)) 
				     phi_ys))))
    (push (/ result (length phi_ys)) (gethash d *TABLE3*))))

(defun accumulate-each-row-of-table3 ()
  (dotimes (d *DMAX*)
    (setf (gethash d *TABLE4*) 
	  (/ (reduce #'+ (gethash d *TABLE3*)) *LAMBDA*))))

(defun compute_f_ac ()
  (dotimes (tau (1+ *LAMBDA*))
    (dotimes (d *DMAX*)
      (let ((sum 0.0))
	(dotimes (l *LAMBDA*)
	  (incf sum (* (- (nth l (gethash d *TABLE3*)) (gethash d *TABLE4*))
		       (- (nth (mod (+ l tau) *LAMBDA*) (gethash d *TABLE3*)) 
			  (gethash d *TABLE4*))))
	  (push (/ sum *LAMBDA*) (gethash d *TABLE5*)))))))

(defun compute_2autocorreltime ()
  "computes 2tauAC"
  (dotimes (d *DMAX*)
    (setf (gethash d *TABLE6*) (/ (reduce #'+ (gethash d *TABLE5*))
				  (nth *LAMBDA* (gethash d *TABLE5*))))))

(defun save-scalar-field-lattice-to-file (outfile)
  "
saves the scalar field lattice information to outfile. The following hash 
tables are saved:

*VERTEX->VERTEX*
*VERTEX->NUM3SX*
*LINK->NUM3SX*

This information is stored in the following format:

V0 N3 (V1 V2 V3) (N3_L01 N3_L02 N3_L02)

where v0 is the vertex, N3 is the number of 3-simplices that this vertex
belongs to, v1, v2, v3 are the vertices to which this link is connected,
N3_L01 is the number of 3-simplices attached to link (V0 V1), etc."

(maphash #'(lambda (k Vs)
	     (let* ((N3_V (gethash k *VERTEX->NUM3SX*))
		    (links (mapcar #'(lambda (k2) (list k k2)) Vs))
		    (N_Ls (mapcar #'(lambda (l) (gethash l *LINK->NUM3SX*))
				  links)))
	       (format outfile "~A ~A ~A ~A" k N3_V Vs N_Ls)))
	 *VERTEX->VERTEX*))

(defun parse-scalar-field-lattice-data-line (line)
  (with-input-from-string (s line)
    (let ((data (loop
		   :for num := (read s nil nil)
		   :while num
		   :collect num)))
      (setf (gethash (nth 0 data) *VERTEX->NUM3SX*) (nth 1 data))
      (setf (gethash (nth 0 data) *VERTEX->VERTEX*) (nth 2 data))
      (mapcar #'(lambda (v n)
		  (setf (gethash (list (nth 0 data) v) *LINK->NUM3SX*) n))
	      (nth 2 data) (nth 3 data)))))

(defun load-scalar-field-lattice-from-file (infile)
  (loop for line = (read-line infile nil)
     while line do (parse-scalar-field-lattice-data-line line)))

(defun generate-data (infilename numthermal numsweeps save-every-n-sweeps)
  "numthermal is the number of sweeps for thermalization and numsweeps is the
total number of sweeps."
  ;; clear all existing tables
  (clrhash *ID->3SIMPLEX*)
  (clrhash *TL2SIMPLEX->ID*)
  (clrhash *SL2SIMPLEX->ID*)
  (clrhash *TL1SIMPLEX->ID*)
  (clrhash *SL1SIMPLEX->ID*)
  (clrhash *VERTEX->NUM3SX*)
  (clrhash *VERTEX->VERTEX*)
  (clrhash *VERTEX->PHI*)
  (clrhash *LINK->NUM3SX*)
  (clrhash *VERTEX->COLOR*)
  (clrhash *VERTEX->DIST*)
  (clrhash *VERTEX->PRED*)
  (clrhash *DIST->VERTICES*)

  ;; load the spacetime data
  (with-open-file (infile infilename :direction :input)
    (load-spacetime-from-file infile))

  ;; populate scalar field tables
  (populate-scalar-field-tables)

  ;; thermalize
  (dotimes (ns numthermal)
    (sweep)
    (when (= 0 (mod ns 100))
      (format t "finished thermalizing sweep ~A of ~A for ~A at ~A~%"
	      ns numthermal infilename (cdt-now-str))))
  

  (let* ((dotpos (position #\. infilename :from-end t))
	 (prefix (subseq infilename 0 dotpos)))
    ;; save the thermalized setup as our initial data point
    (let ((outfilename (concatenate 'string 
				    prefix 
				    "-sfs-new" 
				    (format nil "~A" numthermal)
				    ".sf2p1"))
	  (source-vertex (select-random-site)))
      (breadth-first-search source-vertex)
      (with-open-file (outfile outfilename :direction :output)
	(print-phi-phi-expectation source-vertex outfile)))
    ;; continue sweeping, saving every 10000 sweeps
    (for (ns (1+ numthermal) numsweeps)
	 (sweep)
	 (when (= 0 (mod ns save-every-n-sweeps))
	   (format t "finished data sweep ~A of ~A for ~A at ~A~%"
		   ns numsweeps infilename (cdt-now-str)))
	 (when (= 0 (mod ns save-every-n-sweeps))
	   (let ((outfilename (concatenate 'string 
					   prefix 
					   "-sfs-new" 
					   (format nil "~A" ns)
					   ".sf2p1"))
		 (source-vertex (select-random-site)))
	     (breadth-first-search source-vertex)
	     (with-open-file (outfile outfilename :direction :output)
	       (print-phi-phi-expectation source-vertex outfile)))))))

(defun select-random-site ()
  "selects a random site"
  (let ((marker (random (hash-table-count *VERTEX->NUM3SX*)))
	(counter 0))
    (maphash #'(lambda (k v)
	       (declare (ignore v))
	       (if (= counter marker)
		   (return-from select-random-site k)
		   (incf counter)))
	     *VERTEX->NUM3SX*)))

(defun delta-action (v psi)
  "computes the change in the scalar field action if we change the field value
at site v from phi to psi"
  (let* ((phi (gethash v *VERTEX->PHI*))
	 (psisqr (* psi psi))
	 (phisqr (* phi phi))
	 ; N3_v is the number of 3-simplices associated with site v
	 (N3_v (gethash v *VERTEX->NUM3SX*))
	 ; l_v is the links associated with site v; l_v is a list of lists
	 (l_v (mapcar #'(lambda (k) (list v k))
		      (gethash v *VERTEX->VERTEX*)))
	 ; N3_ls is the number of 3-simplices associated with each link in l_v
         (N3_ls (mapcar #'(lambda (link)
			    (gethash link *LINK->NUM3SX*)) l_v))
	 ;; phijs is a list that holds the value of the scalar field at
	 ;; sites that are connected to v.
	 (phijs (mapcar #'(lambda (k) (gethash k *VERTEX->PHI*))
			(gethash v *VERTEX->VERTEX*)))
	 ;; if (f a1 a2) is supplied to mapcar, along with two lists l1 ans l2
	 ;; a1 is drawn from l1 and a2 is drawn from l2. jkterms holds each
	 ;; term of the first summation
	 (jkterms (mapcar #'(lambda (v phij)
			      (* v (- psisqr phisqr (* 2 phij (+ psi phi)))))
			  N3_ls phijs))
	 ;; sum all elements of jkterms, multiply it by a root2/144
	 (term1 (* (/ (* a (sqrt 2)) 144.0) (reduce #'+ jkterms)))
	 ;; term2 is the second term in the 
	 (term2 (* (* m m) (* a a a) (sqrt 2) N3_v (+ phisqr psisqr))))
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
  (for (num 1 (hash-table-count *VERTEX->NUM3SX*))
       (incf NUM-ATTEMPTED)
       (let ((site (select-random-site))
	     (newphi (rndval)))
	 (when (accept-move? site newphi)
           (incf NUM-ACCEPTED)
           (setf (gethash site *VERTEX->PHI*) newphi))))
  (* 100.0 (/ NUM-ACCEPTED NUM-ATTEMPTED)))

;--------1---------2---------3---------4---------5---------6---------7---------8
(defun parse-3simplex (3sx)
  (let ((p0 (nth-point 3sx 0))
	(p1 (nth-point 3sx 1))
	(p2 (nth-point 3sx 2))
	(p3 (nth-point 3sx 3)))
    ;
    (incf (gethash p0 *VERTEX->NUM3SX* 0))
    (incf (gethash p1 *VERTEX->NUM3SX* 0))
    (incf (gethash p2 *VERTEX->NUM3SX* 0))
    (incf (gethash p3 *VERTEX->NUM3SX* 0))
    ;
    (setf (gethash p0 *VERTEX->PHI*) (rndval))
    (setf (gethash p1 *VERTEX->PHI*) (rndval))
    (setf (gethash p2 *VERTEX->PHI*) (rndval))
    (setf (gethash p3 *VERTEX->PHI*) (rndval))
    ;
    (pushnew p1 (gethash p0 *VERTEX->VERTEX*))
    (pushnew p2 (gethash p0 *VERTEX->VERTEX*))
    (pushnew p3 (gethash p0 *VERTEX->VERTEX*))
    (pushnew p0 (gethash p1 *VERTEX->VERTEX*))
    (pushnew p2 (gethash p1 *VERTEX->VERTEX*))
    (pushnew p3 (gethash p1 *VERTEX->VERTEX*))
    (pushnew p0 (gethash p2 *VERTEX->VERTEX*))
    (pushnew p1 (gethash p2 *VERTEX->VERTEX*))
    (pushnew p3 (gethash p2 *VERTEX->VERTEX*))
    (pushnew p0 (gethash p3 *VERTEX->VERTEX*))
    (pushnew p1 (gethash p3 *VERTEX->VERTEX*))
    (pushnew p2 (gethash p3 *VERTEX->VERTEX*))
    ;
    (incf (gethash `(,p0 ,p1) *LINK->NUM3SX* 0))
    (incf (gethash `(,p0 ,p2) *LINK->NUM3SX* 0))
    (incf (gethash `(,p0 ,p3) *LINK->NUM3SX* 0))
    (incf (gethash `(,p1 ,p2) *LINK->NUM3SX* 0))
    (incf (gethash `(,p1 ,p3) *LINK->NUM3SX* 0))
    (incf (gethash `(,p2 ,p3) *LINK->NUM3SX* 0))))

(defun populate-scalar-field-tables ()
  "assumes spacetime has been loaded"
  (maphash #'(lambda (k v)
	       (declare (ignore k))
	       (parse-3simplex v))
	   *ID->3SIMPLEX*))

(defun print-phi-phi-expectation-console (s)
  "s is the source vertex"
  (maphash #'(lambda (k v)
	       (format t "~A,~A~%"
		       k (compute-phi-phi s v)))
	   *DIST->VERTICES*))

(defun generate-data-console (infilename numthermal numsweeps)
  "numthermal is the number of sweeps for thermalization and numsweeps is the
total number of sweeps."
  ;; clear all existing tables
  (clrhash *ID->3SIMPLEX*)
  (clrhash *TL2SIMPLEX->ID*)
  (clrhash *SL2SIMPLEX->ID*)
  (clrhash *TL1SIMPLEX->ID*)
  (clrhash *SL1SIMPLEX->ID*)
  (clrhash *VERTEX->NUM3SX*)
  (clrhash *VERTEX->VERTEX*)
  (clrhash *VERTEX->PHI*)
  (clrhash *LINK->NUM3SX*)
  (clrhash *VERTEX->COLOR*)
  (clrhash *VERTEX->DIST*)
  (clrhash *VERTEX->PRED*)
  (clrhash *DIST->VERTICES*)

  ;; load the spacetime data
  (with-open-file (infile infilename :direction :input)
    (load-spacetime-from-file infile))

  ;; populate scalar field tables
  (populate-scalar-field-tables)

  ;; thermalize
  (dotimes (ns numthermal)
    (sweep)
    (when (= 0 (mod ns 10))
      (format t "finished thermalizing sweep ~A of ~A for ~A at ~A~%"
	      ns numthermal infilename (cdt-now-str))))
  
  (let ((source-vertex (select-random-site)))
    (breadth-first-search source-vertex)
    (print-phi-phi-expectation-console source-vertex))

    ;; continue sweeping, saving every 10000 sweeps
    (for (ns (1+ numthermal) numsweeps)
	 (sweep)
	 (when (= 0 (mod ns 10))
	   (format t "finished data sweep ~A of ~A for ~A at ~A~%"
		   ns numsweeps infilename (cdt-now-str))
	   (let ((source-vertex (select-random-site)))
	     (breadth-first-search source-vertex)
	     (print-phi-phi-expectation-console source-vertex)))))