; cdt-3plus1-initialization-pbc.lisp --- R x S3 initialization with periodic boundary conditions

(defun initialize-S3-triangulation-pbc (num-time-intervals)

  (setf NUM-T num-time-intervals)
  (setf BCTYPE 'periodic)

;  (do ((n 0 (1+ n))) ((>= n (/ num-time-intervals 2)))

  (for (n 0 (1- (/ num-time-intervals 2)))
       ;; t = 1/2, 5/2, 9/2, ...
       (dolist (fivepts S3-1/2-41)
	 (make-4simplex 4 (+ 1/2 (* 2 n))
			(+ (* n 10) (first fivepts))
			(+ (* n 10) (second fivepts))
			(+ (* n 10) (third fivepts))
			(+ (* n 10) (fourth fivepts))
			(+ (* n 10) (fifth fivepts))))
       (dolist (fivepts S3-1/2-32)
	 (make-4simplex 3 (+ 1/2 (* 2 n))
			(+ (* n 10) (first fivepts))
			(+ (* n 10) (second fivepts))
			(+ (* n 10) (third fivepts))
			(+ (* n 10) (fourth fivepts))
			(+ (* n 10) (fifth fivepts))))
       (dolist (fivepts S3-1/2-23)
	 (make-4simplex 2 (+ 1/2 (* 2 n))
			(+ (* n 10) (first fivepts))
			(+ (* n 10) (second fivepts))
			(+ (* n 10) (third fivepts))
			(+ (* n 10) (fourth fivepts))
			(+ (* n 10) (fifth fivepts))))
       (dolist (fivepts S3-1/2-14)
	 (make-4simplex 1 (+ 1/2 (* 2 n))
			(+ (* n 10) (first fivepts))
			(+ (* n 10) (second fivepts))
			(+ (* n 10) (third fivepts))
			(+ (* n 10) (fourth fivepts))
			(+ (* n 10) (fifth fivepts))))
       
       ;; t = 3/2, 7/2, 11/2, ...
       (dolist (fivepts S3-1/2-41)
	 (make-4simplex 1 (+ 3/2 (* 2 n))
			(+ (* n 10) (fifth fivepts))
			(+ (* (+ n 1) 10) (first fivepts))
			(+ (* (+ n 1) 10) (second fivepts))
			(+ (* (+ n 1) 10) (third fivepts))
			(+ (* (+ n 1) 10) (fourth fivepts))))
       
       (dolist (fivepts S3-1/2-32)
	 (make-4simplex 2 (+ 3/2 (* 2 n))
			(+ (* n 10) (fourth fivepts))
			(+ (* n 10) (fifth fivepts))
			(+ (* (+ n 1) 10) (first fivepts))
			(+ (* (+ n 1) 10) (second fivepts))
			(+ (* (+ n 1) 10) (third fivepts))))
       
       (dolist (fivepts S3-1/2-23)
	 (make-4simplex 3 (+ 3/2 (* 2 n))
			(+ (* n 10) (third fivepts))
			(+ (* n 10) (fourth fivepts))
			(+ (* n 10) (fifth fivepts))
			(+ (* (+ n 1) 10) (first fivepts))
			(+ (* (+ n 1) 10) (second fivepts))))
       
       (dolist (fivepts S3-1/2-14)
	 (make-4simplex 4 (+ 3/2 (* 2 n))
			(+ (* n 10) (second fivepts))
			(+ (* n 10) (third fivepts))
			(+ (* n 10) (fourth fivepts))
			(+ (* n 10) (fifth fivepts))
			(+ (* (+ n 1) 10) (first fivepts)))))
  
  
  ;; for periodic boundary conditions, the points in the NUM-T slice need to be identified with the points
  ;; in the 0 slice; this can be done by subtracting (* 5 NUM-T) from the the hi-points of all the
  ;; simplices in the ts=(- NUM-T 1/2) sandwich
  (maphash #'(lambda (k v)
	       (declare (ignore k))
	       (if (= (- NUM-T 1/2) (second v))
		   (cond ((= 1 (first v))
			  (decf (second (third v)) (* 5 NUM-T))
			  (decf (third (third v)) (* 5 NUM-T))
			  (decf (fourth (third v)) (* 5 NUM-T))
			  (decf (fifth (third v)) (* 5 NUM-T)))
			 ((= 2 (first v))
			  (decf (third (third v)) (* 5 NUM-T))
			  (decf (fourth (third v)) (* 5 NUM-T))
			  (decf (fifth (third v)) (* 5 NUM-T)))
			 ((= 3 (first v))
			  (decf (fourth (third v)) (* 5 NUM-T))
			  (decf (fifth (third v)) (* 5 NUM-T)))
			 ((= 4 (first v))
			  (decf (fifth (third v)) (* 5 NUM-T))))))
	   *4SIMPLEX-STORE*)

  ;; connect the 4simplices in all the sandwiches including the last one; T+1/2 mod T equals 1/2. this
  ;; pairs up the T-1/2 sandwich with the 1/2 sandwich, which is what we want for pbc
  (for (ts 1/2 (- NUM-T 1/2))
       (connect-4simplices-in-slice ts)
       (connect-4simplices-in-adjacent-slices ts (+ ts 1)))
  
  ;; if num-t = 4, we have t=0,1,2,3 and use up points upto 20
  ;; if num-t = 6, we have t=0,1,2,3,4,5 and use up points upto 30
  ;; if num-t = 8, .... use up points upto 40
  (setf NEXT-PT (* NUM-T 5))
  
  ;; now to connect the spatial tetrahedra to the (1,4) and (4,1) simplicies
  ;; when ts = NUM-T - 1/2, (mod (+ ts 1/2) NUM-T) equals 0, which is what we want for periodic b.c.
  (for (ts 1/2 (- NUM-T 1/2))
       (connect-3simplices-to-4simplices (get-3simplices-at-time (- ts 1/2))
					 (get-4simplices-of-type-at-time 4 ts))
       (connect-3simplices-to-4simplices (get-3simplices-at-time (mod (+ ts 1/2) NUM-T))
					 (get-4simplices-of-type-at-time 1 ts)))

  ;; need to set the f-vector correctly
)


(defun initialize-T-slices-with-V-volume-pbc (&key num-time-slices target-volume)
  (setf SPATIAL-TOPOLOGY "S3")
  (initialize-S3-triangulation-pbc num-time-slices)

  (format t "initial count = ~A~%" (count-4simplices-of-all-types))
  (finish-output)
  ;; choice of (2,8) (2,4)v1 (2,4)v2 and (4,6) to increase the volume 
  (while (<= (N4) target-volume)
    
    ;; increase the number of (1,4) and (4,1) simplices
    (dolist (id41 (get-4simplices-of-type 4))
      (when (get-4simplex id41)
	(46-move id41)))

    ;; increase the number of (2,3) simplices
    (dolist (id23 (get-4simplices-of-type 2))
      (when (get-4simplex id23)
	(24-move-v1 id23)))
    
    ;; increase the number of (2,3) and (3,2) simplices
    (dolist (id23 (get-4simplices-of-type 2))
      (when (get-4simplex id23)
    	(24-move-v2 id23)))
    
    ;; increase the number of (3,2) simplices
    (dolist (id32 (get-4simplices-of-type 3))
      (when (get-4simplex id32)
	(24-move-v1 id32)))

    ;; increase the number if (1,4) and (4,1) simplices
    (dolist (id14 (get-4simplices-of-type 1))
      (when (get-4simplex id14)
	(46-move id14)))
    
    (format t "intermediate count ~A~%" (count-4simplices-of-all-types))
    (finish-output))
    
  (format t "final count = ~A~%" (count-4simplices-of-all-types))

  (setf N-INIT (N4)))