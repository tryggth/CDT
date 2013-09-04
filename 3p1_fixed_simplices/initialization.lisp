(declaim (optimize (speed 3)
		   (compilation-speed 0)
		   (debug 0)
		   (safety 0)))
;; cdt-3plus1-initialization-pbc.lisp --- R x S3 initialization with periodic 
;; boundary conditions

(defun initialize-S3-triangulation (num-time-intervals)
  (setf NUM-T num-time-intervals)
  (for (n 0 (1- (/ num-time-intervals 2)))
       (dolist (fivepts S3-1/2-41)
	 (if (= (length fivepts) 5)
	   (make-4simplex-v3 4 (* 2 n) (1+ (* 2 n))
	  		   (+ (* 2 n N0-PER-SLICE) (first fivepts))
			   (+ (* 2 n N0-PER-SLICE) (second fivepts))
			   (+ (* 2 n N0-PER-SLICE) (third fivepts)) 
			   (+ (* 2 n N0-PER-SLICE) (fourth fivepts))
			   (+ (* 2 n N0-PER-SLICE) (fifth fivepts)))
	   (make-4simplex-v3f 4 (* 2 n) (1+ (* 2 n))
	  		   (+ (* 2 n N0-PER-SLICE) (first fivepts))
			   (+ (* 2 n N0-PER-SLICE) (second fivepts))
			   (+ (* 2 n N0-PER-SLICE) (third fivepts)) 
			   (+ (* 2 n N0-PER-SLICE) (fourth fivepts))
			   (+ (* 2 n N0-PER-SLICE) (fifth fivepts)))))
       (dolist (fivepts S3-1/2-32)
	 (if (= (length fivepts) 5)
	   (make-4simplex-v3 3 (* 2 n) (1+ (* 2 n))
			   (+ (* 2 n N0-PER-SLICE) (first fivepts))
			   (+ (* 2 n N0-PER-SLICE) (second fivepts))
			   (+ (* 2 n N0-PER-SLICE) (third fivepts))
			   (+ (* 2 n N0-PER-SLICE) (fourth fivepts))
			   (+ (* 2 n N0-PER-SLICE) (fifth fivepts)))
	   (make-4simplex-v3f 3 (* 2 n) (1+ (* 2 n))
			   (+ (* 2 n N0-PER-SLICE) (first fivepts))
			   (+ (* 2 n N0-PER-SLICE) (second fivepts))
			   (+ (* 2 n N0-PER-SLICE) (third fivepts))
			   (+ (* 2 n N0-PER-SLICE) (fourth fivepts))
			   (+ (* 2 n N0-PER-SLICE) (fifth fivepts)))))
       (dolist (fivepts S3-1/2-23)
	 (if (= (length fivepts) 5)
	   (make-4simplex-v3 2 (* 2 n) (1+ (* 2 n))
	  		     (+ (* 2 n N0-PER-SLICE) (first fivepts))
			     (+ (* 2 n N0-PER-SLICE) (second fivepts))
			     (+ (* 2 n N0-PER-SLICE) (third fivepts)) 
			     (+ (* 2 n N0-PER-SLICE) (fourth fivepts))
			     (+ (* 2 n N0-PER-SLICE) (fifth fivepts)))
	   (make-4simplex-v3f 2 (* 2 n) (1+ (* 2 n))
	  		     (+ (* 2 n N0-PER-SLICE) (first fivepts))
			     (+ (* 2 n N0-PER-SLICE) (second fivepts))
			     (+ (* 2 n N0-PER-SLICE) (third fivepts)) 
			     (+ (* 2 n N0-PER-SLICE) (fourth fivepts))
			     (+ (* 2 n N0-PER-SLICE) (fifth fivepts)))))
       (dolist (fivepts S3-1/2-14)
	 (if (= (length fivepts) 5)
	   (make-4simplex-v3 1 (* 2 n) (1+ (* 2 n))
	  		     (+ (* 2 n N0-PER-SLICE) (first fivepts))
			     (+ (* 2 n N0-PER-SLICE) (second fivepts))
			     (+ (* 2 n N0-PER-SLICE) (third fivepts)) 
			     (+ (* 2 n N0-PER-SLICE) (fourth fivepts))
			     (+ (* 2 n N0-PER-SLICE) (fifth fivepts)))
	   (make-4simplex-v3f 1 (* 2 n) (1+ (* 2 n))
	  		     (+ (* 2 n N0-PER-SLICE) (first fivepts))
			     (+ (* 2 n N0-PER-SLICE) (second fivepts))
			     (+ (* 2 n N0-PER-SLICE) (third fivepts)) 
			     (+ (* 2 n N0-PER-SLICE) (fourth fivepts))
			     (+ (* 2 n N0-PER-SLICE) (fifth fivepts)))))
       (dolist (fivepts S3-1/2-41)
	 (if (= (length fivepts) 5)
	 (make-4simplex-v3 1 (1+ (* 2 n)) (2+ (* 2 n))
			   (+ (* 2 n N0-PER-SLICE) (fifth fivepts))
			   (+ (* 2 (+ n 1) N0-PER-SLICE) (first fivepts))
			   (+ (* 2 (+ n 1) N0-PER-SLICE) (second fivepts))
			   (+ (* 2 (+ n 1) N0-PER-SLICE) (third fivepts))
			   (+ (* 2 (+ n 1) N0-PER-SLICE) (fourth fivepts)))
	 (make-4simplex-v3f 1 (1+ (* 2 n)) (2+ (* 2 n))
			   (+ (* 2 n N0-PER-SLICE) (fifth fivepts))
			   (+ (* 2 (+ n 1) N0-PER-SLICE) (first fivepts))
			   (+ (* 2 (+ n 1) N0-PER-SLICE) (second fivepts))
			   (+ (* 2 (+ n 1) N0-PER-SLICE) (third fivepts))
			   (+ (* 2 (+ n 1) N0-PER-SLICE) (fourth fivepts)))))
       (dolist (fivepts S3-1/2-32)
	 (if (= (length fivepts) 5)
	 (make-4simplex-v3 2 (1+ (* 2 n)) (2+ (* 2 n))
			   (+ (* 2 n N0-PER-SLICE) (fourth fivepts))
			   (+ (* 2 n N0-PER-SLICE) (fifth fivepts))
			   (+ (* 2 (+ n 1) N0-PER-SLICE) (first fivepts))
			   (+ (* 2 (+ n 1) N0-PER-SLICE) (second fivepts))
			   (+ (* 2 (+ n 1) N0-PER-SLICE) (third fivepts)))
	 (make-4simplex-v3f 2 (1+ (* 2 n)) (2+ (* 2 n))
			   (+ (* 2 n N0-PER-SLICE) (fourth fivepts))
			   (+ (* 2 n N0-PER-SLICE) (fifth fivepts))
			   (+ (* 2 (+ n 1) N0-PER-SLICE) (first fivepts))
			   (+ (* 2 (+ n 1) N0-PER-SLICE) (second fivepts))
			   (+ (* 2 (+ n 1) N0-PER-SLICE) (third fivepts)))))
       (dolist (fivepts S3-1/2-23)
	 (if (= (length fivepts) 5)
	 (make-4simplex-v3 3 (1+ (* 2 n)) (2+ (* 2 n))
			   (+ (* 2 n N0-PER-SLICE) (third fivepts))
			   (+ (* 2 n N0-PER-SLICE) (fourth fivepts))
			   (+ (* 2 n N0-PER-SLICE) (fifth fivepts))
			   (+ (* 2 (+ n 1) N0-PER-SLICE) (first fivepts))
			   (+ (* 2 (+ n 1) N0-PER-SLICE) (second fivepts)))
	 (make-4simplex-v3f 3 (1+ (* 2 n)) (2+ (* 2 n))
			   (+ (* 2 n N0-PER-SLICE) (third fivepts))
			   (+ (* 2 n N0-PER-SLICE) (fourth fivepts))
			   (+ (* 2 n N0-PER-SLICE) (fifth fivepts))
			   (+ (* 2 (+ n 1) N0-PER-SLICE) (first fivepts))
			   (+ (* 2 (+ n 1) N0-PER-SLICE) (second fivepts)))))
       (dolist (fivepts S3-1/2-14)
	 (if (= (length fivepts) 5)
	 (make-4simplex-v3 4 (1+ (* 2 n)) (2+ (* 2 n))
			   (+ (* 2 n N0-PER-SLICE) (second fivepts))
			   (+ (* 2 n N0-PER-SLICE) (third fivepts))
			   (+ (* 2 n N0-PER-SLICE) (fourth fivepts))
			   (+ (* 2 n N0-PER-SLICE) (fifth fivepts))
			   (+ (* 2 (+ n 1) N0-PER-SLICE) (first fivepts)))
	 (make-4simplex-v3f 4 (1+ (* 2 n)) (2+ (* 2 n))
			   (+ (* 2 n N0-PER-SLICE) (second fivepts))
			   (+ (* 2 n N0-PER-SLICE) (third fivepts))
			   (+ (* 2 n N0-PER-SLICE) (fourth fivepts))
			   (+ (* 2 n N0-PER-SLICE) (fifth fivepts))
			   (+ (* 2 (+ n 1) N0-PER-SLICE) (first fivepts))))))
  
  (when (string= BCTYPE "PERIODIC")


    (for (ts 0 (- NUM-T 1))
	 (connect-simplices-in-sandwich ts (1+ ts))
	 (connect-simplices-in-adjacent-sandwiches ts (+ ts 1) (+ ts 2)))


    (set-last-used-pt (* NUM-T N0-PER-SLICE))


    (set-f-vector (* NUM-T N0-PER-SLICE) ;; N0
		  (* NUM-T N1-SL-PER-SLICE) ;; N1-SL
		  (* NUM-T N1-TL-PER-SLICE) ;; N1-TL
		  (* NUM-T N2-SL-PER-SLICE) ;; N2-SL
		  (* NUM-T N2-TL-PER-SLICE) ;; N2-TL
		  (* NUM-T N3-SL-PER-SLICE) ;; N3-SL
		  (* NUM-T (+ N3-TL-13-PER-SLICE N3-TL-31-PER-SLICE));; N3-TL-13
		  (* NUM-T N3-TL-22-PER-SLICE);; N3-TL-22
		  (* NUM-T (+ N4-TL-14-PER-SLICE N4-TL-41-PER-SLICE))
		  (* NUM-T (+ N4-TL-23-PER-SLICE N4-TL-32-PER-SLICE))))


  (when (string= BCTYPE "OPEN")
    (error "open boundary conditions not yet implemented")))


(defun initialize-T-slices-with-V-volume (&key num-time-slices target-volume 
					  boundary-conditions)
  (setf STOPOLOGY "S3")
  (setf BCTYPE (string-upcase boundary-conditions))
  (initialize-S3-triangulation num-time-slices)


  (format t "initial count = ~A~%" (count-simplices-of-all-types))


  ;; choice of (2,8) (2,4)v1 (2,4)v2 and (4,6) to increase the volume 


  (loop named tv
     do
     ;; increase the number of (1,4) and (4,1) simplices
       (dolist (id41 (get-simplices-of-type 4))
	 (let ((movedata nil))
	   (when (setf movedata (try-4->6 id41))
	     (3plus1move movedata)))
	 (if (>= (N4) target-volume)
	     (return-from tv)))
       
     ;; increase the number of (2,3) simplices
       (dolist (id23 (get-simplices-of-type 3))
	 (let ((movedata nil))
	   (when (setf movedata (try-2->4-v1 id23))
	     (3plus1move movedata)))
	 (if (>= (N4) target-volume)
	     (return-from tv)))
       
     ;; increase the number of (3,2) simplices
       (dolist (id32 (get-simplices-of-type 2))
	 (let ((movedata nil))
	   (when (setf movedata (try-2->4-v1 id32))
	     (3plus1move movedata)))
	 (if (>= (N4) target-volume)
	     (return-from tv)))


     ;; increase the number if (1,4) and (4,1) simplices
       (dolist (id14 (get-simplices-of-type 1))
	 (let ((movedata nil))
	   (when (setf movedata (try-2->8 id14))
	     (3plus1move movedata)))
	 (if (>= (N4) target-volume)
	     (return-from tv))))
    
  (format t "final count = ~A~%" (count-simplices-of-all-types))
  (setf N-INIT (N4)))

