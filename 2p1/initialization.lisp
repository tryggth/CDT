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

       (dolist (id26 (get-simplices-of-type 3))
	 (let ((movedata nil))
	   (when (setf movedata (try-2->6 id26))
	     (2plus1move movedata)))
	 (if (> (N3) target-volume)
	     (return-from tv)))

       (dolist (id23 (get-simplices-of-type 2))
	 (let ((movedata nil))
	   (when (setf movedata (try-2->3 id23))
	     (2plus1move movedata)))
	 (if (> (N3) target-volume)
	     (return-from tv))))

  (format t "final count = ~A~%" (count-simplices-of-all-types))

  (setf N-INIT (N3)))
