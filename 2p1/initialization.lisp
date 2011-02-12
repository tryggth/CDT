;...........................................................................................................
; cdt-2plus1-initialization.lisp
;...........................................................................................................
(defun initialize-S2-triangulation (num-time-intervals)
  (setf NUM-T num-time-intervals)
  (for (n 0 (1- (/ num-time-intervals 2)))
       (dolist (fourpts S2-1/2-31)
	 (make-3simplex-v3 3 (* 2 n) (+ (* 2 n) 1)                 ;-----o------------ t = 1
			   (+ (* n 44) (first fourpts))            ;    / \
			   (+ (* n 44) (second fourpts))           ;   /   \
			   (+ (* n 44) (third fourpts))            ;  /     \
			   (+ (* n 44) (fourth fourpts))))         ;-o-------o-------- t = 0
       (dolist (fourpts S2-1/2-22)
	 (make-3simplex-v3 2 (* 2 n) (+ (* 2 n) 1)                 ;-----o-----o------ t = 1
			   (+ (* n 44) (first fourpts))            ;    /     / 
			   (+ (* n 44) (second fourpts))           ;   /     /     
			   (+ (* n 44) (third fourpts))            ;  /     /  
			   (+ (* n 44) (fourth fourpts))))         ;-o-----o---------- t = 0
       (dolist (fourpts S2-1/2-13)
	 (make-3simplex-v3 1 (* 2 n) (+ (* 2 n) 1)                 ;-o-------o-------- t = 1
			   (+ (* n 44) (first fourpts))            ;  \     /
			   (+ (* n 44) (second fourpts))           ;   \   /
			   (+ (* n 44) (third fourpts))            ;    \ /
			   (+ (* n 44) (fourth fourpts))))         ;-----o------------ t = 0
       (dolist (fourpts S2-1/2-31)
	 (make-3simplex-v3 1 (+ 1 (* 2 n)) (+ 2 (* 2 n))           ;-o-------o-------- t = 2
			   (+ (* n 44) (fourth fourpts))           ;  \     /
			   (+ (* (+ n 1) 44) (first fourpts))      ;   \   /
			   (+ (* (+ n 1) 44) (second fourpts))     ;    \ /
			   (+ (* (+ n 1) 44) (third fourpts))))    ;-----o------------ t = 1
       (dolist (fourpts S2-1/2-22)
	 (make-3simplex-v3 2 (+ 1 (* 2 n)) (+ 2 (* 2 n))           ;-----o-----o------ t = 2
			   (+ (* n 44) (third fourpts))            ;      \     \  
			   (+ (* n 44) (fourth fourpts))           ;       \     \
			   (+ (* (+ n 1) 44) (first fourpts))      ;        \     \
			   (+ (* (+ n 1) 44) (second fourpts))))   ;---------o-----o-- t = 1
       (dolist (fourpts S2-1/2-13)
	 (make-3simplex-v3 3 (+ 1 (* 2 n)) (+ 2 (* 2 n))           ;-----o------------ t = 2
			   (+ (* n 44) (second fourpts))           ;    / \
			   (+ (* n 44) (third fourpts))            ;   /   \
			   (+ (* n 44) (fourth fourpts))           ;  /     \
			   (+ (* (+ n 1) 44) (first fourpts)))))   ;-o-------o-------- t = 1
  
  
  ;; for periodic boundary conditions, the points in the NUM-T slice need to be identified with the points
  ;; in the 0 slice; this can be done by subtracting (* 22 NUM-T) from the the hi-points of all the
  ;; simplices in the ts=(- NUM-T 1/2) sandwich
  (when (string= BCTYPE "PERIODIC")

    (for (ts 0 (- NUM-T 1))
	 (connect-simplices-in-sandwich ts (1+ ts) )
	 (connect-simplices-in-adjacent-sandwiches ts (+ ts 1) (+ ts 2)))
  
    (set-last-used-pt (* NUM-T 22))
    
    ;; need to figure out N1-TL and N2-TL
    (set-f-vector (* NUM-T 22)                                                  ; N0
		  (* NUM-T 6)                                                  ; N1-SL
		  0                                                            ; N1-TL
		  (* NUM-T 4)                                                  ; N2-SL
		  0                                                            ; N2-TL
		  (+ (count-simplices-of-type 1) (count-simplices-of-type 3))  ; N3-TL-31
		  (count-simplices-of-type 2)))                                ; N3-TL-22

  (when (string= BCTYPE "OPEN")
    (error "open boundary conditions not yet implemented")))

(defun initialize-T2-triangulation (num-time-intervals)
  (setf NUM-T num-time-intervals)
  (error "IxT2 and S1xT2 spacetimes not yet implemented"))

(defun initialize-T-slices-with-V-volume (&key 
					  num-time-slices 
					  target-volume
					  spatial-topology
					  boundary-conditions)
  (setf BCTYPE (string-upcase boundary-conditions))
  (setf STOPOLOGY (string-upcase spatial-topology))

  (when (string= STOPOLOGY "S2")
    (initialize-S2-triangulation num-time-slices))
  (when (string= STOPOLOGY "T2")
    (initialize-T2-triangulation num-time-slices))

  (when (and (string= STOPOLOGY "S2") (string= BCTYPE "PERIODIC"))
    (setf (symbol-function 'action) #'action-S1xS2))

  (format t "initial count = ~A~%" (count-simplices-of-all-types))
  (loop named repeat
     do 
       (dolist (id23 (get-simplices-of-type 2))
	 (when (get-3simplex id23)
	   (let ((movedata nil))
	     (when (setf movedata (try-2->3 id23))
	       (2plus1move movedata))))
	 (if (>= (N3) target-volume)
	     (return-from repeat)))
       
       (dolist (id26 (get-simplices-of-type 1))
	 (when (get-3simplex id26)
	   (let ((movedata nil))
	     (when (setf movedata (try-2->6 id26))
	       (2plus1move movedata))))
	 (if (>= (N3) target-volume)
	     (return-from repeat)))
       
       (dolist (id23 (get-simplices-of-type 2))
	 (when (get-3simplex id23)
	   (let ((movedata nil))
	     (when (setf movedata (try-2->3 id23))
	       (2plus1move movedata))))
	 (if (>= (N3) target-volume)
	     (return-from repeat))))

  (format t "final count = ~A~%" (count-simplices-of-all-types))

  (setf N-INIT (N3)))
