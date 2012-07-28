(load "cdt2p1.lisp")

;;;; Constants
;;;;------------------------------------------------------------------------
(defvar *target-volume* 30850)
;;;;------------------------------------------------------------------------

;;;; Functions
;;;;------------------------------------------------------------------------
(defun extract-points-from-simplices (simplex-list)
  (mapcar #'(lambda (x) (first (rest x))) simplex-list))

(defun triangulate-between-s2-slices-special (original-sheet new-sheet 
					      t0 t1)

  (let* ((adjusted-new-sheet new-sheet)
	 (original-pseudo-faces-and-edges 
	  (get-s2-pseudo-faces-and-edges original-sheet))
	 (new-pseudo-faces-and-edges      
	  (get-s2-pseudo-faces-and-edges adjusted-new-sheet))

	 ;;make convenient bindings for the points, faces and edges
	 (original-points (first  original-pseudo-faces-and-edges))
	 (original-faces  (second original-pseudo-faces-and-edges))
	 (original-edges  (third  original-pseudo-faces-and-edges))
	 (new-points      (first  new-pseudo-faces-and-edges))
	 (new-faces       (second new-pseudo-faces-and-edges))
	 (new-edges       (third  new-pseudo-faces-and-edges))

	 ;;verticies
	 (point1 (first  original-points))
	 (point2 (second original-points))
	 (point3 (third  original-points))
	 (point4 (fourth original-points))
	 (point5 (first  new-points))
	 (point6 (second new-points))
	 (point7 (third  new-points))
	 (point8 (fourth new-points))

	 ;;faces
	 (face1 (first  original-faces))
	 (face2 (second original-faces))
	 (face3 (third  original-faces))
	 (face4 (fourth original-faces))
	 (face5 (first  new-faces))
	 (face6 (second new-faces))
	 (face7 (third  new-faces))
	 (face8 (fourth new-faces))
	 
	 ;;name edges in terms of which points they connect
	 (edge12 (first  original-edges))
	 (edge13 (second original-edges))
	 (edge14 (third  original-edges))
	 (edge23 (fourth original-edges))
	 (edge24 (fifth  original-edges))
	 (edge34 (sixth  original-edges))
	 (edge56 (first  new-edges))
	 (edge57 (second new-edges))
	 (edge58 (third  new-edges))
	 (edge67 (fourth new-edges))
	 (edge68 (fifth  new-edges))
	 (edge78 (sixth  new-edges)))

    ;;use the point/face/edge bindings and the triangulation between
    ;;tetrahedra scheme to come up with the complices that correctly
    ;;triangulate between the two sheets
    (concatenate 'list

		 ;;(3,1) complices:
		 (triangulate-31-complex face4 point5 t0 t1)
		 (triangulate-31-complex face1 point6 t0 t1)
		 (triangulate-31-complex face2 point7 t0 t1)
		 (triangulate-31-complex face3 point8 t0 t1)

		 ;;(2,2) complices:
		 (triangulate-22-complex edge12 edge58 t0 t1)
		 (triangulate-22-complex edge23 edge56 t0 t1)
		 (triangulate-22-complex edge13 edge57 t0 t1)
		 (triangulate-22-complex edge34 edge67 t0 t1)
		 (triangulate-22-complex edge24 edge68 t0 t1)
		 (triangulate-22-complex edge14 edge78 t0 t1)

		 ;;(1,3) complices:
		 (triangulate-13-complex point1 face6 t0 t1)
		 (triangulate-13-complex point2 face7 t0 t1)
		 (triangulate-13-complex point3 face8 t0 t1)
		 (triangulate-13-complex point4 face5 t0 t1))))



(defun initialize-spacetime-and-identify-ends (num-time-slices 
					       &optional
						 initial-spatial-geometry)
					
  (setf NUM-T num-time-slices)
   
  ;;set up the numbers of initial simplices according to the S2
  ;;initialization of the geometry
  (defparameter N0-PER-SLICE 4)
  (defparameter N1-SnL-PER-SLICE 6)
  (defparameter N1-TL-PER-SLICE 12)
  (defparameter N2-SL-PER-SLICE 4)
  (defparameter N2-TL-PER-SLICE 24)
  (defparameter N3-TL-13-PER-SLICE 4)
  (defparameter N3-TL-22-PER-SLICE 6)
  ; here (3,1) means just (3,1), not (3,1)+(1,3)
  (defparameter N3-TL-31-PER-SLICE 4) 
  (defparameter S2-1/2-31 '((1 2 3 5) (2 3 4 6) (3 4 1 7) (4 1 2 8)))
  (defparameter S2-1/2-22 '((1 2 5 8) (2 3 5 6) (3 1 5 7) (3 4 6 7) 
			    (4 2 6 8) (4 1 7 8)))
  (defparameter S2-1/2-13 '((1 5 7 8) (2 5 6 8) (3 5 6 7) (4 6 7 8)))

  (for (n 0 (1- (/ NUM-T 2)))
       (dolist (fourpts S2-1/2-31)
	 (make-3simplex-v3 3 (* 2 n) (1+ (* 2 n))                      
			   (+ (* 2 n N0-PER-SLICE) (first fourpts))    
			   (+ (* 2 n N0-PER-SLICE) (second fourpts))   
			   (+ (* 2 n N0-PER-SLICE) (third fourpts))    
			   (+ (* 2 n N0-PER-SLICE) (fourth fourpts)))) 
       (dolist (fourpts S2-1/2-22)
	 (make-3simplex-v3 2 (* 2 n) (1+ (* 2 n))                       
			   (+ (* 2 n N0-PER-SLICE) (first fourpts))    
			   (+ (* 2 n N0-PER-SLICE) (second fourpts))   
			   (+ (* 2 n N0-PER-SLICE) (third fourpts))    
			   (+ (* 2 n N0-PER-SLICE) (fourth fourpts)))) 
       (dolist (fourpts S2-1/2-13)
	 (make-3simplex-v3 1 (* 2 n) (1+ (* 2 n))                       
			   (+ (* 2 n N0-PER-SLICE) (first fourpts))    
			   (+ (* 2 n N0-PER-SLICE) (second fourpts))    
			   (+ (* 2 n N0-PER-SLICE) (third fourpts))     
			   (+ (* 2 n N0-PER-SLICE) (fourth fourpts))))  
       (dolist (fourpts S2-1/2-31)
	 (make-3simplex-v3 1 (1+ (* 2 n)) (2+ (* 2 n))                      
			   (+ (* 2 n N0-PER-SLICE) (fourth fourpts))        
			   (+ (* 2 (+ n 1) N0-PER-SLICE) (first fourpts))   
			   (+ (* 2 (+ n 1) N0-PER-SLICE) (second fourpts))  
			   (+ (* 2 (+ n 1) N0-PER-SLICE) (third fourpts)))) 
       (dolist (fourpts S2-1/2-22)
	 (make-3simplex-v3 2 (1+ (* 2 n)) (2+ (* 2 n))                      
			   (+ (* 2 n N0-PER-SLICE) (third fourpts))         
			   (+ (* 2 n N0-PER-SLICE) (fourth fourpts))        
			   (+ (* 2 (+ n 1) N0-PER-SLICE) (first fourpts))   
			   (+ (* 2 (+ n 1) N0-PER-SLICE) (second fourpts))))
       (dolist (fourpts S2-1/2-13)
	 (make-3simplex-v3 3 (1+ (* 2 n)) (2+ (* 2 n))                      
			   (+ (* 2 n N0-PER-SLICE) (second fourpts))        
			   (+ (* 2 n N0-PER-SLICE) (third fourpts))         
			   (+ (* 2 n N0-PER-SLICE) (fourth fourpts))        
			   (+ (* 2 (+ n 1) N0-PER-SLICE) (first fourpts)))))

  ;; Replaces the correct time-slices
  (when initial-spatial-geometry
    (let ((existing-triangles ()))
      ;; Remove any triangles that exist in the initial slice and
      ;; make sure that the triangles they're connected to in the
      ;; next time slice are recorded.
      (dolist (3sxid (get-simplices-in-sandwich 0 1))
	(let* ((3sx  (get-3simplex 3sxid))
	       (pts (3sx-points 3sx)))
	  (when (= 1 (3sx-type 3sx))
	    (push (list (second pts) (third pts) (fourth pts)) 
		  existing-triangles)))
	(remove-3simplex 3sxid))
      ;; Remove any subsimplices on the zeroth time slice
      (remove-sl2simplices (get-spacelike-triangles-at-time 0))
      (remove-sl1simplices (get-spacelike-links-at-time 0))
      ;; Remove any subsimplices in the first time sandwich
      (remove-tl2simplices (get-timelike-triangles-in-sandwich 0 1))
      (remove-tl1simplices (get-timelike-links-in-sandwich 0 1))
      ;; Make simplexes in the first slice based on the initial
      ;; spatial geometry.
      (map 'list 
	   #'(lambda (x) (make-3simplex-v3 (first x) (second x) (third x) 
					   (fourth x) (fifth x) (sixth x) 
					   (seventh x)))
	   (triangulate-between-slices existing-triangles 
				       (load-triangles-from-file 
					initial-spatial-geometry)
				       1 0 *LAST-USED-POINT*))))

  (print "First algorithm complete!")

  (let ((existing-triangles ()))
      ;; Remove simplices in the final sandwich.
      (dolist (3sxid (get-simplices-in-sandwich (1- NUM-T) NUM-T))
	(let* ((3sx  (get-3simplex 3sxid))
	       (pts (3sx-points 3sx)))
	  (when (= 3 (3sx-type 3sx))
	    (push (list (first pts) (second pts) (third pts))
		  existing-triangles)))
	(remove-3simplex 3sxid))
      ;; Remove any subsimplices on the NUM-T time slice
      (remove-sl2simplices (get-spacelike-triangles-at-time NUM-T))
      (remove-sl1simplices (get-spacelike-links-at-time NUM-T))
      ;; Remove any subsimplices in the final time sandwich
      (remove-tl2simplices (get-timelike-triangles-in-sandwich
			    (1- NUM-T) NUM-T))
      (remove-tl1simplices (get-timelike-links-in-sandwich
			    (1- NUM-T) NUM-T))
      (setf BCTYPE "PERIODIC")
      ;; Map the initial slice and the final slice to each other.

      (map 'list 
	   #'(lambda (x) (make-3simplex-v3 (first x) (second x) 0
					   (fourth x) (fifth x) (sixth x)
					   (seventh x)))
	   (triangulate-between-s2-slices-special
	    existing-triangles 
	    (extract-points-from-simplices 
	     (get-spacelike-triangles-at-time 0))
	     (1- NUM-T) NUM-T))))

(defun initialize-f-vector ()
       (set-f-vector 
      ;; nSIMPLICES where n<3
      (count-points-in-spacetime)              ;; N0
      (count-spacelike-links-in-spacetime)     ;; N1-SL 
      (count-timelike-links-in-spacetime)      ;; N1-TL
      (count-spacelike-triangles-in-spacetime) ;; N2-SL
      (count-timelike-triangles-in-spacetime)  ;; N2-TL
      
      ;;N3 = N3-TL-13 + N3-TL-31
      (+ (count-simplices-of-type 3)
	 (count-simplices-of-type 1))
      
      ;;N3-TL-22
      (count-simplices-of-type 2)))

(defun grow-spacetime (target-volume)
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
  
  (format t "breakdown by location = ~a~%" (count-boundary-vs-bulk))
  
  (setf N-INIT (N3)))
;;;;------------------------------------------------------------------------



;;;; Script
;;;;------------------------------------------------------------------------
;; Sets sweeps
(setf NUM-SWEEPS 50000)

;; Sets topology
(setf STOPOLOGY "S2")

;; Make the spacetime have open boundary conditions for initialization
(setf BCTYPE "OPEN")

;; Initializes spacetime (sort of)
(initialize-spacetime-and-identify-ends 64 "tetra.txt")

;; Forces identification between initial and final slices
(for (ts 0  (1- NUM-T))
  (connect-simplices-in-sandwich ts (1+ ts))
  (connect-simplices-in-adjacent-sandwiches ts (+ ts 1) (+ ts 2)))

(initialize-f-vector)

;; Set's the B-Vector to zeroes everywhere, ensuring it doesn't affect
;; the dynamics.
(set-b-vector 0 0 0 0 0 0)

;; Resets boundary type to periodic
(setf BCTYPE "PERIODIC")

;; Tell me whats going on
(format t "initial count = ~A~%" (count-simplices-of-all-types))

;; Grow spacetime to the chosen volume
(grow-spacetime *targ et-volume*)

;; Set coupling constants
(set-k0-k3-alpha 1.0 0.7577 -1)

;; Gather data
(generate-spacetime-and-movie-data)
