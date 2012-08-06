;............................................................................
; cdt-2plus1-initialization.lisp
;............................................................................

;; Authors:
;; ------- Rajesh Kommu
;; ------- David Kamensky
;; ------- Jonah Miller (jonah.maxwell.miller@gmail.com)

;;function to find neighbors of a given triangle in a list of
;;triangles (triangles are stored as 3-tuples of vertices)
(defun get-nbors-of-triangle (key-triangle list-of-triangles)
  (let ((retval ()))
    (dolist (tri list-of-triangles)
      (when (= 2 (length (intersection tri key-triangle)))
	(setf retval (push tri retval))))
    retval))

;;find all triangles containing a particular point.  triangles are
;;represented in a list of 3-tuples
(defun triangles-around-point (point list-of-triangles)
  (let ((retval ()))
    (dolist (tri list-of-triangles)
      (when (intersection tri (list point))
	(setf retval (push tri retval))))
    retval))

;;function to retrieve "pseudo-faces" and "pseudo-edges" from a list
;; of triangle verticies, of the form ((id1 id2 id3) ...).

;;the return value is of the form 

;;  ((point1 point2 point3 point4) (face1 face2 face3 face4) (edge12
;;  edge13 edge14 edge23 edge24 edge34))

;;where faces are lists of 3-tuples of vertices and edges are lists of
;;pairs of vertices.  the numbers in the edge names refer to the
;;numbers of the end points.  the numbers of points are the same as
;;the numbers of the faces opposite the points.  point* are just the
;;point ids of the points that are the verticies of the "tetrahedra"
;;that we are pretending the s2 triangulation is.

(defun get-s2-pseudo-faces-and-edges (triangle-sheet)

  ;;TODO come up with a more even way of splitting the sheet into faces

  (let* (;;pick a random triangle for face1
	 (face1 (list (nth (random (length triangle-sheet)) triangle-sheet)))
	 
	 ;; choose a random neighbor of that triangle to be face2.  
	 ;; the ;vertex of face2 not shared with face1 becomes point1.  
	 ;; the vertex of face1 not shared with face2 becomes point2
	 (face2  (list (first (get-nbors-of-triangle 
			       (first face1) triangle-sheet))))
	 (point1 (first (set-difference (first face2) (first face1))))
	 (point2 (first (set-difference (first face1) (first face2))))
		 
	 ;;choose all other triangles meeting at a random endpoint of
	 ;;the border line segment for face3.  the endpoint selected
	 ;;becomes point4, the other endpoint becomes point3.  the
	 ;;edge along the rounded part of this "triangle fan" is the
	 ;;only non-trivial edge (edge12) to construct.
	 (border-segment (intersection (first face1) (first face2)))
	 (point3 (first  border-segment))
	 (point4 (second border-segment))
	 ;(assuming x and y have the same cardinality)
	 (face3 (set-difference
		 (triangles-around-point point4 triangle-sheet)
		 (concatenate 'list face1 face2)
		 :test #'(lambda (x y) (not (set-difference x y))))) 
	 (edge12 (map 'list 
		      #'(lambda (x) (set-difference x (list point4))) face3))

	 ;;come up with other edges
	 (edge13 (list (list point1 point3)))
	 (edge14 (list (list point1 point4)))
	 (edge23 (list (list point2 point3)))
	 (edge24 (list (list point2 point4)))
	 (edge34 (list (list point3 point4)))
	 
	 ;;all other triangles form face4
	 (face4 (set-difference triangle-sheet 
				(concatenate 'list face1 face2 face3) 
				:test #'(lambda (x y) 
					  (not (set-difference x y))))))

    ;;get all of the values into a list to return
    (list (list point1 point2 point3 point4)
	  (list face1  face2  face3  face4 )
	  (list edge12 edge13 edge14 edge23 edge24 edge34))))

;;an analogous function to get-s2-pseudo-faces-and-edges, but for
;;triangle-sheet of t2 topology
(defun get-t2-pseudo-faces-and-edges (triangle-sheet)
  (declare (ignore triangle-sheet))
  (error "t2 topology not implemented yet"))

;;functions/macros to triangulate different types of complices into
;;many simplices of the same type.  the simplices are returned as a
;;list of lists, each with the form (type tmlo tmhi id1 id2 id3 id4),
;;for easy use with make-simplex-v3.
; to prevent errors
(defun triangulate-31-complex (face point tmlo tmhi)
  (if (> tmlo tmhi) ;This causes errors but is harmless.
      (triangulate-13-complex point face tmhi tmlo) 
      (map 'list #'(lambda (x) (list 3 tmlo tmhi 
				     (first x) 
				     (second x) 
				     (third x) point)) face)))

(defun triangulate-13-complex (point face tmlo tmhi)
  (if (> tmlo tmhi)
      (triangulate-31-complex face point tmhi tmlo)
      (map 'list #'(lambda (x) (list 1 tmlo tmhi point 
				     (first x) 
				     (second x) 
				     (third x))) face)))

(defun triangulate-22-complex (edge1 edge2 tmlo tmhi)
  (if (> tmlo tmhi)
      (triangulate-22-complex edge2 edge1 tmhi tmlo)
      (let ((retval ()))
	(dolist (lo-edge edge1)
	  (dolist (hi-edge edge2)
	    (push (list 2 tmlo tmhi 
			(first lo-edge) (second lo-edge) 
			(first hi-edge) (second hi-edge))
		  retval)))
	retval)))

;; original-sheet and new-sheet are lists of 3-tuples of pt ids, to be
;; connected by tetrahedra.  the function will return a list of
;; tetrahedra connecting the two sheets of triangles, such that all
;; tetrahedra meet at time-like faces with one another.  both sheets
;; are assumed to have topology s2.

;; the list of returned tetrahedra is of the form 

;; ((type tmlo tmhi id1 id2 id3 id4) ... ), 

;; for use with make-simplex-v3.  the call to make-simplex will
;; automatically update the last-used point id as tethrahedra are
;; added to the hash table.

;; point ids in original-sheet will not be changed.  values in
;; new-sheet will be set to reflect connectivity with original-sheet.
;; t0 will become the new time value of all verticies in
;; original-sheet and t1 will be the new time value for verticies in
;; new-sheet.  last-used-pt-id is the last used id of points and will
;; be added to all id's in new-sheet to prevent spurious connections
;; to other geometry

(defun triangulate-between-s2-slices (original-sheet new-sheet 
				      t0 t1 last-used-pt-id)

  ;;the approach here is to break each s2-topology sheet of triangles
  ;;into four "pseudo-faces", connected analogously to those of a
  ;;single tetrahedron.  we can then follow the program for filling
  ;;the space between two tetrahedra, with the slight complication
  ;;that the 3-simplices of the inter-tetrahedral filling are replaced
  ;;with analogous complices, each of which may be decomposed into
  ;;many simplices of the same type.

  ;;add last-used-pt-id to all verticies of new-sheet
  (let* ((adjusted-new-sheet 
	  (map 
	   'list 
	   #'(lambda (x) (map 
			  'list 
			  #'(lambda (y) (+ y last-used-pt-id)) x)) new-sheet))
	 
	 ;;first, we must decide on the pseudo-faces and pseudo-edges.
	 ;;Store the pseudo-faces as lists of triangles, and
	 ;;pseudo-edges as lists of pairs of points.

	 ;;get the points, faces, ad edges
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

	 ;;come up with numbered names for easier correspondence with
	 ;;the inter-tetrahedral filling scheme (vertex numbers
	 ;;opposite face numbers)

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

;;anaologous function for triangle sheets of t2 topology
(defun triangulate-between-t2-slices (original-sheet new-sheet 
				      t0 t1 last-used-pt-id)
  (declare (ignore original-sheet)
	   (ignore new-sheet)
	   (ignore t0)
	   (ignore t1)
	   (ignore last-used-pt-id))
  (error "t2 topology not yet implemented"))

;;macro to automatically choose the right triangulation between slices
(defmacro triangulate-between-slices (original-sheet new-sheet 
				      t0 t1 last-used-pt-id)
  `(cond 
     ((string= STOPOLOGY "S2")
      (triangulate-between-s2-slices ,original-sheet ,new-sheet 
				     ,t0 ,t1 ,last-used-pt-id))
     ((string= STOPOLOGY "T2")
      (triangulate-between-t2-slices ,original-sheet ,new-sheet 
				     ,t0 ,t1 ,last-used-pt-id))
     (t (error "unrecognized spatial topology"))))

;;test function to generate an s2 triangulation of arbitrary size.
;;this particular method focuses a lot of curvature at a few points,
;;so it's not particularly spherical, but it is a useful test when
;;only topology matters
(defun generate-s2-triangulation-of-size (n)
  
  (let* ((triangulation (list (list 4 3 2) (list 1 3 4)
			      (list 1 4 2) (list 2 3 1)))
	 (number-of-triangles 4)
	 (current-triangle () )
	 (rest-of-list     () )
	 (last-used-point 4))
    (while (< number-of-triangles n)
      (setf current-triangle triangulation)
      (while (and (< number-of-triangles n) 
		  (setf rest-of-list (cdr current-triangle)))
	(let* ((triangle (car current-triangle))
	       (p1 (first triangle)) (p2 (second triangle)) 
	       (p3 (third triangle)) (p4 (incf last-used-point)))
	  (setf (car current-triangle) (list p1 p4 p3))
	  (setf (cdr current-triangle) 
		(concatenate 'list (list (list p1 p2 p4) 
					 (list p2 p3 p4)) rest-of-list))
	  (setf current-triangle (cdr (cdr (cdr current-triangle))))
	  (setf number-of-triangles (+ 2 number-of-triangles)))))
    triangulation))

;;function to load a set of triangles from a file.  the triangles are
;;assumed to be stored as a lisp readable list of lists of 3 integers
;;each.  The 3 integers represent point ids.  The point id's can just
;;start at 1, as they will be fixed to match the rest of the spacetime
;;later.
(defun load-triangles-from-file (filename)
  (with-open-file (f filename)
    (read f)))

;;function to connect simplices that have been generated
;;by initialization functions
(defun connect-existing-simplices (&optional initial-spatial-geometry 
				     final-spatial-geometry)
  
  (cond
    ((string= BCTYPE "PERIODIC")

     ;;the code in "simplex.lisp" will automatically modulo high time-stamps
     ;; by NUM-T, wrapping the geometry around at the end
     (for (ts 0 (- NUM-T 1))
	  (connect-simplices-in-sandwich ts (1+ ts) )
	  (connect-simplices-in-adjacent-sandwiches ts (+ ts 1) (+ ts 2)))

     ;; JM: For periodic boundary conditions, the last used point is
     ;; this. However, this relation doesn't hold if the boundaries
     ;; are arbitrary.
     
     ;; JM 2: I have set make-3simplex-v3 to keep track of
     ;; *LAST-USED-POINT* thus it  we don't need this command
     ;; anymore. Commenting it out, but leaving it for consistency reasons.
     ;; (set-last-used-pt (* NUM-T N0-PER-SLICE))

     ;;periodicity gives a 1-to-1 correspondence between spatial
     ;;sheets and 3-d sandwiches, so all per-slice quantities are
     ;;just multiplied by NUM-T
     (set-f-vector (* NUM-T N0-PER-SLICE)    ; N0
		   (* NUM-T N1-SL-PER-SLICE) ; N1-SL
		   (* NUM-T N1-TL-PER-SLICE) ; N1-TL
		   (* NUM-T N2-SL-PER-SLICE) ; N2-SL
		   (* NUM-T N2-TL-PER-SLICE) ; N2-TL
                   ; N3-TL-13 + N3-TL-31
		   (* NUM-T (+ N3-TL-13-PER-SLICE N3-TL-31-PER-SLICE))
		   (* NUM-T N3-TL-22-PER-SLICE)); N3-TL-22
     (set-b-vector 0 0 0 0 0 0)) 
    
        
    ((string= BCTYPE "OPEN")

     ;;set the geometries of the initial and final slices according to
     ;;the optional file parameters, defaulting to single tetrahedra
     ;;when no filenames are given

     ;; JM: The initialial time slice is T=0, and we want to keep it
     ;; that way, so we replace the geometry of the initial time
     ;; slice.
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
     (when final-spatial-geometry
       (let ((existing-triangles ()))
	 ;; Remove simplices in the final sandwich.  

	 ;; JM: Something strange is going on here. David's code
	 ;; deleted simplexes between the final time slice (indexed as
	 ;; NUM-T-1) and the first time slice (if we have periodic
	 ;; boundary conditions) or no time slice at all since the
	 ;; time slice indexed as NUM-T doesn't exist in the periodic
	 ;; case. I don't think David did this on purpose, and I
	 ;; suspect that it was causing a strange disconnect between
	 ;; the initial slice. I am changing this so it replaces time slice
	 ;; NUM-T - 1.

	 ;; Remove any triangles that exist in the initial slice and
	 ;; make sure that the triangles they're connected to in the
	 ;; next time slice are recorded.
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
	 ;; Make simplexes in the first slice based on the initial
	 ;; spatial geometry.	 
	 (map 'list 
	      #'(lambda (x) (make-3simplex-v3 (first x) (second x) (third x)
					      (fourth x) (fifth x) (sixth x)
					      (seventh x)))
	      (triangulate-between-slices existing-triangles 
					     (load-triangles-from-file
					      final-spatial-geometry)
					     (1- NUM-T) NUM-T
					     *LAST-USED-POINT*))))

     ;;connect simplices inside of slices
     (for (ts 0  (1- NUM-T))
        (connect-simplices-in-sandwich ts (1+ ts))
        (connect-simplices-in-adjacent-sandwiches ts (+ ts 1) (+ ts 2)))

     ;; JM: We need to count up the number of points we've used so far
     ;; so that we don't re-use or re-assign points.
     
     ;; JM 2: I have set make-3simplex-v3 to keep track of
     ;; *LAST-USED-POINT* thus we don't need this command
     ;; anymore. Commenting it out, but leaving it for consistency reasons.
     ;; (set-last-used-pt (count-points-in-spacetime))

     ;;set the f-vector.  make sure to consider the arbitrary boundary sheets
     ;;of triangles, computing and adding their contributions separately
     
     ;; JM: Although we can use *-PER-SLICE* in the periodic boundary
     ;; conditions system, this isn't sufficient for when boundary
     ;; conditions are not minimal.
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
      (count-simplices-of-type 2))

     ;;set the initial values of the b-vector (numbers of things on a
     ;;boundary)
     (set-b-vector (count-spacelike-links-at-time NUM-T) ; N1-SL-TOP
		   ; N3-22-TOP
		   (count-simplices-in-sandwich-of-type (1- NUM-T) NUM-T 2)
		   ; N3-31-TOP
		   (count-simplices-in-sandwich-of-type (1- NUM-T) NUM-T 1)
		   (count-spacelike-links-at-time 0) ; N1-SL-BOT
		   ; N3-22-BOT
		   (count-simplices-in-sandwich-of-type 0 1 2)
		   ; N3-31-BOT
		   (count-simplices-in-sandwich-of-type 0 1 3))
		   
)

    (t (error "unrecognized boundary condition type"))))


(defun initialize-S2-triangulation (num-time-slices boundary-conditions 
				    &optional 
				      initial-spatial-geometry
				      final-spatial-geometry)

  (setf NUM-T   num-time-slices
        BCTYPE  (string-upcase boundary-conditions))

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
        ;-----o------------ t = 1
        ;    / \
        ;   /   \
        ;  /     \
        ;-o-------o-------- t = 0
       (dolist (fourpts S2-1/2-31)
	 (make-3simplex-v3 3 (* 2 n) (1+ (* 2 n))                      
			   (+ (* 2 n N0-PER-SLICE) (first fourpts))    
			   (+ (* 2 n N0-PER-SLICE) (second fourpts))   
			   (+ (* 2 n N0-PER-SLICE) (third fourpts))    
			   (+ (* 2 n N0-PER-SLICE) (fourth fourpts)))) 
       ;-----o-----o------ t = 1
       ;    /     / 
       ;   /     /     
       ;  /     /  
       ;-o-----o---------- t = 0
       (dolist (fourpts S2-1/2-22)
	 (make-3simplex-v3 2 (* 2 n) (1+ (* 2 n))                       
			   (+ (* 2 n N0-PER-SLICE) (first fourpts))    
			   (+ (* 2 n N0-PER-SLICE) (second fourpts))   
			   (+ (* 2 n N0-PER-SLICE) (third fourpts))    
			   (+ (* 2 n N0-PER-SLICE) (fourth fourpts)))) 
       ;-o-------o-------- t = 1
       ;  \     /
       ;   \   /
       ;    \ /
       ;-----o------------ t = 0
       (dolist (fourpts S2-1/2-13)
	 (make-3simplex-v3 1 (* 2 n) (1+ (* 2 n))                       
			   (+ (* 2 n N0-PER-SLICE) (first fourpts))    
			   (+ (* 2 n N0-PER-SLICE) (second fourpts))    
			   (+ (* 2 n N0-PER-SLICE) (third fourpts))     
			   (+ (* 2 n N0-PER-SLICE) (fourth fourpts))))  
       ;-o-------o-------- t = 2
       ;  \     /
       ;   \   /
       ;    \ /
       ;-----o------------ t = 1
       (dolist (fourpts S2-1/2-31)
	 (make-3simplex-v3 1 (1+ (* 2 n)) (2+ (* 2 n))                      
			   (+ (* 2 n N0-PER-SLICE) (fourth fourpts))        
			   (+ (* 2 (+ n 1) N0-PER-SLICE) (first fourpts))   
			   (+ (* 2 (+ n 1) N0-PER-SLICE) (second fourpts))  
			   (+ (* 2 (+ n 1) N0-PER-SLICE) (third fourpts)))) 
       ;-----o-----o------ t = 2
       ;      \     \  
       ;       \     \
       ;        \     \
       ;---------o-----o-- t = 1
       (dolist (fourpts S2-1/2-22)
	 (make-3simplex-v3 2 (1+ (* 2 n)) (2+ (* 2 n))                      
			   (+ (* 2 n N0-PER-SLICE) (third fourpts))         
			   (+ (* 2 n N0-PER-SLICE) (fourth fourpts))        
			   (+ (* 2 (+ n 1) N0-PER-SLICE) (first fourpts))   
			   (+ (* 2 (+ n 1) N0-PER-SLICE) (second fourpts))))
       ;-----o------------ t = 2
       ;    / \
       ;   /   \
       ;  /     \
       ;-o-------o-------- t = 1
       (dolist (fourpts S2-1/2-13)
	 (make-3simplex-v3 3 (1+ (* 2 n)) (2+ (* 2 n))                      
			   (+ (* 2 n N0-PER-SLICE) (second fourpts))        
			   (+ (* 2 n N0-PER-SLICE) (third fourpts))         
			   (+ (* 2 n N0-PER-SLICE) (fourth fourpts))        
			   (+ (* 2 (+ n 1) N0-PER-SLICE) (first fourpts)))))

  ;;match triangular faces between 3-simplices
  (connect-existing-simplices initial-spatial-geometry 
			      final-spatial-geometry))

(defun initialize-T2-triangulation (num-time-slices boundary-conditions
				    &optional 
				    initial-spatial-geometry
				    final-spatial-geometry)

  (setf NUM-T   num-time-slices
        BCTYPE  (string-upcase boundary-conditions))

  ;;this code assumes that each toroidal space is being triangulated
  ;;with 32 triangles.  this is not minimal, but it allows for a
  ;;simpler initialization.
  ;;
  ;;  the first torus is triangulated as below
  ;;
  ;;   1   5  9  13   1
  ;;    o--o--o--o--o
  ;;    | /| /| /| /|        each subsequent 
  ;;    |/ |/ |/ |/ |        torus has its 
  ;; 2  o--o--o--o--o 2      points indexed
  ;;    | /| /| /| /|        16 higher than
  ;;    |/ |/ |/ |/ |        the previous
  ;; 3  o--o--o--o--o 3      torus
  ;;    | /| /| /| /|
  ;;    |/ |/ |/ |/ |        
  ;; 4  o--o--o--o--o 4      the TL simplicial filling
  ;;    | /| /| /| /|        alternates between cells
  ;;    |/ |/ |/ |/ |        in a checker-board pattern
  ;;    o--o--o--o--o        so that TL links never cross
  ;;   1   5  9  13   1
  ;;

  ;;set up the numbers of initial simplices according to the T2
  ;;initialization of the geometry
  ; count one corner of each cell
  (defparameter N0-PER-SLICE 16)        
  ; (3 per cell) x (16 cells)
  (defparameter N1-SL-PER-SLICE 48)     
  ; (1 internal + 2 external per cell) x (16 cells)
  (defparameter N1-TL-PER-SLICE 48)     
  ; (2 per cell) x (16 cells)
  (defparameter N2-SL-PER-SLICE 32)     
  ; (2 internal + 4 external per cell) x (16 cells)
  (defparameter N2-TL-PER-SLICE 96)     
  (defparameter N3-TL-13-PER-SLICE 32) ; <---|
  (defparameter N3-TL-22-PER-SLICE 32) ; <---|---| two of each type per cell
  (defparameter N3-TL-31-PER-SLICE 32) ; <---|

  ;;iterate over time slices
  (for (time-slice 0 (1- NUM-T))

    ;;for each time slice, iterate over rows and columns
    ;;of the flattened torus rectangle
    (for (n 0 3) (for (m 0 3)

      (let* (;;define point id's for current time slice
             (p0 (+ 1 (* 16 time-slice) m (* 4 n)))
             (p1 (+ 1 (* 16 time-slice) m (* 4 (mod (1+ n) 4))))
             (p2 (+ 1 (* 16 time-slice) (mod (1+ m) 4) (* 4 n)))
             (p3 (+ 1 (* 16 time-slice) (mod (1+ m) 4) (* 4 (mod (1+ n) 4))))
             ;;define corresponding point id's for next time slice
             (p0-next (+ 16 p0))
             (p1-next (+ 16 p1))
             (p2-next (+ 16 p2))
             (p3-next (+ 16 p3)))

         ;;create simplices for current cell of flattened torus

         ;;enforce checkerboard pattern to prevent crossing of TL links
         (if (evenp (+ n m))

           (progn ;;create 6 3-simplices for the current cell

             (make-3simplex-v3 3 time-slice (1+ time-slice) 
			       p0 p1 p2 p1-next)
             (make-3simplex-v3 2 time-slice (1+ time-slice) 
			       p0 p2 p1-next p2-next)
             (make-3simplex-v3 1 time-slice (1+ time-slice) 
			       p0 p0-next p1-next p2-next)

             (make-3simplex-v3 3 time-slice (1+ time-slice)
			       p1 p2 p3 p1-next)
             (make-3simplex-v3 2 time-slice (1+ time-slice)
			       p2 p3 p1-next p2-next)
             (make-3simplex-v3 1 time-slice (1+ time-slice)
			       p3 p1-next p2-next p3-next))

           (progn ;;otherwise, invert the simplicial breakdown (along
		  ;;the time-direction)

             (make-3simplex-v3 1 time-slice (1+ time-slice) 
			       p1 p0-next p1-next p2-next)
             (make-3simplex-v3 2 time-slice (1+ time-slice) 
			       p1 p2 p0-next p2-next)
             (make-3simplex-v3 3 time-slice (1+ time-slice) 
			       p0 p1 p2 p0-next)

             (make-3simplex-v3 1 time-slice (1+ time-slice) 
			       p1 p1-next p2-next p3-next)
             (make-3simplex-v3 2 time-slice (1+ time-slice)
			       p1 p2 p2-next p3-next)
             (make-3simplex-v3 3 time-slice (1+ time-slice) 
			       p1 p2 p3 p3-next)))))))

  ;;match triangular faces between 3-simplices
  (connect-existing-simplices initial-spatial-geometry 
			      final-spatial-geometry))
            

;; JM: what follows are two methods to increase the volume up to the
;; desired size. The system has yet to thermalize, so algorithm
;; should be irrelevant. The first algorithm is written by David
;; Kamensky. The second is written by Rajesh Kommu.

(defun increase-volume-v1 (target-volume)
  "try volume-increasing moves on random simplices until the desired
   volume is reached."
  (while (< (N3) target-volume)
    ; the range of type-chooser affects 23 / 13 / 31 balance
    (let* ((type-chooser (random 6))
	   (simplex-chooser (random *LAST-USED-3SXID*))
	   (movedata (try-move simplex-chooser 
			       (if (< type-chooser 1) 0 1))))
      (when movedata (2plus1move movedata)))))
	
(defun increase-volume-v2 (target-volume)
  "Use moves to increase the number of simplices until the target
   volume is reached."
  (loop named tv
     do
       
     ;; This is for debugging. Comment it out in general.
;;       (let ((duplicates (list-vals-with-trait 
;;			  #'contains-an-identical-pair *ID->3SIMPLEX* 3)))
;;	 (format t "Duplicate simplices: ~%~S" duplicates))

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
	     (return-from tv)))))


(defun initialize-T-slices-with-V-volume (&key 
					  num-time-slices
					  target-volume
					  spatial-topology
					  boundary-conditions
					  initial-spatial-geometry
					  final-spatial-geometry)

  ;;set global variables according to parameters
  (setf STOPOLOGY  (string-upcase spatial-topology))

  ;;perform initialization based on type of spatial topology
  (cond 
    ((string= STOPOLOGY "S2")
     (initialize-S2-triangulation num-time-slices boundary-conditions 
				  initial-spatial-geometry
				  final-spatial-geometry))
    ((string= STOPOLOGY "T2") 
     (initialize-T2-triangulation num-time-slices boundary-conditions
				  initial-spatial-geometry
				  final-spatial-geometry))
    (t (error "unrecognized spatial topology")))
  
  (format t "final count = ~A~%" (count-simplices-of-all-types))
  (format t "breakdown by location = ~a~%" (count-boundary-vs-bulk))
  (increase-volume-v2 target-volume)
  (format t "~%Plot of 3-simplices as a function of proper time:~%")
  (plot-3-simplices-of-time)
  (format t "~%Plot of spacelike triangles as a function of proper time:~%")
  (plot-spacelike-triangles-of-time)

  (setf N-INIT (N3)))
