;; try-a->b methods return (new4sxids nbors old4sxids old3sxids fvector) 
;; IFF the move can be successfully made 

(defun 3plus1move (sxdata) ;;new4sxdata 4sxnbors old4sxids old3sxids fvector)
  (let ((new4sxids (make-4simplices-in-bulk (first sxdata))))
    (connect-4simplices-within-list new4sxids)
    (connect-4simplices-across-lists new4sxids (second sxdata))
    (remove-4simplices (third sxdata))
    (remove-3simplices (fourth sxdata))
    (update-f-vector (fifth sxdata))))

;; cdt-3plus1-moves.lisp
;;;--------------------------------------------------------------------------[0]
;;; 2->8
;;;
;;;     6      
;;; t+1 -----------
;;;     2 3 4 5 7
;;; t   -----------
;;;     1      
;;; t-1 -----------
;;;
;;; input: either a (1,4) or a (4,1) 4-simplex
;;;
;;; (1,4) + (4,1) -> 4(1,4) + 4(4,1)
;;;
;;; (1|2 3 4 5) 
;;; + 
;;; (2 3 4 5|6) 
;;; ->
;;; (1|2 3 4 7) + (1|2 3 7 5) + (1|2 7 4 5) + (1|7 3 4 5)  
;;; +
;;; (2 3 4 7|6) + (2 3 7 5|6) + (2 7 4 5|6) + (7 3 4 5|6)
;;;--------------------------------------------------------------------------[0]
(defun 2->8-subcomplex (sxid)
  "returns ((14id 41id))"
  (let ((sx (get-4simplex sxid))
	(subcomplex '()))
    (when sx
      (cond ((= 1 (4sx-type sx))
	     (unless (= 0 (nth 0 (4sx-sx4ids sx)))
	       (setf subcomplex (list (list sxid (nth 0 (4sx-sx4ids sx)))))))
	    ((= 4 (4sx-type sx))
	     (unless (= 0 (nth 4 (4sx-sx4ids sx)))
	       (setf subcomplex (list (list (nth 4 (4sx-sx4ids sx)) sxid)))))))
    subcomplex))

(defun try-2->8 (sxid)
  (dolist (curr (2->8-subcomplex sxid))
    (let* ((14id (first curr))
	   (41id (second curr))
	   (14sx (get-4simplex 14id))
	   (41sx (get-4simplex 41id))
	   (14tmlo (4sx-tmlo 14sx)) (14tmhi (4sx-tmhi 14sx))
	   (41tmlo (4sx-tmlo 41sx)) (41tmhi (4sx-tmhi 41sx))
	   (nbors (set-difference (union (4sx-sx4ids 14sx) (4sx-sx4ids 41sx)) (list 0 14id 41id)))
	   (newpt (next-pt))
	   (oldint3sx (nth 0 (4sx-sx3ids 14sx)))
	   (tetra (3sx-points (get-3simplex oldint3sx)))
	   (lopt (nth-point 14sx 0)) (hipt (nth-point 41sx 4))
	   (newsxdata (list (list 1 14tmlo 14tmhi (list lopt newpt 
							(first tetra) (second tetra) (third tetra)))
			    (list 1 14tmlo 14tmhi (list lopt newpt 
							(second tetra) (third tetra) (fourth tetra)))
			    (list 1 14tmlo 14tmhi (list lopt newpt 
							(third tetra) (fourth tetra) (first tetra)))
			    (list 1 14tmlo 14tmhi (list lopt newpt 
							(fourth tetra) (first tetra) (second tetra)))
			    (list 4 41tmlo 41tmhi (list (first tetra) (second tetra) (third tetra) 
							newpt hipt))
			    (list 4 41tmlo 41tmhi (list (second tetra) (third tetra) (fourth tetra) 
							newpt hipt))
			    (list 4 41tmlo 41tmhi (list (third tetra) (fourth tetra) (first tetra) 
							newpt hipt))
			    (list 4 41tmlo 41tmhi (list (fourth tetra) (first tetra) (second tetra) 
							newpt hipt)))))
      (return-from try-2->8 (list newsxdata nbors curr (list oldint3sx) DF28)))))
;;;------------------------------------------------------------------------------------------------------[9]
;;; 8->2
;;;
;;;     6
;;; t+1 -----------
;;;     2 3 4 5 7
;;; t   -----------
;;;     1
;;; t-1 -----------
;;;
;;; input: EITHER a (1,4) OR a (4,1) simplex
;;;
;;; 4(1,4) + 4(4,1) -> (1,4) + (4,1)
;;;
;;; (1|2 3 4 7) + (1|2 3 7 5) + (1|2 7 4 5) + (1|7 3 4 5)  
;;; +
;;; (2 3 4 7|6) + (2 3 7 5|6) + (2 7 4 5|6) + (7 3 4 5|6)
;;; ->
;;; (1|2 3 4 5) 
;;; + 
;;; (2 3 4 5|6) 
;;;------------------------------------------------------------------------------------------------------[9]
(defun 8->2-subcomplex (sxid)
  "returns ((1id1 1id2 1id3 1id4 4id4 4id3 4id2 4id1))"
  (let ((sx nil)(subcmplx nil))
    (when (setf sx (get-4simplex sxid))
      (cond ((= 1 (4sx-type sx))
	     (let ((41id (nth 0 (4sx-sx4ids sx)))(41sx nil))
	       (when (setf 41sx (get-4simplex 41id))
		 (let ((14nbors (neighbors-of-type sx 1)))
		   (unless (< (length 14nbors) 3)
		     (do-tuples/c (id1 id2 id3) 14nbors
		       (let ((sx1 nil) (sx2 nil) (sx3 nil))
			 (when (and (setf sx1 (get-4simplex id1))(setf sx2 (get-4simplex id2))
				    (setf sx3 (get-4simplex id3))(4simplices-connected? id1 id2)
				    (4simplices-connected? id2 id3)(4simplices-connected? id3 id1)
				    (4simplices-connected? (nth 0 (4sx-sx4ids sx1))(nth 0 (4sx-sx4ids sx2)))
				    (4simplices-connected? (nth 0 (4sx-sx4ids sx2))(nth 0 (4sx-sx4ids sx3)))
				    (4simplices-connected? (nth 0 (4sx-sx4ids sx3))(nth 0 (4sx-sx4ids sx1)))
				    (4simplices-connected? (nth 0 (4sx-sx4ids sx1)) 41id)
				    (4simplices-connected? (nth 0 (4sx-sx4ids sx2)) 41id)
				    (4simplices-connected? (nth 0 (4sx-sx4ids sx3)) 41id))
			   (pushnew (list sxid id1 id2 id3 (nth 0 (4sx-sx4ids sx3)) 
					  (nth 0 (4sx-sx4ids sx2))(nth 0 (4sx-sx4ids sx1)) 41id)
				    subcmplx :test #'set-equal?)))))))))
	    ((= 4 (4sx-type sx))
	     (let ((14id (nth 4 (4sx-sx4ids sx)))(14sx nil))
	       (when (setf 14sx (get-4simplex 14id))
		 (let ((41nbors (neighbors-of-type sx 4)))
		   (unless (< (length 41nbors) 3)
		     (do-tuples/c (id1 id2 id3) 41nbors
		       (let ((sx1 nil) (sx2 nil) (sx3 nil))
			 (when (and (setf sx1 (get-4simplex id1)) (setf sx2 (get-4simplex id2))
				    (setf sx3 (get-4simplex id3)) (4simplices-connected? id1 id2)
				    (4simplices-connected? id2 id3) (4simplices-connected? id3 id1)
				    (4simplices-connected? (nth 4 (4sx-sx4ids sx1))(nth 4 (4sx-sx4ids sx2)))
				    (4simplices-connected? (nth 4 (4sx-sx4ids sx2))(nth 4 (4sx-sx4ids sx3)))
				    (4simplices-connected? (nth 4 (4sx-sx4ids sx3))(nth 4 (4sx-sx4ids sx1)))
				    (4simplices-connected? (nth 4 (4sx-sx4ids sx1)) 14id)
				    (4simplices-connected? (nth 4 (4sx-sx4ids sx2)) 14id)
				    (4simplices-connected? (nth 4 (4sx-sx4ids sx3)) 14id))
			   (pushnew (list 14id (nth 4 (4sx-sx4ids sx3)) (nth 4 (4sx-sx4ids sx2))
					  (nth 4 (4sx-sx4ids sx1)) id1 id2 id3 sxid)
				    subcmplx :test #'set-equal?)))))))))))
    subcmplx))

(defun try-8->2 (sxid)
  (dolist (subcx (8->2-subcomplex sxid))
    (let* ((1id1 (first subcx)) (1id2 (second subcx)) (1id3 (third subcx)) (1id4 (fourth subcx))
	   (4id4 (fifth subcx)) (4id3 (sixth subcx)) (4id2 (seventh subcx)) (4id1 (eighth subcx))
	   (1sx1 (get-4simplex 1id1)) (1sx2 (get-4simplex 1id2)) 
	   (1sx3 (get-4simplex 1id3)) (1sx4 (get-4simplex 1id4))
	   (4sx1 (get-4simplex 4id1)) (4sx2 (get-4simplex 4id2))
	   (4sx3 (get-4simplex 4id3)) (4sx4 (get-4simplex 4id4))
	   (14tmlo (4sx-tmlo 1sx1)) (14tmhi (4sx-tmhi 1sx1))
	   (41tmlo (second 4sx1)) (41tmhi (4sx-tmhi 4sx1))
	   (nbors (set-difference (unions (4sx-sx4ids 1sx1) (4sx-sx4ids 1sx2) (4sx-sx4ids 1sx3) 
					  (4sx-sx4ids 1sx4) (4sx-sx4ids 4sx1) (4sx-sx4ids 4sx2) 
					  (4sx-sx4ids 4sx3) (4sx-sx4ids 4sx4))
				   (list 0 1id1 1id2 1id3 1id4 4id1 4id2 4id3 4id4)))
	   (lopt (nth-point 1sx1 0)) (hipt (nth-point 4sx1 4))
	   (oldint3sxs (list (link-id 1id1 1id2)(link-id 1id1 1id3)(link-id 1id1 1id4)(link-id 1id1 4id1)
			     (link-id 1id2 1id3)(link-id 1id2 1id4)(link-id 1id2 4id2)(link-id 1id3 1id4)
			     (link-id 1id3 4id3)(link-id 1id4 4id4)(link-id 4id1 4id2)(link-id 4id1 4id3)
			     (link-id 4id1 4id4)(link-id 4id2 4id3)(link-id 4id2 4id4)(link-id 4id3 4id4)))
	   (newint3sx (list 0 14tmlo 14tmhi
			    (set-difference 
			     (unions (3sx-points (get-3simplex (fourth oldint3sxs)))
				     (3sx-points (get-3simplex (seventh oldint3sxs)))
				     (3sx-points (get-3simplex (ninth oldint3sxs)))
				     (3sx-points (get-3simplex (tenth oldint3sxs))))
			     (intersections (3sx-points (get-3simplex (fourth oldint3sxs)))
					    (3sx-points (get-3simplex (seventh oldint3sxs)))
					    (3sx-points (get-3simplex (ninth oldint3sxs)))
					    (3sx-points (get-3simplex (tenth oldint3sxs)))))))
	   (newsxdata nil))
      (unless (gethash newint3sx *3SIMPLEX->ID*)
	(setf newsxdata (list (list 1 14tmlo 14tmhi (cons lopt (3sx-points newint3sx)))
			      (list 4 41tmlo 41tmhi (append (3sx-points newint3sx)(list hipt)))))
	(return-from try-8->2 (list newsxdata nbors subcx oldint3sxs DF82))))))
;;;------------------------------------------------------------------------------------------------------[1]
;;; 4->6
;;;
;;;     7
;;; t+1 -----------
;;;     2 3 4 5 6
;;; t   -----------
;;;     1
;;; t-1 -----------
;;;
;;; input: EITHER a (1,4) OR a (4,1) simplex
;;;
;;; 2(1,4) + 2(4,1) -> 3(1,4) + 3(4,1)
;;;
;;; (1|2 3 4 5) + (1|3 4 5 6) 
;;; + 
;;; (2 3 4 5|7) + (3 4 5 6|7)
;;; ->
;;; (1|2 3 4 6) + (1|2 4 5 6) + (1|2 3 5 6)
;;; +
;;; (2 3 4 6|7) + (2 4 5 6|7) + (2 3 5 6|7)
;;;
;;;------------------------------------------------------------------------------------------------------[1]
(defun 4->6-subcomplex (sxid)
  "returns ((14id1 14id2 41id2 41id1) (14id1 14id3 41id3 41id1) ... )"
  (let ((28subcmplx (2->8-subcomplex sxid))
	(subcomplex '()))
    (when 28subcmplx
       (let* ((14id1 (first (first 28subcmplx)))
	     (41id1 (second (first 28subcmplx)))
	     (14sx1 (get-4simplex 14id1))
	     (41sx1 (get-4simplex 41id1))
	     (14-nbors-of-14 (neighbors-of-type 14sx1 1))
	     (41-nbors-of-41 (neighbors-of-type 41sx1 4)))
	(dolist (14nborof14 14-nbors-of-14)
	  (let ((14sx (get-4simplex 14nborof14)))
	    (when (and 14sx (find (nth 0 (4sx-sx4ids 14sx)) 41-nbors-of-41))
	      (pushnew (list 14id1 14nborof14 (nth 0 (4sx-sx4ids 14sx)) 41id1) 
		       subcomplex :test #'set-equal?))))))
    subcomplex))

(defun try-4->6 (sxid)
  (dolist (subcx (4->6-subcomplex sxid))
    (let* ((14id1 (first subcx)) (14id2 (second subcx)) (41id2 (third subcx)) (41id1 (fourth subcx))
	   (14sx1 (get-4simplex 14id1)) (14sx2 (get-4simplex 14id2))
	   (41sx1 (get-4simplex 41id1)) (41sx2 (get-4simplex 41id2))
	   (14tmlo (4sx-tmlo 14sx1)) (14tmhi (4sx-tmhi 14sx2))
	   (41tmlo (4sx-tmlo 41sx2)) (41tmhi (4sx-tmhi 41sx1))
	   (nbors (set-difference (unions (4sx-sx4ids 14sx1) (4sx-sx4ids 14sx2) 
					  (4sx-sx4ids 41sx1) (4sx-sx4ids 41sx2))
				  (list 14id1 14id2 41id1 41id2 0)))
	   (oldint3sxids (list (link-id 14id1 14id2) (link-id 41id1 41id2)
			       (link-id 14id1 41id1) (link-id 14id2 41id2)))
	   (pts26 (set-exclusive-or (3sx-points (get-3simplex (third oldint3sxids))) 
				    (3sx-points (get-3simplex (fourth oldint3sxids)))))
	   (pts345 (3sx-hipts (get-3simplex (first oldint3sxids))))
	   (pt1 (3sx-lopts (get-3simplex (first oldint3sxids))))
	   (pt7 (3sx-hipts (get-3simplex (second oldint3sxids))))
	   (newint3sx1 (list 0 14tmlo 14tmhi (append pts26 (firstl pts345) (secondl pts345))))
	   (newint3sx2 (list 0 14tmlo 14tmhi (append pts26 (secondl pts345) (thirdl pts345))))
	   (newint3sx3 (list 0 14tmlo 14tmhi (append pts26 (thirdl pts345) (firstl pts345))))
	   (newint3sx4 (list 1 14tmlo 14tmhi (append pt1 pts26 (firstl pts345))))
	   (newint3sx5 (list 1 14tmlo 14tmhi (append pt1 pts26 (secondl pts345))))
	   (newint3sx6 (list 1 14tmlo 14tmhi (append pt1 pts26 (thirdl pts345))))
	   (newint3sx7 (list 3 41tmlo 41tmhi (append pts26 (firstl pts345) pt7)))
	   (newint3sx8 (list 3 41tmlo 41tmhi (append pts26 (secondl pts345) pt7)))
	   (newint3sx9 (list 3 41tmlo 41tmhi (append pts26 (thirdl pts345) pt7)))
	   (newsxdata nil))
      (unless (or (gethash newint3sx1 *3SIMPLEX->ID*) (gethash newint3sx2 *3SIMPLEX->ID*)
		  (gethash newint3sx3 *3SIMPLEX->ID*) (gethash newint3sx4 *3SIMPLEX->ID*)
		  (gethash newint3sx5 *3SIMPLEX->ID*) (gethash newint3sx6 *3SIMPLEX->ID*)
		  (gethash newint3sx7 *3SIMPLEX->ID*) (gethash newint3sx8 *3SIMPLEX->ID*)
		  (gethash newint3sx9 *3SIMPLEX->ID*))
	(setf newsxdata (list (list 1 14tmlo 14tmhi (append pt1 (3sx-points newint3sx1)))
			      (list 1 14tmlo 14tmhi (append pt1 (3sx-points newint3sx2)))
			      (list 1 14tmlo 14tmhi (append pt1 (3sx-points newint3sx3)))
			      (list 4 41tmlo 41tmhi (append (3sx-points newint3sx1) pt7))
			      (list 4 41tmlo 41tmhi (append (3sx-points newint3sx2) pt7))
			      (list 4 41tmlo 41tmhi (append (3sx-points newint3sx3) pt7))))
	(return-from try-4->6 (list newsxdata nbors subcx oldint3sxids DF46))))))
;;;------------------------------------------------------------------------------------------------------[8]
;;; 6->4
;;;
;;;     7
;;; t+1 -----------
;;;     2 3 4 5 6
;;; t   -----------
;;;     1
;;; t-1 -----------
;;;
;;; input: EITHER a (1,4) OR a (4,1) simplex
;;;
;;; 3(1,4) + 3(4,1) -> 2(1,4) + 2(4,1)
;;;
;;; (1|2 3 4 6) + (1|2 4 5 6) + (1|2 3 5 6)
;;; +
;;; (2 3 4 6|7) + (2 4 5 6|7) + (2 3 5 6|7)
;;; ->
;;; (1|2 3 4 5) + (1|3 4 5 6) 
;;; + 
;;; (2 3 4 5|7) + (3 4 5 6|7)
;;;------------------------------------------------------------------------------------------------------[8]
(defun 6->4-subcomplex (sxid)
  "returns ((1id1 1id2 1id3 4id3 4id2 4id1))"
  (let ((sx nil)(subcmplx nil))
    (when (setf sx (get-4simplex sxid))
      (cond ((= 1 (4sx-type sx))
	     (let ((41id (nth 0 (4sx-sx4ids sx)))(41sx nil))
	       (when (setf 41sx (get-4simplex 41id))
		 (let ((14nbors (neighbors-of-type sx 1)))
		   (unless (< (length 14nbors) 2)
		     (do-tuples/c (id1 id2) 14nbors
		       (let ((sx1 nil) (sx2 nil))
			 (when (and (setf sx1 (get-4simplex id1))(setf sx2 (get-4simplex id2))
				    (4simplices-connected? id1 id2)
				    (4simplices-connected? (nth 0 (4sx-sx4ids sx1))(nth 0 (4sx-sx4ids sx2)))
				    (4simplices-connected? (nth 0 (4sx-sx4ids sx1)) 41id)
				    (4simplices-connected? (nth 0 (4sx-sx4ids sx2)) 41id))
			   (pushnew (list sxid id1 id2 
					  (nth 0 (4sx-sx4ids sx2))(nth 0 (4sx-sx4ids sx1)) 41id)
				    subcmplx :test #'set-equal?)))))))))
	    ((= 4 (4sx-type sx))
	     (let ((14id (nth 4 (4sx-sx4ids sx)))(14sx nil))
	       (when (setf 14sx (get-4simplex 14id))
		 (let ((41nbors (neighbors-of-type sx 4)))
		   (unless (< (length 41nbors) 2)
		     (do-tuples/c (id1 id2) 41nbors
		       (let ((sx1 nil) (sx2 nil))
			 (when (and (setf sx1 (get-4simplex id1)) (setf sx2 (get-4simplex id2))
				    (4simplices-connected? id1 id2)
				    (4simplices-connected? (nth 4 (4sx-sx4ids sx1))(nth 4 (4sx-sx4ids sx2)))
				    (4simplices-connected? (nth 4 (4sx-sx4ids sx1)) 14id)
				    (4simplices-connected? (nth 4 (4sx-sx4ids sx2)) 14id))
			   (pushnew (list 14id (nth 4 (4sx-sx4ids sx2)) (nth 4 (4sx-sx4ids sx1)) 
					  id1 id2 sxid)
				    subcmplx :test #'set-equal?)))))))))))
    subcmplx))

(defun try-6->4 (sxid)
  (dolist (subcx (6->4-subcomplex sxid))
    (let* ((14id1 (first subcx)) (14id2 (second subcx)) (14id3 (third subcx))
	   (41id3 (fourth subcx)) (41id2 (fifth subcx)) (41id1 (sixth subcx))
	   (14sx1 (get-4simplex 14id1)) (14sx2 (get-4simplex 14id2)) (14sx3 (get-4simplex 14id3))
	   (41sx1 (get-4simplex 41id1)) (41sx2 (get-4simplex 41id2)) (41sx3 (get-4simplex 41id3))
	   (14tmlo (4sx-tmlo 14sx1)) (14tmhi (4sx-tmhi 14sx2))
	   (41tmlo (4sx-tmlo 41sx2)) (41tmhi (4sx-tmhi 41sx1))
	   (nbors (set-difference (unions (4sx-sx4ids 14sx1) (4sx-sx4ids 14sx2) (4sx-sx4ids 14sx3)
					  (4sx-sx4ids 41sx1) (4sx-sx4ids 41sx2) (4sx-sx4ids 41sx3))
				  (list 14id1 14id2 14id3 41id1 41id2 41id3 0)))
	   (oldint3sxids (list (link-id 14id1 41id1) (link-id 14id2 41id2) (link-id 14id3 41id3) 
			       (link-id 14id1 14id3) (link-id 14id1 14id2) (link-id 14id2 14id3)
			       (link-id 41id1 41id3) (link-id 41id1 41id2) (link-id 41id2 41id3)))
	   (pts26 (intersections (3sx-points (get-3simplex (first oldint3sxids)))
				 (3sx-points (get-3simplex (second oldint3sxids)))
				 (3sx-points (get-3simplex (third oldint3sxids)))))
	   (pts345 (set-difference (unions (3sx-points (get-3simplex (first oldint3sxids)))
					   (3sx-points (get-3simplex (second oldint3sxids)))
					   (3sx-points (get-3simplex (third oldint3sxids))))
				   pts26))
	   (pt1 (4sx-lopts 14sx1)) (pt7 (4sx-hipts 41sx1))
	   (newint3sx1 (list 0 14tmlo 14tmhi (append (firstl pts26) pts345)))
	   (newint3sx2 (list 0 14tmlo 14tmhi (append (secondl pts26) pts345)))
	   (newint3sx3 (list 1 14tmlo 14tmhi (append pt1 pts345)))
	   (newint3sx4 (list 3 41tmlo 41tmhi (append pts345 pt7)))
	   (newsxdata nil))
      (unless (or (gethash newint3sx1 *3SIMPLEX->ID*) (gethash newint3sx2 *3SIMPLEX->ID*)
		  (gethash newint3sx3 *3SIMPLEX->ID*) (gethash newint3sx4 *3SIMPLEX->ID*))
	(setf newsxdata (list (list 1 14tmlo 14tmhi (append pt1 (3sx-points newint3sx1)))
			      (list 1 14tmlo 14tmhi (append pt1 (3sx-points newint3sx2)))
			      (list 4 41tmlo 41tmhi (append (3sx-points newint3sx1) pt7))
			      (list 4 41tmlo 41tmhi (append (3sx-points newint3sx2) pt7))))
	(return-from try-6->4 (list newsxdata nbors subcx oldint3sxids DF64))))))
;;;------------------------------------------------------------------------------------------------------[2]
;;; 2->4-v1
;;;
;;;     3 4 5 6           5 6
;;; t+1 -----------   t+1 -----------
;;;     1 2               1 2 3 4
;;; t   -----------     t -----------
;;;
;;; input: EITHER A (1,4) or a (4,1) or a (2,3) or a (3,2) simplex
;;;
;;; 1(1,4) + 1(2,3) -> 1(1,4) + 3(2,3)
;;; 1(4,1) + 1(3,2) -> 1(4,1) + 3(3,2)
;;;
;;; (1|3 4 5 6)
;;; +
;;; (1 2|3 4 5)
;;; ->
;;; (2|3 4 5 6)
;;; +
;;; (1 2|4 5 6) + (1 2|5 6 3) + (1 2|6 3 4)
;;;
;;; or
;;;
;;; (3 4 5 6|1)
;;; +
;;; (3 4 5|1 2)
;;; ->
;;; (3 4 5 6|2)
;;; +
;;; (4 5 6|1 2) + (5 3 6|1 2) + (3 4 6|1 2) 
;;;------------------------------------------------------------------------------------------------------[2]
(defun 2->4-v1-subcomplex (sxid)
  "returns ((12or43 1or4id1 2or3id1) (12or43 1or4id2 2or3id2) ... ) etc."
  (let ((sx '())
	(subcomplex '()))
    (when (setf sx (get-4simplex sxid))
      (cond ((= 1 (first sx))
	     (dolist (23nbor (neighbors-of-type sx 2))
	       (pushnew (list 12 sxid 23nbor) subcomplex :test #'set-equal?)))
	    ((= 2 (first sx))
	     (dolist (14nbor (neighbors-of-type sx 1))
	       (pushnew (list 12 14nbor sxid) subcomplex :test #'set-equal?)))
	    ((= 3 (first sx))
	     (dolist (41nbor (neighbors-of-type sx 4))
	       (pushnew (list 43 41nbor sxid) subcomplex :test #'set-equal?)))
	    ((= 4 (first sx))
	     (dolist (32nbor (neighbors-of-type sx 3))
	       (pushnew (list 43 sxid 32nbor) subcomplex :test #'set-equal?)))))
    subcomplex))

(defun 2->4-v1-internal-12 (14id 23id)
  (let ((14sx nil) (23sx nil))
    (when (and (setf 14sx (get-4simplex 14id)) (setf 23sx (get-4simplex 23id)))
      (let* ((oldint3sxid (link-id 14id 23id))
	     (pts345 (3sx-hipts (get-3simplex oldint3sxid)))
	     (pt6 (set-difference (4sx-hipts 14sx) pts345))
	     (pts12 (4sx-lopts 23sx))
	     (pt2 (set-difference pts12 (4sx-lopts 14sx)))
	     (tmlo (4sx-tmlo 14sx)) (tmhi (4sx-tmhi 23sx))
	     (nbors (set-difference (union (4sx-sx4ids 14sx) (4sx-sx4ids 23sx))
				    (list 14id 23id 0)))
	     (newint3sx1 (list 1 tmlo tmhi (append pt2 (firstl pts345) (secondl pts345) pt6)))
	     (newint3sx2 (list 1 tmlo tmhi (append pt2 (secondl pts345) (thirdl pts345) pt6)))
	     (newint3sx3 (list 1 tmlo tmhi (append pt2 (thirdl pts345) (firstl pts345) pt6)))
	     (newint3sx4 (list 2 tmlo tmhi (append pts12 (firstl pts345) pt6)))
	     (newint3sx5 (list 2 tmlo tmhi (append pts12 (secondl pts345) pt6)))
	     (newint3sx6 (list 2 tmlo tmhi (append pts12 (thirdl pts345) pt6)))
	     (newsxdata nil))
	(unless (or (gethash newint3sx1 *3SIMPLEX->ID*) (gethash newint3sx2 *3SIMPLEX->ID*)
		    (gethash newint3sx3 *3SIMPLEX->ID*) (gethash newint3sx4 *3SIMPLEX->ID*)
		    (gethash newint3sx5 *3SIMPLEX->ID*) (gethash newint3sx6 *3SIMPLEX->ID*))
	  (setf newsxdata (list (list 1 tmlo tmhi (append pt2 (4sx-hipts 14sx)))
				(list 2 tmlo tmhi (append pts12 (3sx-hipts newint3sx1)))
				(list 2 tmlo tmhi (append pts12 (3sx-hipts newint3sx2)))
				(list 2 tmlo tmhi (append pts12 (3sx-hipts newint3sx3)))))
	  (return-from 2->4-v1-internal-12 (list newsxdata nbors (list 14id 23id) (list oldint3sxid) DF24)))))))

(defun 2->4-v1-internal-43 (41id 32id)
  (let ((41sx nil) (32sx nil))
    (when (and (setf 41sx (get-4simplex 41id)) (setf 32sx (get-4simplex 32id)))
      (let* ((oldint3sxid (link-id 41id 32id))
	     (pts345 (3sx-lopts (get-3simplex oldint3sxid)))
	     (pt6 (set-difference (4sx-lopts 41sx) pts345))
	     (pts12 (4sx-hipts 32sx))
	     (pt2 (set-difference pts12 (4sx-hipts 41sx)))
	     (tmlo (4sx-tmlo 41sx)) (tmhi (4sx-tmhi 32sx))
	     (nbors (set-difference (union (4sx-sx4ids 41sx) (4sx-sx4ids 32sx))
				    (list 41id 32id 0)))
	     (newint3sx1 (list 3 tmlo tmhi (append (firstl pts345) (secondl pts345) pt6 pt2)))
	     (newint3sx2 (list 3 tmlo tmhi (append (secondl pts345) (thirdl pts345) pt6 pt2)))
	     (newint3sx3 (list 3 tmlo tmhi (append (thirdl pts345) (firstl pts345) pt6 pt2)))
	     (newint3sx4 (list 2 tmlo tmhi (append (firstl pts345) pt6 pts12)))
	     (newint3sx5 (list 2 tmlo tmhi (append (secondl pts345) pt6 pts12)))
	     (newint3sx6 (list 2 tmlo tmhi (append (thirdl pts345) pt6 pts12)))
	     (newsxdata nil))
	(unless (or (gethash newint3sx1 *3SIMPLEX->ID*) (gethash newint3sx2 *3SIMPLEX->ID*)
		    (gethash newint3sx3 *3SIMPLEX->ID*) (gethash newint3sx4 *3SIMPLEX->ID*)
		    (gethash newint3sx5 *3SIMPLEX->ID*) (gethash newint3sx6 *3SIMPLEX->ID*))
	  (setf newsxdata (list (list 4 tmlo tmhi (append (4sx-lopts 41sx) pt2))
				(list 3 tmlo tmhi (append (3sx-lopts newint3sx1) pts12))
				(list 3 tmlo tmhi (append (3sx-lopts newint3sx2) pts12))
				(list 3 tmlo tmhi (append (3sx-lopts newint3sx3) pts12))))
	  (return-from 2->4-v1-internal-43 (list newsxdata nbors (list 41id 32id) (list oldint3sxid) DF24)))))))

(defun try-2->4-v1 (sxid)
  (let ((subcmplx (2->4-v1-subcomplex sxid))
	(movedata nil))
    (unless (null subcmplx)
      (dolist (curr subcmplx)
	(cond ((= 12 (first curr))
	       (setf movedata (2->4-v1-internal-12 (second curr) (third curr)))
	       (when movedata
		 (return-from try-2->4-v1 movedata)))
	      ((= 43 (first curr))
	       (setf movedata (2->4-v1-internal-43 (second curr) (third curr)))
	       (when movedata
		 (return-from try-2->4-v1 movedata))))))))
;;;------------------------------------------------------------------------------------------------------[7]
;;; 4->2-v1
;;;
;;;     3 4 5 6
;;; t+1 -----------
;;;     1 2
;;; t   -----------
;;;
;;; input: EITHER a (1,4) OR a (4,1) OR a (2,3) OR a (3,2) simplex
;;;
;;; 1(1,4) + 3(2,3) -> 1(1,4) + 1(2,3)
;;; 1(4,1) + 3(3,2) -> 1(4,1) + 1(3,2) 
;;;
;;; (1 2|3 4 5) + (1 2|3 4 6) + (1 2|3 5 6)
;;; +
;;; (2|3 4 5 6)
;;; ->
;;; (1|3 4 5 6)
;;; +
;;; (1 2|4 5 6)
;;;------------------------------------------------------------------------------------------------------[7]
(defun 4->2-v1-subcomplex (sxid)
  "returns ((1222or4333 1or4id 2or3id1 2or3id2 2or3id3)...)"
  (let ((sx '())
	(subcmplx '()))
    (when (setf sx (get-4simplex sxid))
      (cond ((= 1 (4sx-type sx))
	     (let ((23nbors (neighbors-of-type sx 2)))
	       (unless (< (length 23nbors) 3)
		 (do-tuples/c (id1 id2 id3) 23nbors
			      (let ((sx1 nil) (sx2 nil) (sx3 nil))
				(when (and (setf sx1 (get-4simplex id1)) (setf sx2 (get-4simplex id2))
					   (setf sx3 (get-4simplex id3)) (4simplices-connected? id1 id2)
					   (4simplices-connected? id2 id3) (4simplices-connected? id3 id1))
				  (pushnew (list 1222 sxid id1 id2 id3) subcmplx :test #'set-equal?)))))))
	    ((= 2 (4sx-type sx))
	     (let ((23nbors (neighbors-of-type sx 2))
		   (14nbors (neighbors-of-type sx 1)))
	       (unless (or (< (length 23nbors) 2) (< (length 14nbors) 1))
		 (do-tuples/c (id1 id2) 23nbors
			      (let ((sx1 nil) (sx2 nil))
				(when (and (setf sx1 (get-4simplex id1)) (setf sx2 (get-4simplex id2))
					   (4simplices-connected? id1 id2))
				  (dolist (14nbor 14nbors)
				    (when (and (4simplices-connected? id1 14nbor)
					       (4simplices-connected? id2 14nbor))
				      (pushnew (list 1222 14nbor sxid id1 id2) subcmplx :test #'set-equal?)))))))))
	    ((= 3 (4sx-type sx))
	     (let ((32nbors (neighbors-of-type sx 3))
		   (41nbors (neighbors-of-type sx 4)))
	       (unless (or (< (length 32nbors) 2) (< (length 41nbors) 1))
		 (do-tuples/c (id1 id2) 32nbors
			      (let ((sx1 nil) (sx2 nil))
				(when (and (setf sx1 (get-4simplex id1)) (setf sx2 (get-4simplex id2))
					   (4simplices-connected? id1 id2))
				  (dolist (41nbor 41nbors)
				    (when (and (4simplices-connected? id1 41nbor)
					       (4simplices-connected? id2 41nbor))
				      (pushnew (list 4333 41nbor sxid id1 id2) subcmplx :test #'set-equal?)))))))))
	    ((= 4 (4sx-type sx))
	     (let ((32nbors (neighbors-of-type sx 3)))
	       (unless (< (length 32nbors) 3)
		 (do-tuples/c (id1 id2 id3) 32nbors
			      (let ((sx1 nil) (sx2 nil) (sx3 nil))
				(when (and (setf sx1 (get-4simplex id1)) (setf sx2 (get-4simplex id2))
					   (setf sx3 (get-4simplex id3)) (4simplices-connected? id1 id2)
					   (4simplices-connected? id2 id3) (4simplices-connected? id3 id1))
				  (pushnew (list 4333 sxid id1 id2 id3) subcmplx :test #'set-equal?)))))))))
    subcmplx))

(defun 4->2-v1-internal-1222 (14id 23id1 23id2 23id3)
  (let ((14sx nil) (23sx1 nil) (23sx2 nil) (23sx3 nil))
    (when (and (setf 14sx (get-4simplex 14id)) (setf 23sx1 (get-4simplex 23id1))
	       (setf 23sx2 (get-4simplex 23id2)) (setf 23sx3 (get-4simplex 23id3)))
      (let* ((oldint3sxids (list (link-id 14id 23id1) (link-id 14id 23id2) (link-id 14id 23id3)
				 (link-id 23id1 23id2) (link-id 23id2 23id3) (link-id 23id3 23id1)))
	     (pt3 (intersections (4sx-hipts 23sx1) (4sx-hipts 23sx2) (4sx-hipts 23sx3)))
	     (pts456 (set-difference (4sx-hipts 14sx) pt3))
	     (pts12 (4sx-lopts 23sx1))
	     (pt1 (set-difference pts12 (4sx-lopts 14sx)))
	     (tmlo (4sx-tmlo 14sx)) (tmhi (4sx-tmhi 23sx1))
	     (nbors (set-difference (unions (4sx-sx4ids 14sx)
					    (4sx-sx4ids 23sx1) (4sx-sx4ids 23sx2) (4sx-sx4ids 23sx3))
				    (list 14id 23id1 23id2 23id3 0)))
	     (newint3sx (list 1 tmlo tmhi (append pt1 pts456)))
	     (newsxdata nil))
	(unless (gethash newint3sx *3SIMPLEX->ID*)
	  (setf newsxdata (list (list 1 tmlo tmhi (append pt1 (4sx-hipts 14sx)))
				(list 2 tmlo tmhi (append pts12 (3sx-hipts newint3sx)))))
	  (return-from 4->2-v1-internal-1222 (list newsxdata nbors (list 14id 23id1 23id2 23id3)
						   oldint3sxids DF42)))))))

(defun 4->2-v1-internal-4333 (41id 32id1 32id2 32id3)
  (setf CURRENT-MOVE-IDENTIFIER "4->2-v1-4333") (incf CURRENT-MOVE-NUMBER)
  (let ((41sx nil) (32sx1 nil) (32sx2 nil) (32sx3 nil))
    (when (and (setf 41sx (get-4simplex 41id)) (setf 32sx1 (get-4simplex 32id1))
	       (setf 32sx2 (get-4simplex 32id2)) (setf 32sx3 (get-4simplex 32id3)))
      (let* ((oldint3sxids (list (link-id 41id 32id1) (link-id 41id 32id2) (link-id 41id 32id3)
				 (link-id 32id1 32id2) (link-id 32id2 32id3) (link-id 32id3 32id1)))
	     (pt3 (intersections (4sx-lopts 32sx1) (4sx-lopts 32sx2) (4sx-lopts 32sx3)))
	     (pts456 (set-difference (4sx-lopts 41sx) pt3))
	     (pts12 (4sx-hipts 32sx1))
	     (pt1 (set-difference pts12 (4sx-hipts 41sx)))
	     (tmlo (4sx-tmlo 41sx)) (tmhi (4sx-tmhi 32sx1))
	     (nbors (set-difference (unions (4sx-sx4ids 41sx)
					    (4sx-sx4ids 32sx1) (4sx-sx4ids 32sx2) (4sx-sx4ids 32sx3))
				    (list 41id 32id1 32id2 32id3 0)))
	     (newint3sx (list 3 tmlo tmhi (append pts456 pt1)))
	     (newsxdata nil))
	(unless (gethash newint3sx *3SIMPLEX->ID*)
	  (setf newsxdata (list (list 4 tmlo tmhi (append (4sx-lopts 41sx) pt1))
				(list 3 tmlo tmhi (append (3sx-lopts newint3sx) pts12))))
	  (return-from 4->2-v1-internal-4333 (list newsxdata nbors (list 41id 32id1 32id2 32id3)
						   oldint3sxids DF42)))))))

(defun try-4->2-v1 (sxid)
  (let ((subcmplx (4->2-v1-subcomplex sxid))
	(movedata nil))
    (unless (null subcmplx)
      (dolist (curr subcmplx)
	(cond ((= 1222 (first curr))
	       (setf movedata (4->2-v1-internal-1222 (second curr) (third curr) (fourth curr) (fifth curr)))
	       (when movedata
		 (return-from try-4->2-v1 movedata)))
	      ((= 4333 (first curr))
	       (setf movedata (4->2-v1-internal-4333 (second curr) (third curr) (fourth curr) (fifth curr)))
	       (when movedata
		 (return-from try-4->2-v1 movedata))))))))
;;;------------------------------------------------------------------------------------------------------[3]
;;; 2->4-v2
;;;
;;;     4 5 6
;;; t+1 -----------
;;;     1 2 3
;;; t   -----------
;;;
;;; INPUT: EITHER a (2,3) OR a (3,2) simplex
;;;
;;; 1(2,3) + 1(3,2) -> 2(2,3) + 2(3,2)
;;;
;;; (1 2| 4 5 6)
;;; +
;;; (1 2 3|4 5)
;;; ->
;;; (1 3|4 5 6) + (2 3|4 5 6)
;;; +
;;; (1 2 3|4 6) + (1 2 3|5 6)
;;;------------------------------------------------------------------------------------------------------[3]
(defun 2->4-v2-subcomplex (sxid)
  "returns ((23id 32id1) (23id 32id2) (23id 32id3)) or ((23id1 32id) (23id2 32id) (23id3 32id))"
  (let ((sx '())
	(subcmplx '()))
    (when (setf sx (get-4simplex sxid))
      (cond ((= 2 (4sx-type sx))
	     (dolist (32nbor (neighbors-of-type sx 3))
	       (pushnew (list sxid 32nbor) subcmplx :test #'set-equal?)))
	    ((= 3 (4sx-type sx))
	     (dolist (23nbor (neighbors-of-type sx 2))
	       (pushnew (list 23nbor sxid) subcmplx :test #'set-equal?)))))
    subcmplx))

(defun try-2->4-v2 (sxid)
  (dolist (subcx (2->4-v2-subcomplex sxid))
    (let* ((23id (first subcx)) (32id (second subcx))
	   (23sx (get-4simplex 23id)) (32sx (get-4simplex 32id))
	   (tmlo (4sx-tmlo 23sx)) (tmhi (4sx-tmhi 32sx))
	   (nbors (set-difference (union (4sx-sx4ids 23sx) (4sx-sx4ids 32sx)) 
				  (list 23id 32id 0)))
	   (oldint3sxid (link-id 23id 32id))
	   (pts12 (4sx-lopts 23sx)) (pts456 (4sx-hipts 23sx))
	   (pts123 (4sx-lopts 32sx)) (pts45 (4sx-hipts 32sx))
	   (pt3 (set-difference pts123 pts12)) (pt6 (set-difference pts456 pts45))
	   (newint3sx1 (list 1 tmlo tmhi (append pt3 pts456)))
	   (newint3sx2 (list 3 tmlo tmhi (append pts123 pt6)))
	   (newint3sx3 (list 2 tmlo tmhi (append (firstl pts12) pt3 (firstl pts45) pt6)))
	   (newint3sx4 (list 2 tmlo tmhi (append (firstl pts12) pt3 (secondl pts45) pt6)))
	   (newint3sx5 (list 2 tmlo tmhi (append (secondl pts12) pt3 (firstl pts45) pt6)))
	   (newint3sx6 (list 2 tmlo tmhi (append (secondl pts12) pt3 (secondl pts45) pt6)))
	   (newsxdata nil))
      (unless (or (gethash newint3sx1 *3SIMPLEX->ID*) (gethash newint3sx2 *3SIMPLEX->ID*)
		  (gethash newint3sx3 *3SIMPLEX->ID*) (gethash newint3sx4 *3SIMPLEX->ID*)
		  (gethash newint3sx5 *3SIMPLEX->ID*) (gethash newint3sx6 *3SIMPLEX->ID*))
	(setf newsxdata (list (list 2 tmlo tmhi (append (3sx-lopts newint3sx3) pts456))
			      (list 2 tmlo tmhi (append (3sx-lopts newint3sx5) pts456))
			      (list 3 tmlo tmhi (append pts123 (3sx-hipts newint3sx3)))
			      (list 3 tmlo tmhi (append pts123 (3sx-hipts newint3sx4)))))
	(return-from try-2->4-v2 (list newsxdata nbors subcx (list oldint3sxid) DF24)))))) 
;;;------------------------------------------------------------------------------------------------------[6]
;;; 4->2-v2
;;;
;;;     4 5 6
;;; t+1 -----------
;;;     1 2 3
;;; t   -----------
;;;
;;; INPUT: EITHER a (2,3) OR a (3,2) simplex
;;;
;;; 2(2,3) + 2(3,2) -> 1(2,3) + 1(3,2)
;;;
;;; (1 3|4 5 6) + (2 3|4 5 6)
;;; +
;;; (1 2 3|4 6) + (1 2 3|5 6)
;;; ->
;;; (1 2| 4 5 6)
;;; +
;;; (1 2 3|4 5)
;;;------------------------------------------------------------------------------------------------------[6]
(defun 4->2-v2-subcomplex (sxid)
  "returns ((23id 23id 32id 32id)...)"
  (let ((sx '())
	(subcmplx '()))
    (when (setf sx (get-4simplex sxid))
      (cond ((= 2 (4sx-type sx))
	     (let ((23nbors (neighbors-of-type sx 2))
		   (32nbors (neighbors-of-type sx 3)))
	       (unless (or (< (length 23nbors) 1) (< (length 32nbors) 2))
		 (do-tuples/c (id1 id2) 32nbors
			      (let ((sx1 nil) (sx2 nil))
				(when (and (setf sx1 (get-4simplex id1)) (setf sx2 (get-4simplex id2))
					   (4simplices-connected? id1 id2))
				  (dolist (23nbor 23nbors)
				    (when (and (4simplices-connected? id1 23nbor)
					       (4simplices-connected? id2 23nbor))
				      (pushnew (list sxid 23nbor id1 id2) subcmplx 
					       :test #'set-equal?)))))))))
	    ((= 3 (4sx-type sx))
	     (let ((23nbors (neighbors-of-type sx 2))
		   (32nbors (neighbors-of-type sx 3)))
	       (unless (or (< (length 23nbors) 2) (< (length 32nbors) 1))
		 (do-tuples/c (id1 id2) 23nbors
			      (let ((sx1 nil) (sx2 nil))
				(when (and (setf sx1 (get-4simplex id1)) (setf sx2 (get-4simplex id2))
					   (4simplices-connected? id1 id2))
				  (dolist (32nbor 32nbors)
				    (when (and (4simplices-connected? id1 32nbor)
					       (4simplices-connected? id2 32nbor))
				      (pushnew (list id1 id2 32nbor sxid) subcmplx 
					       :test #'set-equal?)))))))))))
    subcmplx))

(defun try-4->2-v2 (sxid)
  (dolist (subcx (4->2-v2-subcomplex sxid))
    (let* ((23id1 (first subcx)) (23id2 (second subcx))
	   (32id1 (third subcx)) (32id2 (fourth subcx))
	   (23sx1 (get-4simplex 23id1)) (23sx2 (get-4simplex 23id2))
	   (32sx1 (get-4simplex 32id1)) (32sx2 (get-4simplex 32id2))
	   (tmlo (4sx-tmlo 23sx1)) (tmhi (4sx-tmhi 32sx2))
	   (nbors (set-difference (unions (4sx-sx4ids 23sx1) (4sx-sx4ids 23sx2)
					  (4sx-sx4ids 32sx1) (4sx-sx4ids 32sx2))
				  (list 23id1 23id2 32id1 32id2 0)))
	   (oldint3sxids (list (link-id 23id1 23id2) (link-id 32id1 32id2) (link-id 23id1 32id1) 
			       (link-id 23id1 32id2) (link-id 23id2 32id1) (link-id 23id2 32id2)))
	   (pts13 (4sx-lopts 23sx1)) (pts23 (4sx-lopts 23sx2)) (pts456 (4sx-hipts 23sx1))
	   (pts123 (4sx-lopts 32sx1)) (pts46 (4sx-hipts 32sx1)) (pts56 (4sx-hipts 32sx2))
	   (newint3sx (list 2 tmlo tmhi (append (set-exclusive-or pts13 pts23) 
						(set-exclusive-or pts46 pts56))))
	   (newsxdata nil))
      (unless (gethash newint3sx *3SIMPLEX->ID*)
	(setf newsxdata (list (list 2 tmlo tmhi (append (3sx-lopts newint3sx) pts456))
			      (list 3 tmlo tmhi (append pts123 (3sx-hipts newint3sx)))))
	(return-from try-4->2-v2 (list newsxdata nbors subcx oldint3sxids DF42))))))
;;;------------------------------------------------------------------------------------------------------[4]
;;; 3->3-v1
;;;
;;;     3 4 5 6           5 6
;;; t+1 -----------   t+1 ----------
;;;     1 2               1 2 3 4
;;; t   -----------     t ----------
;;;
;;; INPUT: either a (1,4) or a (4,1) or a (2,3) or a (3,2) simplex
;;;
;;; 1(1,4) + 2(2,3) -> 1(1,4) + 2(2,3)
;;; 1(4,1) + 2(3,2) -> 1(4,1) + 2(3,2)
;;;
;;; (1|3 4 5 6)
;;; +
;;; (1 2|3 4 5) + (1 2|3 4 6)
;;; ->
;;; (2|3 4 5 6)
;;; +
;;; (1 2|3 5 6) + (1 2|4 5 6)
;;; 
;;; (1 2 3 4|5)
;;; +
;;; (1 2 3|5 6) + (1 2 4|5 6)
;;; ->
;;; (1 2 3 4|6)
;;; +
;;; (2 3 4|5 6) + (1 3 4|5 6)
;;;------------------------------------------------------------------------------------------------------[4]
(defun 3->3-v1-subcomplex (sxid)
  "returns ((122or433 14or41id 23or32id1 23or32id2)...)"
  (let ((sx '())
	(subcmplx '()))
    (when (setf sx (get-4simplex sxid))
      (cond ((= 1 (4sx-type sx))
	     (let ((23nbors (neighbors-of-type sx 2)))
	       (unless (< (length 23nbors) 2)
		 (do-tuples/c (id1 id2) 23nbors
			      (let ((sx1 nil) (sx2 nil))
				(when (and (setf sx1 (get-4simplex id1)) (setf sx2 (get-4simplex id2))
					   (4simplices-connected? id1 id2))
				  (pushnew (list 122 sxid id1 id2) subcmplx :test #'set-equal?)))))))
	    ((= 2 (4sx-type sx))
	     (let ((23nbors (neighbors-of-type sx 2))
		   (14nbors (neighbors-of-type sx 1)))
	       (dolist (23nbor 23nbors)
		 (dolist (14nbor 14nbors)
		   (let ((sx1 nil) (sx2 nil))
		     (when (and (setf sx1 (get-4simplex 23nbor)) (setf sx2 (get-4simplex 14nbor))
				(4simplices-connected? 23nbor 14nbor))
		       (pushnew (list 122 14nbor sxid 23nbor) subcmplx :test #'set-equal?)))))))
	    ((= 3 (4sx-type sx))
	     (let ((32nbors (neighbors-of-type sx 3))
		   (41nbors (neighbors-of-type sx 4)))
	       (dolist (32nbor 32nbors)
		 (dolist (41nbor 41nbors)
		   (let ((sx1 nil) (sx2 nil))
		     (when (and (setf sx1 (get-4simplex 32nbor)) (setf sx2 (get-4simplex 41nbor))
				(4simplices-connected? 32nbor 41nbor))
		       (pushnew (list 433 41nbor sxid 32nbor) subcmplx :test #'set-equal?)))))))
	    ((= 4 (4sx-type sx))
	     (let ((32nbors (neighbors-of-type sx 3)))
	       (unless (< (length 32nbors) 2)
		 (do-tuples/c (id1 id2) 32nbors
			      (let ((sx1 nil) (sx2 nil))
				(when (and (setf sx1 (get-4simplex id1)) (setf sx2 (get-4simplex id2))
					   (4simplices-connected? id1 id2))
				  (pushnew (list 433 sxid id1 id2) subcmplx :test #'set-equal?)))))))))
    subcmplx))

(defun 3->3-v1-internal-122 (14id 23id1 23id2)
  (let ((14sx nil) (23sx1 nil) (23sx2 nil))
    (when (and (setf 14sx (get-4simplex 14id)) (setf 23sx1 (get-4simplex 23id1))
	       (setf 23sx2 (get-4simplex 23id2)))
      (let* ((oldint3sxids (list (link-id 14id 23id1) (link-id 14id 23id2) (link-id 23id1 23id2)))
	     (pts12 (4sx-lopts 23sx1))
	     (pt2 (set-difference pts12 (4sx-lopts 14sx)))
	     (pts34 (intersection (4sx-hipts 23sx1) (4sx-hipts 23sx2)))
	     (pts56 (set-exclusive-or (4sx-hipts 23sx1) (4sx-hipts 23sx2)))
	     (tmlo (4sx-tmlo 23sx1)) (tmhi (4sx-tmhi 23sx2))
	     (nbors (set-difference (unions (4sx-sx4ids 14sx)(4sx-sx4ids 23sx1) (4sx-sx4ids 23sx2))
				    (list 14id 23id1 23id2 0)))
	     (newint3sx1 (list 1 tmlo tmhi (append pt2 (firstl pts34) pts56)))
	     (newint3sx2 (list 1 tmlo tmhi (append pt2 (secondl pts34) pts56)))
	     (newint3sx3 (list 2 tmlo tmhi (append pts12 pts56)))
	     (newsxdata nil))
	(unless (or (gethash newint3sx1 *3SIMPLEX->ID*)
		    (gethash newint3sx2 *3SIMPLEX->ID*)
		    (gethash newint3sx3 *3SIMPLEX->ID*))
	  (setf newsxdata (list (list 1 tmlo tmhi (append pt2 (4sx-hipts 14sx)))
				(list 2 tmlo tmhi (append pts12 (3sx-hipts newint3sx1)))
				(list 2 tmlo tmhi (append pts12 (3sx-hipts newint3sx2)))))
	  (return-from 3->3-v1-internal-122 (list newsxdata nbors (list 14id 23id1 23id2) oldint3sxids DF33)))))))

(defun 3->3-v1-internal-433 (41id 32id1 32id2)
  (let ((41sx nil) (32sx1 nil) (32sx2 nil))
    (when (and (setf 41sx (get-4simplex 41id)) (setf 32sx1 (get-4simplex 32id1))
	       (setf 32sx2 (get-4simplex 32id2)))
      (let* ((oldint3sxids (list (link-id 41id 32id1) (link-id 41id 32id2) (link-id 32id1 32id2)))
	     (pts12 (4sx-hipts 32sx1))
	     (pt2 (set-difference pts12 (4sx-hipts 41sx)))
	     (pts34 (intersection (4sx-lopts 32sx1) (4sx-lopts 32sx2)))
	     (pts56 (set-exclusive-or (4sx-lopts 32sx1) (4sx-lopts 32sx2)))
	     (tmlo (4sx-tmlo 32sx2)) (tmhi (4sx-tmhi 32sx1))
	     (nbors (set-difference (unions (4sx-sx4ids 41sx)(4sx-sx4ids 32sx1) (4sx-sx4ids 32sx2))
				    (list 41id 32id1 32id2 0)))
	     (newint3sx1 (list 3 tmlo tmhi (append pts56 (firstl pts34) pt2)))
	     (newint3sx2 (list 3 tmlo tmhi (append pts56 (secondl pts34) pt2)))
	     (newint3sx3 (list 2 tmlo tmhi (append pts56 pts12)))
	     (newsxdata nil))
	(unless (or (gethash newint3sx1 *3SIMPLEX->ID*)
		    (gethash newint3sx2 *3SIMPLEX->ID*)
		    (gethash newint3sx3 *3SIMPLEX->ID*))
	  (setf newsxdata (list (list 4 tmlo tmhi (append (4sx-lopts 41sx) pt2))
				(list 3 tmlo tmhi (append (3sx-lopts newint3sx1) pts12))
				(list 3 tmlo tmhi (append (3sx-lopts newint3sx2) pts12))))
	  (return-from 3->3-v1-internal-433 (list newsxdata nbors (list 41id 32id1 32id2) oldint3sxids DF33)))))))

(defun try-3->3-v1 (sxid)
  (let ((subcmplx (3->3-v1-subcomplex sxid))
	(movedata nil))
    (unless (null subcmplx)
      (dolist (curr subcmplx)
	(cond ((= 122 (first curr))
	       (setf movedata (3->3-v1-internal-122 (second curr) (third curr) (fourth curr)))
	       (when movedata
		 (return-from try-3->3-v1 movedata)))
	      ((= 433 (first curr))
	       (setf movedata (3->3-v1-internal-433 (second curr) (third curr) (fourth curr)))
	       (when movedata
		 (return-from try-3->3-v1 movedata))))))))
;;;------------------------------------------------------------------------------------------------------[5]
;;; 3->3-v2
;;;
;;;     4 5 6
;;; t+1 -----------
;;;     1 2 3
;;; t   -----------
;;;
;;; INPUT: either a (2,3) or a (3,2) simplex
;;;
;;; 2(2,3) + 1(3,2) -> 1(2,3) + 2(3,2)
;;; 1(2,3) + 2(3,2) -> 2(2,3) + 1(3,2)
;;;
;;; (1 2|4 5 6) + (1 3|4 5 6)
;;; +
;;; (1 2 3|4 5)
;;; ->
;;; (2 3|4 5 6)
;;; +
;;; (1 2 3|4 6) + (1 2 3|5 6)
;;;------------------------------------------------------------------------------------------------------[5]
(defun 3->3-v2-subcomplex (sxid)
  "returns ((223or332 23or32id1 23or32id2 32or23id) ...)"
  (let ((sx '())
	(subcmplx '()))
    (when (setf sx (get-4simplex sxid))
      (let ((23nbors (neighbors-of-type sx 2))
	    (32nbors (neighbors-of-type sx 3)))
	(cond ((= 2 (4sx-type sx))
	       (dolist (23nbor 23nbors)
		 (dolist (32nbor 32nbors)
		   (when (4simplices-connected? 23nbor 32nbor)
		     (pushnew (list 223 sxid 23nbor 32nbor) subcmplx :test #'set-equal?))))
	       (when (>= (length 32nbors) 2)
		 (do-tuples/c (id1 id2) 32nbors
			      (when (4simplices-connected? id1 id2)
				(pushnew (list 332 id1 id2 sxid) subcmplx :test #'set-equal?)))))
	      ((= 3 (4sx-type sx))
	       (dolist (32nbor 32nbors)
		 (dolist (23nbor 23nbors)
		   (when (4simplices-connected? 32nbor 23nbor)
		     (pushnew (list 332 sxid 32nbor 23nbor) subcmplx :test #'set-equal?))))
	       (when (>= (length 23nbors) 2)
		 (do-tuples/c (id1 id2) 23nbors
			      (when (4simplices-connected? id1 id2)
				(pushnew (list 223 id1 id2 sxid) subcmplx :test #'set-equal?))))))))
    subcmplx))

(defun 3->3-v2-internal-223 (23id1 23id2 32id)
  (let ((23sx1 nil) (23sx2 nil) (32sx nil))
    (when (and (setf 23sx1 (get-4simplex 23id1)) (setf 23sx2 (get-4simplex 23id2))
	       (setf 32sx (get-4simplex 32id)))
      (let* ((oldint3sxids (list (link-id 23id1 23id2) (link-id 23id1 32id) (link-id 23id2 32id)))
	     (pts123 (4sx-lopts 32sx))
	     (pts23 (set-exclusive-or (4sx-lopts 23sx1) (4sx-lopts 23sx2)))
	     (pts456 (4sx-hipts 23sx1))
	     (pts45 (4sx-hipts 32sx))
	     (pt6 (set-difference pts456 pts45))
	     (tmlo (4sx-tmlo 23sx1)) (tmhi (4sx-tmhi 23sx2))
	     (nbors (set-difference (unions (4sx-sx4ids 23sx1)(4sx-sx4ids 23sx2) (4sx-sx4ids 32sx))
				    (list 23id1 23id2 32id 0)))
	     (newint3sx1 (list 2 tmlo tmhi (append pts23 (firstl pts45) pt6)))
	     (newint3sx2 (list 2 tmlo tmhi (append pts23 (secondl pts45) pt6)))
	     (newint3sx3 (list 3 tmlo tmhi (append pts123 pt6)))
	     (newsxdata nil))
	(unless (or (gethash newint3sx1 *3SIMPLEX->ID*)
		    (gethash newint3sx2 *3SIMPLEX->ID*)
		    (gethash newint3sx3 *3SIMPLEX->ID*))
	  (setf newsxdata (list (list 2 tmlo tmhi (append pts23 pts456))
				(list 3 tmlo tmhi (append pts123 (3sx-hipts newint3sx1)))
				(list 3 tmlo tmhi (append pts123 (3sx-hipts newint3sx2)))))
	  (return-from 3->3-v2-internal-223 (list newsxdata nbors (list 23id1 23id2 32id) oldint3sxids DF33)))))))

(defun 3->3-v2-internal-332 (32id1 32id2 23id)
  (let ((32sx1 nil) (32sx2 nil) (23sx nil))
    (when (and (setf 32sx1 (get-4simplex 32id1)) (setf 32sx2 (get-4simplex 32id2))
	       (setf 23sx (get-4simplex 23id)))
      (let* ((oldint3sxids (list (link-id 32id1 32id2) (link-id 32id1 23id) (link-id 32id2 23id)))
	     (pts123 (4sx-hipts 23sx))
	     (pts23 (set-exclusive-or (4sx-hipts 32sx1) (4sx-hipts 32sx2)))
	     (pts456 (4sx-lopts 32sx1))
	     (pts45 (4sx-lopts 23sx))
	     (pt6 (set-difference pts456 pts45))
	     (tmlo (4sx-tmlo 32sx1)) (tmhi (4sx-tmhi 32sx2))
	     (nbors (set-difference (unions (4sx-sx4ids 32sx1)(4sx-sx4ids 32sx2) (4sx-sx4ids 23sx))
				    (list 32id1 32id2 23id 0)))
	     (newint3sx1 (list 2 tmlo tmhi (append (firstl pts45) pt6 pts23)))
	     (newint3sx2 (list 2 tmlo tmhi (append (secondl pts45) pt6 pts23)))
	     (newint3sx3 (list 1 tmlo tmhi (append pt6 pts123)))
	     (newsxdata nil))
	(unless (or (gethash newint3sx1 *3SIMPLEX->ID*)
		    (gethash newint3sx2 *3SIMPLEX->ID*)
		    (gethash newint3sx3 *3SIMPLEX->ID*))
	  (setf newsxdata (list (list 3 tmlo tmhi (append pts456 pts23))
				(list 2 tmlo tmhi (append (3sx-lopts newint3sx1) pts123))
				(list 2 tmlo tmhi (append (3sx-lopts newint3sx2) pts123))))
	  (return-from 3->3-v2-internal-332 (list newsxdata nbors (list 32id1 32id2 23id) oldint3sxids DF33)))))))

(defun try-3->3-v2 (sxid)
  (let ((subcmplx (3->3-v2-subcomplex sxid))
	(movedata nil))
    (unless (null subcmplx)
      (dolist (curr subcmplx)
	(cond ((= 223 (first curr))
	       (setf movedata (3->3-v2-internal-223 (second curr) (third curr) (fourth curr)))
	       (when movedata
		 (return-from try-3->3-v2 movedata)))
	      ((= 332 (first curr))
	       (setf movedata (3->3-v2-internal-332 (second curr) (third curr) (fourth curr)))
	       (when movedata
		 (return-from try-3->3-v2 movedata))))))))
