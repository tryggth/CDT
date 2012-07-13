;; try-a->b methods returns the following list, IFF the move can be
;; successfully made (new3sxids nbors old3sxids oldTL2sxs oldSL2sxs
;; oldTL1sxs oldSL1sx fvector)

(defun 2plus1move (sxdata)
  (let ((new3sxids (make-3simplices-in-bulk (first sxdata))))
    (connect-3simplices-within-list new3sxids)
    (connect-3simplices-across-lists new3sxids (second sxdata))
    (remove-3simplices (third sxdata))
    (remove-tl2simplices (fourth sxdata))
    (remove-sl2simplices (fifth sxdata))
    (remove-tl1simplices (sixth sxdata))
    (remove-sl1simplices (seventh sxdata))
    (update-f-vector (eighth sxdata))
    (update-b-vector (ninth sxdata))))

;;------------------------------------------------------------------------[0]
;; 5
;;-------------- t+1
;; 2 3 4
;;-------------- t
;; 1
;;-------------- t-1
;;
;; input : either a (1,3) or a (3,1) 3-simplex
;;
;; 1(1,3) + 1(3,1) -> 3(1,3) + 3(3,1)
;;
;; (1|2 3 4) + (2 3 4|5) 
;; ->
;; (1|2 3 6) + (2 3 6|5) +
;; (1|2 6 4) + (2 6 4|5) +
;; (1|6 3 4) + (6 3 4|5)
;;------------------------------------------------------------------------[0]
(defun 2->6-subcomplex (sxid)
  (let ((subcmplx nil)
	(sx nil))
    (when (setf sx (get-3simplex sxid))
      (cond ((= 1 (3sx-type sx))
	     (when (/= 0 (nth 0 (3sx-sx3ids sx)))
	       (push (list sxid (nth 0 (3sx-sx3ids sx))) subcmplx)))
	    ((= 3 (3sx-type sx))
	     (when (/= 0 (nth 3 (3sx-sx3ids sx)))
	       (push (list (nth 3 (3sx-sx3ids sx)) sxid) subcmplx)))))
    subcmplx))
				  
(defun try-2->6 (sxid)
;;  (setf CURRENT-MOVE-IDENTIFIER "try-2->6")
  (dolist (curr (2->6-subcomplex sxid))
    (let* ((sx13 (get-3simplex (first curr)))
	   (sx31 (get-3simplex (second curr)))
	   (nbors (set-difference 
		   (union (3sx-sx3ids sx13) (3sx-sx3ids sx31)) 
		   curr))
	   (lopt (nth-point sx13 0))
	   (hipt (nth-point sx31 3))
	   (bigtr (3sx-lopts sx31))
	   (13tmlo (3sx-tmlo sx13)) 
	   (13tmhi (3sx-tmhi sx13))
	   (31tmlo (3sx-tmlo sx31)) 
	   (31tmhi (3sx-tmhi sx31))
	   (newpt (next-pt))
	   (oldSL2sxs `((,31tmlo ,bigtr)))
	   (newsxdata 
	    `((1 ,13tmlo ,13tmhi (,lopt ,newpt ,@(circular-subseq bigtr 0 2)))
	      (1 ,13tmlo ,13tmhi (,lopt ,newpt ,@(circular-subseq bigtr 1 2)))
	      (1 ,13tmlo ,13tmhi (,lopt ,newpt ,@(circular-subseq bigtr 2 2)))
	      (3 ,31tmlo ,31tmhi (,@(circular-subseq bigtr 0 2) ,newpt ,hipt))
	      (3 ,31tmlo ,31tmhi (,@(circular-subseq bigtr 1 2) ,newpt ,hipt))
	      (3 ,31tmlo ,31tmhi (,@(circular-subseq bigtr 2 2) 
				    ,newpt ,hipt)))))
      (return-from try-2->6 
	(list newsxdata nbors curr 
	      nil oldSL2sxs
	      nil nil
	      DF26 DB26)))))
;;------------------------------------------------------------------------[1]
;; 5
;;-------------- t+1
;; 2 3 4
;;-------------- t
;; 1
;;-------------- t-1
;;
;; input : either a (1,3) or a (3,1) 3-simplex
;;
;; 3(1,3) + 3(3,1) -> 1(1,3) + 1(3,1)
;;
;; (1|2 3 6) + (2 3 6|5) +
;; (1|2 6 4) + (2 6 4|5) +
;; (1|6 3 4) + (6 3 4|5)
;; ->
;; (1|2 3 4) + (2 3 4|5) 
;;------------------------------------------------------------------------[1]
(defun 6->2-subcomplex (sxid)
  "returns a list of the form ((13id1 13id2 13id3 31id3 31id2 31id1)...)"
  (let ((sx nil)
	(subcmplx nil))
    (when (setf sx (get-3simplex sxid))
      (cond ((= 1 (3sx-type sx))
	     (let ((31id (nth 0 (3sx-sx3ids sx)))
		   (31sx nil))
	       (when (setf 31sx (get-3simplex 31id))
		 (let ((13nbors (neighbors-of-type sx 1)))
		   (unless (< (length 13nbors) 2)
		     (do-tuples/c (currid nextid) 13nbors
		       (let ((curr nil) (next nil))
			 (when (and (setf curr (get-3simplex currid)) 
				    (setf next (get-3simplex nextid))
				    (3simplices-connected? currid nextid)
				    (3simplices-connected? 
				     (nth 0 (3sx-sx3ids curr))
				     (nth 0 (3sx-sx3ids next)))
				    (3simplices-connected? 
				     (nth 0 (3sx-sx3ids curr)) 31id)
				    (3simplices-connected? 
				     (nth 0 (3sx-sx3ids next)) 31id))
			   (pushnew (list sxid currid nextid 
					  (nth 0 (3sx-sx3ids next)) 
					  (nth 0 (3sx-sx3ids curr))
					  31id)
				    subcmplx :test #'set-equal?)))))))))
	    ((= 3 (3sx-type sx))
	     (let ((13id (nth 3 (3sx-sx3ids sx)))
		   (13sx nil))
	       (when (setf 13sx (get-3simplex 13id))
		 (let ((31nbors (neighbors-of-type sx 3)))
		   (unless (< (length 31nbors) 2)
		     (do-tuples/c (currid nextid) 31nbors
		       (let ((curr nil) (next nil))
			 (when (and (setf curr (get-3simplex currid)) 
				    (setf next (get-3simplex nextid))
				    (3simplices-connected? currid nextid)
				    (3simplices-connected? 
				     (nth 3 (3sx-sx3ids curr))
				     (nth 3 (3sx-sx3ids next)))
				    (3simplices-connected? 
				     (nth 3 (3sx-sx3ids curr)) 13id)
				    (3simplices-connected? 
				     (nth 3 (3sx-sx3ids next)) 13id))
			   (pushnew (list 13id 
					  (nth 3 (3sx-sx3ids next)) 
					  (nth 3 (3sx-sx3ids curr))
					  currid nextid sxid) 
				    subcmplx :test #'set-equal?)))))))))))
    subcmplx))
	    
(defun try-6->2 (sxid)
  ;;  (setf CURRENT-MOVE-IDENTIFIER "try-8->2")
  (dolist (subcx (6->2-subcomplex sxid))
    (let* ((1id1 (first subcx)) 
	   (1id2 (second subcx)) 
	   (1id3 (third subcx)) 
	   (3id3 (fourth subcx)) 
	   (3id2 (fifth subcx)) 
	   (3id1 (sixth subcx)) 
	   (1sx1 (get-3simplex 1id1)) 
	   (1sx2 (get-3simplex 1id2)) 
	   (1sx3 (get-3simplex 1id3)) 
	   (3sx1 (get-3simplex 3id1)) 
	   (3sx2 (get-3simplex 3id2))
	   (3sx3 (get-3simplex 3id3)) 
	   (13tmlo (3sx-tmlo 1sx1)) 
	   (13tmhi (3sx-tmhi 1sx1))
	   (31tmlo (3sx-tmlo 3sx1)) 
	   (31tmhi (3sx-tmhi 3sx1))
	   (nbors (set-difference 
		   (unions (3sx-sx3ids 1sx1) (3sx-sx3ids 1sx2) 
			   (3sx-sx3ids 1sx3) (3sx-sx3ids 3sx1) 
			   (3sx-sx3ids 3sx2) (3sx-sx3ids 3sx3))
		   (list 0 1id1 1id2 1id3 3id1 3id2 3id3)))
	   (pt1 (3sx-lopts 1sx1)) ; (1)
	   (pt5 (3sx-hipts 3sx1)) ; (5)
	   (pt6 (intersections (3sx-hipts 1sx1) (3sx-hipts 1sx2)
			       (3sx-hipts 1sx3))) ; (6)
	   (pts234 (set-difference (unions (3sx-hipts 1sx1) (3sx-hipts 1sx2)
					   (3sx-hipts 1sx3))
				   pt6)) ; (2 3 4)
	   (oldSL2sxs `((,13tmhi (,@pt6 ,@(circular-subseq pts234 0 2)))
			(,13tmhi (,@pt6 ,@(circular-subseq pts234 1 2)))
			(,13tmhi (,@pt6 ,@(circular-subseq pts234 2 2)))))
	   (oldTL2sxs `((1 ,13tmlo (,@pt1 ,(first pts234) ,@pt6))
			(1 ,13tmlo (,@pt1 ,(second pts234) ,@pt6))
			(1 ,13tmlo (,@pt1 ,(third pts234) ,@pt6))
			(2 ,31tmlo (,(first pts234) ,@pt6 ,@pt5))
			(2 ,31tmlo (,(second pts234) ,@pt6 ,@pt5))
			(2 ,31tmlo (,(third pts234) ,@pt6 ,@pt5))))
	   (oldSL1sxs `((,13tmhi (,(first pts234) ,@pt6))
			(,13tmhi (,(second pts234) ,@pt6))
			(,13tmhi (,(third pts234) ,@pt6))))
	   (oldTL1sxs `((1 ,13tmlo (,@pt1 ,@pt6))
			(1 ,31tmlo (,@pt6 ,@pt5))))
	   (newsxdata nil))
      (unless (gethash `(,13tmhi ,pts234) *SL2SIMPLEX->ID*)
	(setf newsxdata 
	      `((1 ,13tmlo ,13tmhi (,@pt1 ,@pts234))
		(3 ,31tmlo ,31tmhi (,@pts234 ,@pt5))))
	(return-from try-6->2 
	  (list newsxdata nbors subcx
		oldTL2sxs oldSL2sxs
		oldTL1sxs oldSL1sxs
		DF62 DB62))))))

;;------------------------------------------------------------------------[2]
;; 6
;;-------------- t+1
;; 2 3 4 5
;;-------------- t
;; 1
;;-------------- t-1
;;
;; input : either a (1,3) or a (3,1) 3-simplex
;;
;; 2(1,3) + 2(3,1) -> 2(1,3) + 2(3,1)
;;
;; (1|2 3 4) + (1|3 4 5) +
;; (2 3 4|6) + (3 4 5|6)
;; ->
;; (1|2 4 5) + (1|2 3 5) +
;; (2 4 5|6) + (2 3 5|6)
;;------------------------------------------------------------------------[2]
(defun 4->4-subcomplex (sxid)
  "returns a list of the form ((13id1 13id2 31id2 31id1)...)"
  (let ((sx nil)
	(subcmplx nil))
    (when (setf sx (get-3simplex sxid))
      (cond ((= 1 (3sx-type sx))
	     (let ((31id (nth 0 (3sx-sx3ids sx))) (31sx nil))
	       (when (setf 31sx (get-3simplex 31id))
		 (let ((13nbors (neighbors-of-type sx 1)))
		   (dolist (13nbor 13nbors)
		     (when (3simplices-connected? 
			    (nth 0 (3sx-sx3ids (get-3simplex 13nbor))) 31id)
		       (pushnew (list sxid 13nbor 
				      (nth 0 (3sx-sx3ids 
					      (get-3simplex 13nbor)))
				      31id)
				subcmplx :test #'set-equal?)))))))
	    ((= 3 (3sx-type sx))
	     (let ((13id (nth 3 (3sx-sx3ids sx))) (13sx nil))
	       (when (setf 13sx (get-3simplex 13id))
		 (let ((31nbors (neighbors-of-type sx 3)))
		   (dolist (31nbor 31nbors)
		     (when (3simplices-connected? 
			    (nth 3 (3sx-sx3ids (get-3simplex 31nbor))) 13id)
		       (pushnew (list 13id 
				      (nth 3 (3sx-sx3ids 
					      (get-3simplex 31nbor)))
				      31nbor sxid)
				subcmplx :test #'set-equal?)))))))))
    subcmplx))

(defun try-4->4 (sxid)
;;  (setf CURRENT-MOVE-IDENTIFIER "try-4->4")
  (dolist (subcx (4->4-subcomplex sxid))
    (let* ((1id1 (first subcx)) 
	   (1id2 (second subcx)) 
	   (3id2 (third subcx)) 
	   (3id1 (fourth subcx)) 
	   (1sx1 (get-3simplex 1id1)) 
	   (1sx2 (get-3simplex 1id2)) 
	   (3sx1 (get-3simplex 3id1)) 
	   (3sx2 (get-3simplex 3id2))
	   (13tmlo (3sx-tmlo 1sx1)) 
	   (13tmhi (3sx-tmhi 1sx1))
	   (31tmlo (3sx-tmlo 3sx1)) 
	   (31tmhi (3sx-tmhi 3sx1))
	   (nbors (set-difference 
		   (unions (3sx-sx3ids 1sx1) (3sx-sx3ids 1sx2) 
			   (3sx-sx3ids 3sx1) (3sx-sx3ids 3sx2))
		   (list 0 1id1 1id2 3id1 3id2)))
	   (pt1 (3sx-lopts 1sx1)) ; (1)
	   (pt6 (3sx-hipts 3sx1)) ; (6)
	   (pts25 (set-exclusive-or (3sx-hipts 1sx1) (3sx-hipts 1sx2)));(2 5)
	   (pts34 (intersection (3sx-lopts 3sx1) (3sx-lopts 3sx2)));(3 4)
	   (oldTL2sxs `((1 ,13tmlo (,@pt1 ,@pts34))
			(2 ,31tmlo (,@pts34 ,@pt6))))
	   (oldSL2sxs `((,13tmhi (,(first pts25) ,@pts34))
			(,31tmlo (,@pts34 ,(second pts25)))))
	   (oldSL1sxs `((,13tmhi ,pts34)))
	   (newsxdata nil))
      (unless (gethash `(,31tmlo ,pts25) *SL1SIMPLEX->ID*)
	(setf newsxdata 
	      `((1 ,13tmlo ,13tmhi (,@pt1 ,@pts25 ,(first pts34)))
		(1 ,13tmlo ,13tmhi (,@pt1 ,@pts25 ,(second pts34)))
		(3 ,31tmlo ,31tmhi (,@pts25 ,(first pts34) ,@pt6))
		(3 ,31tmlo ,31tmhi (,@pts25 ,(second pts34) ,@pt6))))
	(return-from try-4->4 
	  (list newsxdata nbors subcx 
		oldTL2sxs oldSL2sxs
		nil oldSL1sxs
		DF44 DB44))))))
;;------------------------------------------------------------------------[3]
;; 2 3 4
;;-------------- t+1
;; 1 5
;;-------------- t
;;
;; input : either a (1,3) or a (3,1) or a (2,2) 3-simplex
;;
;; 1(1,3) + 1(2,2) -> 1(1,3) + 2(2,2)
;;
;; (1|2 3 4) + (1 5|3 4) -> (5|2 3 4) + (1 5|2 3) + (1 5|2 4)
;; or
;; (2 3 4|1) + (3 4|1 5) -> (2 3 4|5) + (2 3|1 5) + (2 4|1 5)
;;------------------------------------------------------------------------[3]
(defun 2->3-subcomplex (sxid)
  "returns a list of the form ((1or3 13or31id 22id)...)
   where the first number 1 or 3 tells us about the
   type of the simplex participating in the move"
  (let ((sx nil)
	(subcmplx nil))
    (when (setf sx (get-3simplex sxid))
      (cond ((or (= 1 (3sx-type sx)) (= 3 (3sx-type sx)))
	     (let ((22nbors (neighbors-of-type sx 2)))
	       (dolist (22nbor 22nbors)
		 (pushnew (list (3sx-type sx) sxid 22nbor) 
			  subcmplx :test #'set-equal?))))
	    ((= 2 (3sx-type sx))
	     (let ((13nbors (append (neighbors-of-type sx 1))))
	       (dolist (13nbor 13nbors)
		 (pushnew (list 1 13nbor sxid) subcmplx :test #'set-equal?)))
	     (let ((31nbors (append (neighbors-of-type sx 3))))
	       (dolist (31nbor 31nbors)
		 (pushnew (list 3 31nbor sxid) 
			  subcmplx :test #'set-equal?))))))
    subcmplx))

(defun 2->3-move-internal-12 (13id 22id)
;;  (setf CURRENT-MOVE-IDENTIFIER "2->3-move-internal-12")
  (let ((13sx nil) (22sx nil))
    (when (and (setf 13sx (get-3simplex 13id)) 
	       (setf 22sx (get-3simplex 22id)))
      (let* ((pts234 (3sx-hipts 13sx));(2 3 4)
	     (pts34 (3sx-hipts 22sx));(3 4)
	     (pts15 (3sx-lopts 22sx));(1 5)
	     (pt1 (3sx-lopts 13sx));(1)
	     (pt5 (set-difference pts15 pt1));(5)
	     (pt2 (set-difference pts234 pts34));(2)
	     (nbors (set-difference 
		     (union (3sx-sx3ids 13sx) (3sx-sx3ids 22sx)) 
		     (list 0 13id 22id)))
	     (tmlo (3sx-tmlo 22sx))
	     (tmhi (3sx-tmhi 22sx))
	     (oldTL2sxs `((1 ,tmlo (,@pt1 ,@pts34))))
	     (newsxdata nil))
	(unless (gethash `(1 ,tmlo (,@pt5 ,@pt2)) *TL1SIMPLEX->ID*)
	  (setf newsxdata 
		`((1 ,tmlo ,tmhi (,@pt5 ,@pts234))
		  (2 ,tmlo ,tmhi (,@pt1 ,@pt5 ,@pt2 ,(first pts34)))
		  (2 ,tmlo ,tmhi (,@pt1 ,@pt5 ,@pt2 ,(second pts34)))))
	  (return-from 2->3-move-internal-12 
	    (list newsxdata nbors (list 13id 22id) 
		  oldTL2sxs nil
		  nil nil
		  DF23 (DB23 13id))))))))

(defun 2->3-move-internal-32 (31id 22id)
;;  (setf CURRENT-MOVE-IDENTIFIER "2->3-move-internal-32")
  (let ((31sx nil) (22sx nil))
    (when (and (setf 31sx (get-3simplex 31id)) 
	       (setf 22sx (get-3simplex 22id)))
      (let* ((pts234 (3sx-lopts 31sx));(2 3 4)
	     (pts34 (3sx-lopts 22sx));(3 4)
	     (pts15 (3sx-hipts 22sx));(1 5)
	     (pt1 (3sx-hipts 31sx));(1)
	     (pt5 (set-difference pts15 pt1));(5)
	     (pt2 (set-difference pts234 pts34));(2)
	     (nbors (set-difference 
		     (union (3sx-sx3ids 31sx) (3sx-sx3ids 22sx)) 
		     (list 0 31id 22id)))
	     (tmlo (3sx-tmlo 22sx))
	     (tmhi (3sx-tmhi 22sx))
	     (oldTL2sxs `((2 ,tmlo (,@pts34 ,@pt1))))
	     (newsxdata nil))
	(unless (gethash `(1 ,tmlo (,@pt5 ,@pt2)) *TL1SIMPLEX->ID*)
	  (setf newsxdata 
		`((3 ,tmlo ,tmhi (,@pts234 ,@pt5))
		  (2 ,tmlo ,tmhi (,@pt2 ,(first pts34) ,@pt1 ,@pt5))
		  (2 ,tmlo ,tmhi (,@pt2 ,(second pts34) ,@pt1 ,@pt5))))
	  (return-from 2->3-move-internal-32 
	    (list newsxdata nbors (list 31id 22id) 
		  oldTL2sxs nil
		  nil nil
		  DF23 (DB23 31id))))))))

(defun try-2->3 (sxid)
  (let ((subcmplx (2->3-subcomplex sxid))
	(movedata nil))
    (unless (null subcmplx)
      (dolist (curr subcmplx)
	(cond ((= 1 (first curr))
	       (setf movedata (2->3-move-internal-12 
			       (second curr) (third curr)))
	       (when movedata
		 (return-from try-2->3 movedata)))
	      ((= 3 (first curr))
	       (setf movedata (2->3-move-internal-32 
			       (second curr) (third curr)))
	       (when movedata
		 (return-from try-2->3 movedata))))))))
;;------------------------------------------------------------------------[4]
;; 2 3 4
;;-------------- t+1
;; 1 5
;;-------------- t
;;
;; input : either a (1,3) or a (3,1) or a (2,2) 3-simplex
;;
;; 1(1,3) + 2(2,2) -> 1(1,3) + 1(2,2)
;;
;; (5|2 3 4) + (1 5|2 3) + (1 5|2 4) -> (1|2 3 4) + (1 5|3 4)
;; or
;; (2 3 4|5) + (2 3|1 5) + (2 4|1 5) -> (2 3 4|1) + (3 4|1 5)
;;------------------------------------------------------------------------[4]
(defun 3->2-subcomplex (sxid)
  "returns a list of the form ((1or3 13or31id 22id1 22id2)...) 
   where the first number 1 or 3 tells us 
   about the type of the simplex participating in the move"
  (let ((sx nil)
	(subcmplx nil))
    (when (setf sx (get-3simplex sxid))
      (cond ((or (= 1 (3sx-type sx)) (= 3 (3sx-type sx)))
	     (let ((22nbors (neighbors-of-type sx 2)))
	       (dolist (22nbor 22nbors)
		 (let ((22sx (get-3simplex 22nbor)))
		   (when 22sx
		     (let ((22nborsof22nbor (neighbors-of-type 22sx 2)))
		       (dolist (22nborof22nbor 22nborsof22nbor)
			 (when (3simplices-connected? 22nborof22nbor sxid)
			   (pushnew (list (3sx-type sx) sxid 22nbor 
					  22nborof22nbor) 
				    subcmplx :test #'set-equal?)))))))))
	    ((= 2 (3sx-type sx))
	     (let ((22nbors (neighbors-of-type sx 2)))
	       (dolist (22nbor 22nbors)
		 (let ((22sx (get-3simplex 22nbor)))
		   (when 22sx
		     (let ((13nborsof22nbor (neighbors-of-type 22sx 1)))
		       (dolist (13nborof22nbor 13nborsof22nbor)
			 (when (3simplices-connected? 13nborof22nbor sxid)
			   (pushnew (list 1 13nborof22nbor sxid 22nbor) 
				    subcmplx :test #'set-equal?))))
		     (let ((31nborsof22nbor (neighbors-of-type 22sx 3)))
		       (dolist (31nborof22nbor 31nborsof22nbor)
			 (when (3simplices-connected? 31nborof22nbor sxid)
			   (pushnew (list 3 31nborof22nbor sxid 22nbor) 
				    subcmplx :test #'set-equal?)))))))))))
    subcmplx))

(defun 3->2-move-internal-122 (13id 22id1 22id2)
  "the (3,2) move performed on a (1,3) simplex attached to two (2,2) 
   simplices"
  (let ((13sx nil) (22sx1 nil) (22sx2 nil))
    (when (and (setf 13sx (get-3simplex 13id)) 
	       (setf 22sx1 (get-3simplex 22id1))
	       (setf 22sx2 (get-3simplex 22id2)))
      (let* ((pts234 (3sx-hipts 13sx));(2 3 4)
	     (pts15 (3sx-lopts 22sx1));(1 5)
	     (pt5 (3sx-lopts 13sx));(5)
	     (pt1 (set-difference pts15 pt5));(1)
	     (pts34 (set-exclusive-or 
		     (3sx-hipts 22sx1) (3sx-hipts 22sx2)));(3 4)
	     (pt2 (set-difference pts234 pts34));(2)
	     (nbors (set-difference 
		     (unions (3sx-sx3ids 13sx) (3sx-sx3ids 22sx1) 
			     (3sx-sx3ids 22sx2)) 
		     (list 0 13id 22id1 22id2)))
	     (tmlo (3sx-tmlo 22sx1))
	     (tmhi (3sx-tmhi 22sx1))
	     (oldTL2sxs `((1 ,tmlo (,@pt5 ,@pt2 ,(first pts34)))
			  (1 ,tmlo (,@pt5 ,@pt2 ,(second pts34)))
			  (2 ,tmlo (,@pts15 ,@pt2))))
	     (oldTL1sxs `((1 ,tmlo (,@pt5 ,@pt2))))
	     (newsxdata nil))
	(unless (gethash `(1 ,tmlo (,@pt1 ,@pts34)) *TL2SIMPLEX->ID*)
	  (setf newsxdata 
		`((1 ,tmlo ,tmhi (,@pt1 ,@pts234))
		  (2 ,tmlo ,tmhi (,@pts15 ,@pts34))))
	  (return-from 3->2-move-internal-122 
	    (list newsxdata nbors (list 13id 22id1 22id2)
		  oldTL2sxs nil
		  oldTL1sxs nil
		  DF32 (DB32 13id))))))))

(defun 3->2-move-internal-322 (31id 22id1 22id2)
  (let ((31sx nil) (22sx1 nil) (22sx2 nil))
    (when (and (setf 31sx (get-3simplex 31id)) 
	       (setf 22sx1 (get-3simplex 22id1))
	       (setf 22sx2 (get-3simplex 22id2)))
      (let* ((pts234 (3sx-lopts 31sx));(2 3 4)
	     (pts15 (3sx-hipts 22sx1));(1 5)
	     (pt5 (3sx-hipts 31sx));(5)
	     (pt1 (set-difference pts15 pt5));(1)
	     (pts34 (set-exclusive-or 
		     (3sx-lopts 22sx1) (3sx-lopts 22sx2)));(3 4)
	     (pt2 (set-difference pts234 pts34))
	     (nbors (set-difference 
		     (unions (3sx-sx3ids 31sx) (3sx-sx3ids 22sx1) 
			     (3sx-sx3ids 22sx2)) 
		     (list 0 31id 22id1 22id2)))
	     (tmlo (3sx-tmlo 22sx1))
	     (tmhi (3sx-tmhi 22sx1))
	     (oldTL2sxs `((2 ,tmlo (,@pt2 ,(first pts34) ,@pt5))
			  (2 ,tmlo (,@pt2 ,(second pts34) ,@pt5))
			  (1 ,tmlo (,@pt2 ,@pts15))))
	     (oldTL1sxs `((1 ,tmlo (,@pt2 ,@pt5))))
	     (newsxdata nil))
	(unless (gethash `(2 ,tmlo (,@pts34 ,@pt1)) *TL2SIMPLEX->ID*)
	  (setf newsxdata 
		`((3 ,tmlo ,tmhi (,@pts234 ,@pt1))
		  (2 ,tmlo ,tmhi (,@pts34 ,@pts15))))
	  (return-from 3->2-move-internal-322 
	    (list newsxdata nbors (list 31id 22id1 22id2)
		  oldTL2sxs nil
		  oldTL1sxs nil
		  DF32 (DB32 31id))))))))

(defun try-3->2 (sxid)
  (let ((subcmplx (3->2-subcomplex sxid))
	(movedata nil))
    (unless (null subcmplx)
      (dolist (curr subcmplx)
	(cond ((= 1 (first curr))
	       (setf movedata (3->2-move-internal-122 
			       (second curr) (third curr) (fourth curr)))
	       (when movedata
		 (return-from try-3->2 movedata)))
	      ((= 3 (first curr))
	       (setf movedata (3->2-move-internal-322 
			       (second curr) (third curr) (fourth curr)))
	       (when movedata 
		 (return-from try-3->2 movedata))))))))
