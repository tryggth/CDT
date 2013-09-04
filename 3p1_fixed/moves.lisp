(declaim (optimize (speed 3)
		   (compilation-speed 0)
		   (debug 0)
		   (safety 0)))
;; try-a->b methods return 
;;
;; (new4sxdata nbors old4sxids oldtl3sxs oldsl3sxs oldtl2sxs oldsl2sxs oldtl1sxs oldsl1sxs fvector) 
;; IFF the move can be successfully made 


(defun 3plus1move (sxdata)
  (let ((new4sxids (make-4simplices-in-bulk (first sxdata))))
    (connect-4simplices-within-list new4sxids)
    (connect-4simplices-across-lists new4sxids (second sxdata))
    (remove-4simplices (third sxdata))
    (remove-tl3simplices (fourth sxdata))
    (remove-sl3simplices (fifth sxdata))
    (remove-tl2simplices (sixth sxdata))
    (remove-sl2simplices (seventh sxdata))
    (remove-tl1simplices (eighth sxdata))
    (remove-sl1simplices (ninth sxdata))
    (update-f-vector (tenth sxdata))))


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
;;  (setf CURRENT-MOVE-IDENTIFIER "2->8")
  (dolist (curr (2->8-subcomplex sxid))
    (let* ((14id (first curr))
	   (41id (second curr))
	   (14sx (get-4simplex 14id))
	   (41sx (get-4simplex 41id))
	   (14tmlo (4sx-tmlo 14sx)) 
	   (14tmhi (4sx-tmhi 14sx))
	   (41tmlo (4sx-tmlo 41sx)) 
	   (41tmhi (4sx-tmhi 41sx))
	   (nbors (set-difference 
		   (union (4sx-sx4ids 14sx) (4sx-sx4ids 41sx)) 
		   (list 0 14id 41id)))
	   (newpt (next-pt))
	   (pts2345 (4sx-hipts 14sx));;(2 3 4 5)
	   (pt1 (4sx-lopts 14sx));;(1)
	   (pt6 (4sx-hipts 41sx));;(6)
	   (oldslint3sxs `((,41tmlo ,pts2345)))
	   (newsxdata `((1 ,14tmlo ,14tmhi 
			   (,@pt1 ,newpt ,(first pts2345) ,(second pts2345)
				  ,(third pts2345)))
			(1 ,14tmlo ,14tmhi 
			   (,@pt1 ,newpt ,(second pts2345) ,(third pts2345)
				  ,(fourth pts2345)))
			(1 ,14tmlo ,14tmhi 
			   (,@pt1 ,newpt ,(third pts2345) ,(fourth pts2345)
				  ,(first pts2345)))
			(1 ,14tmlo ,14tmhi 
			   (,@pt1 ,newpt ,(fourth pts2345) ,(first pts2345)
				  ,(second pts2345)))
			(4 ,41tmlo ,41tmhi 
			   (,newpt ,(first pts2345) ,(second pts2345)
				   ,(third pts2345) ,@pt6))
			(4 ,41tmlo ,41tmhi 
			   (,newpt ,(second pts2345) ,(third pts2345)
				   ,(fourth pts2345) ,@pt6))
			(4 ,41tmlo ,41tmhi 
			   (,newpt ,(third pts2345) ,(fourth pts2345)
				   ,(first pts2345) ,@pt6))
			(4 ,41tmlo ,41tmhi 
			   (,newpt ,(fourth pts2345) ,(first pts2345)
				   ,(second pts2345) ,@pt6)))))
      (unless (or (member 14id *FIXED-4SXID*) (member 41id *FIXED-4SXID*))
	(return-from try-2->8 
	  (list newsxdata nbors curr nil oldslint3sxs nil nil nil nil DF28))))))
;;;--------------------------------------------------------------------------[9]
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
;;;--------------------------------------------------------------------------[9]
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
			 (when (and (setf sx1 (get-4simplex id1))
				    (setf sx2 (get-4simplex id2))
				    (setf sx3 (get-4simplex id3))
				    (4simplices-connected? id1 id2)
				    (4simplices-connected? id2 id3)
				    (4simplices-connected? id3 id1)
				    (4simplices-connected? 
				     (nth 0 (4sx-sx4ids sx1))
				     (nth 0 (4sx-sx4ids sx2)))
				    (4simplices-connected? 
				     (nth 0 (4sx-sx4ids sx2))
				     (nth 0 (4sx-sx4ids sx3)))
				    (4simplices-connected? 
				     (nth 0 (4sx-sx4ids sx3))
				     (nth 0 (4sx-sx4ids sx1)))
				    (4simplices-connected? 
				     (nth 0 (4sx-sx4ids sx1)) 41id)
				    (4simplices-connected? 
				     (nth 0 (4sx-sx4ids sx2)) 41id)
				    (4simplices-connected? 
				     (nth 0 (4sx-sx4ids sx3)) 41id))
			   (pushnew (list sxid id1 id2 id3 
					  (nth 0 (4sx-sx4ids sx3)) 
					  (nth 0 (4sx-sx4ids sx2))
					  (nth 0 (4sx-sx4ids sx1)) 41id)
				    subcmplx :test #'set-equal?)))))))))
	    ((= 4 (4sx-type sx))
	     (let ((14id (nth 4 (4sx-sx4ids sx)))(14sx nil))
	       (when (setf 14sx (get-4simplex 14id))
		 (let ((41nbors (neighbors-of-type sx 4)))
		   (unless (< (length 41nbors) 3)
		     (do-tuples/c (id1 id2 id3) 41nbors
		       (let ((sx1 nil) (sx2 nil) (sx3 nil))
			 (when (and (setf sx1 (get-4simplex id1)) 
				    (setf sx2 (get-4simplex id2))
				    (setf sx3 (get-4simplex id3)) 
				    (4simplices-connected? id1 id2)
				    (4simplices-connected? id2 id3) 
				    (4simplices-connected? id3 id1)
				    (4simplices-connected? 
				     (nth 4 (4sx-sx4ids sx1))
				     (nth 4 (4sx-sx4ids sx2)))
				    (4simplices-connected? 
				     (nth 4 (4sx-sx4ids sx2))
				     (nth 4 (4sx-sx4ids sx3)))
				    (4simplices-connected? 
				     (nth 4 (4sx-sx4ids sx3))
				     (nth 4 (4sx-sx4ids sx1)))
				    (4simplices-connected? 
				     (nth 4 (4sx-sx4ids sx1)) 14id)
				    (4simplices-connected? 
				     (nth 4 (4sx-sx4ids sx2)) 14id)
				    (4simplices-connected? 
				     (nth 4 (4sx-sx4ids sx3)) 14id))
			   (pushnew (list 14id (nth 4 (4sx-sx4ids sx3)) 
					  (nth 4 (4sx-sx4ids sx2))
					  (nth 4 (4sx-sx4ids sx1)) id1 id2 id3 
					  sxid)
				    subcmplx :test #'set-equal?)))))))))))
    subcmplx))


(defun try-8->2 (sxid)
;;  (setf CURRENT-MOVE-IDENTIFIER "8->2")
  (dolist (subcx (8->2-subcomplex sxid))
    (let* ((1id1 (first subcx)) 
	   (1id2 (second subcx)) 
	   (1id3 (third subcx)) 
	   (1id4 (fourth subcx)) 
	   (4id4 (fifth subcx)) 
	   (4id3 (sixth subcx)) 
	   (4id2 (seventh subcx)) 
	   (4id1 (eighth subcx))
	   (1sx1 (get-4simplex 1id1)) 
	   (1sx2 (get-4simplex 1id2)) 
	   (1sx3 (get-4simplex 1id3)) 
	   (1sx4 (get-4simplex 1id4))
	   (4sx1 (get-4simplex 4id1)) 
	   (4sx2 (get-4simplex 4id2))
	   (4sx3 (get-4simplex 4id3)) 
	   (4sx4 (get-4simplex 4id4))
	   (14tmlo (4sx-tmlo 1sx1)) 
	   (14tmhi (4sx-tmhi 1sx1))
	   (41tmlo (4sx-tmlo 4sx1)) 
	   (41tmhi (4sx-tmhi 4sx1))
	   (nbors (set-difference 
		   (unions (4sx-sx4ids 1sx1) (4sx-sx4ids 1sx2) 
			   (4sx-sx4ids 1sx3) (4sx-sx4ids 1sx4) 
			   (4sx-sx4ids 4sx1) (4sx-sx4ids 4sx2) 
			   (4sx-sx4ids 4sx3) (4sx-sx4ids 4sx4))
		   (list 0 1id1 1id2 1id3 1id4 4id1 4id2 4id3 4id4)))
	   (pt1 (4sx-lopts 1sx1));;(1)
	   (pt6 (4sx-hipts 4sx1));;(6)
	   (pt7 (intersections (4sx-hipts 1sx1) (4sx-hipts 1sx2)
			       (4sx-hipts 1sx3) (4sx-hipts 1sx4)));;(7)
	   (pts2345 (set-difference
		     (unions (4sx-hipts 1sx1) (4sx-hipts 1sx2)
			     (4sx-hipts 1sx3) (4sx-hipts 1sx4))
		     pt7));;(2 3 4 5)
	   (oldslint3sxs `((,14tmhi (,@pt7 ,@(circular-subseq pts2345 0 3)))
			   (,14tmhi (,@pt7 ,@(circular-subseq pts2345 1 3)))
			   (,14tmhi (,@pt7 ,@(circular-subseq pts2345 2 3)))
			   (,14tmhi (,@pt7 ,@(circular-subseq pts2345 3 3)))))
	   (oldtlint3sxs `((1 ,14tmlo 
			      (,@pt1 ,@pt7 ,(first pts2345) ,(second pts2345)))
			   (1 ,14tmlo 
			      (,@pt1 ,@pt7 ,(first pts2345) ,(third pts2345)))
			   (1 ,14tmlo 
			      (,@pt1 ,@pt7 ,(first pts2345) ,(fourth pts2345)))
			   (1 ,14tmlo 
			      (,@pt1 ,@pt7 ,(second pts2345) ,(third pts2345)))
			   (1 ,14tmlo 
			      (,@pt1 ,@pt7 ,(second pts2345) ,(fourth pts2345)))
			   (1 ,14tmlo 
			      (,@pt1 ,@pt7 ,(third pts2345) ,(fourth pts2345)))
			   (3 ,41tmlo 
			      (,@pt7 ,(first pts2345) ,(second pts2345) ,@pt6))
			   (3 ,41tmlo 
			      (,@pt7 ,(first pts2345) ,(third pts2345) ,@pt6))
			   (3 ,41tmlo 
			      (,@pt7 ,(first pts2345) ,(fourth pts2345) ,@pt6))
			   (3 ,41tmlo 
			      (,@pt7 ,(second pts2345) ,(third pts2345) ,@pt6))
			   (3 ,41tmlo 
			      (,@pt7 ,(second pts2345) ,(fourth pts2345) ,@pt6))
			   (3 ,41tmlo 
			      (,@pt7 ,(third pts2345) ,(fourth pts2345) ,@pt6))))
	   (oldslint2sxs `((,14tmhi (,@pt7 ,(first pts2345) ,(second pts2345)))
			   (,14tmhi (,@pt7 ,(first pts2345) ,(third pts2345)))
			   (,14tmhi (,@pt7 ,(first pts2345) ,(fourth pts2345)))
			   (,14tmhi (,@pt7 ,(second pts2345) ,(third pts2345)))
			   (,14tmhi (,@pt7 ,(second pts2345) ,(fourth pts2345)))
			   (,14tmhi (,@pt7 ,(third pts2345) ,(fourth pts2345)))))
	   (oldtlint2sxs `((1 ,14tmlo (,@pt1 ,(first pts2345) ,@pt7))
			   (1 ,14tmlo (,@pt1 ,(second pts2345) ,@pt7))
			   (1 ,14tmlo (,@pt1 ,(third pts2345) ,@pt7))
			   (1 ,14tmlo (,@pt1 ,(fourth pts2345) ,@pt7))
			   (2 ,41tmlo (,(first pts2345) ,@pt7 ,@pt6))
			   (2 ,41tmlo (,(second pts2345) ,@pt7 ,@pt6))
			   (2 ,41tmlo (,(third pts2345) ,@pt7 ,@pt6))
			   (2 ,41tmlo (,(fourth pts2345) ,@pt7 ,@pt6))))
	   (oldtlint1sxs `((1 ,14tmlo (,@pt1 ,@pt7))
			   (1 ,41tmlo (,@pt7 ,@pt6))))
	   (oldslint1sxs `((,14tmhi (,(first pts2345) ,@pt7))
			   (,14tmhi (,(second pts2345) ,@pt7))
			   (,41tmlo (,(third pts2345) ,@pt7))
			   (,41tmlo (,(fourth pts2345) ,@pt7))))
	   (newsxdata nil))
      (unless (gethash `(,14tmhi ,pts2345) *SL3SIMPLEX->ID*)
	(setf newsxdata 
	      `((1 ,14tmlo ,14tmhi (,@pt1 ,@pts2345))
		(4 ,41tmlo ,41tmhi (,@pts2345 ,@pt6)))) 
      (unless (or (member 1id1 *FIXED-4SXID*) (member 1id2 *FIXED-4SXID*) (member 1id2 *FIXED-4SXID*) (member 1id3 *FIXED-4SXID*) (member 1id4 *FIXED-4SXID*) (member 4id4 *FIXED-4SXID*) (member 4id3 *FIXED-4SXID*) (member 4id2 *FIXED-4SXID*) (member 4id1 *FIXED-4SXID*))
	(return-from try-8->2 (list newsxdata nbors subcx 
				    oldtlint3sxs oldslint3sxs
				    oldtlint2sxs oldslint2sxs
				    oldtlint1sxs oldslint1sxs
				    DF82)))))))


;;;--------------------------------------------------------------------------[1]
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
;;;--------------------------------------------------------------------------[1]
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
;;  (setf CURRENT-MOVE-IDENTIFIER "4->6")
  (dolist (subcx (4->6-subcomplex sxid))
    (let* ((14id1 (first subcx)) 
	   (14id2 (second subcx)) 
	   (41id2 (third subcx)) 
	   (41id1 (fourth subcx))
	   (14sx1 (get-4simplex 14id1)) 
	   (14sx2 (get-4simplex 14id2))
	   (41sx1 (get-4simplex 41id1)) 
	   (41sx2 (get-4simplex 41id2))
	   (14tmlo (4sx-tmlo 14sx1)) 
	   (14tmhi (4sx-tmhi 14sx2))
	   (41tmlo (4sx-tmlo 41sx2)) 
	   (41tmhi (4sx-tmhi 41sx1))
	   (14sx1pts (4sx-points 14sx1))
	   (14sx2pts (4sx-points 14sx2))
	   (41sx1pts (4sx-points 41sx1))
	   (41sx2pts (4sx-points 41sx2))
	   (nbors (set-difference 
		   (unions (4sx-sx4ids 14sx1) (4sx-sx4ids 14sx2) 
			   (4sx-sx4ids 41sx1) (4sx-sx4ids 41sx2))
		   (list 14id1 14id2 41id1 41id2 0)))
	   (oldTL3sxs `((1 ,14tmlo ,(intersection 14sx1pts 14sx2pts))
			(3 ,41tmlo ,(intersection 41sx1pts 41sx2pts))))
	   (oldSL3sxs `((,14tmhi ,(intersection 14sx1pts 41sx1pts))
			(,41tmlo ,(intersection 14sx2pts 41sx2pts))))
	   (pts26 (set-exclusive-or 
		   (second (first oldSL3sxs))
		   (second (second oldSL3sxs))));(2 6)
	   (pts345 (intersection (4sx-hipts 14sx1) (4sx-hipts 14sx2)));(3 4 5)
	   (pt1 (4sx-lopts 14sx1));;(1)
	   (pt7 (4sx-hipts 41sx1));;(7)
	   (oldSL2sxs `((,14tmhi ,pts345)))
	   (newsxdata nil))
      (unless (gethash `(,14tmhi ,pts26) *SL1SIMPLEX->ID*)
	(setf newsxdata `((1 ,14tmlo ,14tmhi 
			     (,@pt1 ,@pts26 ,(first pts345) ,(second pts345)))
			  (1 ,14tmlo ,14tmhi 
			     (,@pt1 ,@pts26 ,(second pts345) ,(third pts345)))
			  (1 ,14tmlo ,14tmhi 
			     (,@pt1 ,@pts26 ,(third pts345) ,(first pts345)))
			  (4 ,41tmlo ,41tmhi 
			     (,@pts26 ,(first pts345) ,(second pts345) ,@pt7))
			  (4 ,41tmlo ,41tmhi 
			     (,@pts26 ,(second pts345) ,(third pts345) ,@pt7))
			  (4 ,41tmlo ,41tmhi 
			     (,@pts26 ,(third pts345) ,(first pts345) ,@pt7))))
	(unless (or (member 14id1 *FIXED-4SXID*) (member 14id2 *FIXED-4SXID*) (member 41id2 *FIXED-4SXID*) (member 41id1 *FIXED-4SXID*))
	  (return-from try-4->6 
	    (list newsxdata nbors subcx 
		oldTL3sxs oldSL3sxs 
		nil oldSL2sxs
		nil nil 
		DF46)))))))


;;;--------------------------------------------------------------------------[8]
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
;;;--------------------------------------------------------------------------[8]
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
			 (when (and (setf sx1 (get-4simplex id1))
				    (setf sx2 (get-4simplex id2))
				    (4simplices-connected? id1 id2)
				    (4simplices-connected? 
				     (nth 0 (4sx-sx4ids sx1))
				     (nth 0 (4sx-sx4ids sx2)))
				    (4simplices-connected? 
				     (nth 0 (4sx-sx4ids sx1)) 41id)
				    (4simplices-connected? 
				     (nth 0 (4sx-sx4ids sx2)) 41id))
			   (pushnew (list sxid id1 id2 
					  (nth 0 (4sx-sx4ids sx2))
					  (nth 0 (4sx-sx4ids sx1)) 41id)
				    subcmplx :test #'set-equal?)))))))))
	    ((= 4 (4sx-type sx))
	     (let ((14id (nth 4 (4sx-sx4ids sx)))(14sx nil))
	       (when (setf 14sx (get-4simplex 14id))
		 (let ((41nbors (neighbors-of-type sx 4)))
		   (unless (< (length 41nbors) 2)
		     (do-tuples/c (id1 id2) 41nbors
		       (let ((sx1 nil) (sx2 nil))
			 (when (and (setf sx1 (get-4simplex id1)) 
				    (setf sx2 (get-4simplex id2))
				    (4simplices-connected? id1 id2)
				    (4simplices-connected? 
				     (nth 4 (4sx-sx4ids sx1))
				     (nth 4 (4sx-sx4ids sx2)))
				    (4simplices-connected? 
				     (nth 4 (4sx-sx4ids sx1)) 14id)
				    (4simplices-connected? 
				     (nth 4 (4sx-sx4ids sx2)) 14id))
			   (pushnew (list 14id (nth 4 (4sx-sx4ids sx2)) 
					  (nth 4 (4sx-sx4ids sx1)) 
					  id1 id2 sxid)
				    subcmplx :test #'set-equal?)))))))))))
    subcmplx))


(defun try-6->4 (sxid)
;;  (setf CURRENT-MOVE-IDENTIFIER "6->4")
  (dolist (subcx (6->4-subcomplex sxid))
    (let* ((14id1 (first subcx)) 
	   (14id2 (second subcx)) 
	   (14id3 (third subcx))
	   (41id3 (fourth subcx)) 
	   (41id2 (fifth subcx)) 
	   (41id1 (sixth subcx))
	   (14sx1 (get-4simplex 14id1)) 
	   (14sx2 (get-4simplex 14id2)) 
	   (14sx3 (get-4simplex 14id3))
	   (41sx1 (get-4simplex 41id1)) 
	   (41sx2 (get-4simplex 41id2)) 
	   (41sx3 (get-4simplex 41id3))
	   (14sx1pts (4sx-points 14sx1))
	   (14sx2pts (4sx-points 14sx2))
	   (14sx3pts (4sx-points 14sx3))
	   (41sx1pts (4sx-points 41sx1))
	   (41sx2pts (4sx-points 41sx2))
	   (41sx3pts (4sx-points 41sx3))
	   (14tmlo (4sx-tmlo 14sx1)) 
	   (14tmhi (4sx-tmhi 14sx2))
	   (41tmlo (4sx-tmlo 41sx2)) 
	   (41tmhi (4sx-tmhi 41sx1))
	   (nbors (set-difference 
		   (unions (4sx-sx4ids 14sx1) (4sx-sx4ids 14sx2) 
			   (4sx-sx4ids 14sx3) (4sx-sx4ids 41sx1) 
			   (4sx-sx4ids 41sx2) (4sx-sx4ids 41sx3))
		   (list 14id1 14id2 14id3 41id1 41id2 41id3 0)))
	   (oldSL3sxs `((,14tmhi ,(4sx-hipts 14sx1))
			(,14tmhi ,(4sx-hipts 14sx2))
			(,14tmhi ,(4sx-hipts 14sx3))))
	   (oldTL3sxs `((1 ,14tmlo ,(intersection 14sx1pts 14sx2pts))
			(1 ,14tmlo ,(intersection 14sx2pts 14sx3pts))
			(1 ,14tmlo ,(intersection 14sx3pts 14sx1pts))
			(4 ,41tmlo ,(intersection 41sx1pts 41sx2pts))
			(4 ,41tmlo ,(intersection 41sx2pts 41sx3pts))
			(4 ,41tmlo ,(intersection 41sx3pts 41sx1pts))))
	   (pts26 (intersections 
		   (second (first oldSL3sxs))
		   (second (second oldSL3sxs))
		   (second (third oldSL3sxs))));;(2 6)
	   (pts345 (set-difference 
		    (unions (second (first oldSL3sxs))
			    (second (second oldSL3sxs))
			    (second (third oldSL3sxs)))
		    pts26));;(3 4 5)
	   (pt1 (4sx-lopts 14sx1));;(1)
	   (pt7 (4sx-hipts 41sx1));;(7)
	   (oldTL2sxs `((1 ,14tmlo (,@pt1 ,@pts26))
			(2 ,41tmlo (,@pts26 ,@pt7))))
	   (oldSL2sxs `((,41tmlo (,(first pts345) ,@pts26))
			(,41tmlo (,(second pts345) ,@pts26))
			(,41tmlo (,(third pts345) ,@pts26))))
	   (oldSL1sxs `((,41tmlo ,pts26)))
	   (newsxdata nil))
      (unless (gethash `(,14tmhi ,pts345) *SL2SIMPLEX->ID*)
	(setf newsxdata 
	      `((1 ,14tmlo ,14tmhi (,@pt1 ,(first pts26) ,@pts345))
		(1 ,14tmlo ,14tmhi (,@pt1 ,(second pts26) ,@pts345))
		(4 ,41tmlo ,41tmhi (,(first pts26) ,@pts345 ,@pt7))
		(4 ,41tmlo ,41tmhi (,(second pts26) ,@pts345 ,@pt7))))
	(unless (or (member 14id1 *FIXED-4SXID*) (member 14id2 *FIXED-4SXID*) (member 14id3 *FIXED-4SXID*) (member 41id3 *FIXED-4SXID*) (member 41id2 *FIXED-4SXID*) (member 41id1 *FIXED-4SXID*))
	  (return-from try-6->4 
	    (list newsxdata nbors subcx 
		oldTL3sxs oldSL3sxs
		oldTL2sxs oldSL2sxs
		nil oldSL1sxs
		DF64)))))))


;;;--------------------------------------------------------------------------[2]
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
;;;--------------------------------------------------------------------------[2]
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
;;  (setf CURRENT-MOVE-IDENTIFIER "2->4-v1-internal-12")
  (let ((14sx nil) (23sx nil))
    (when (and (setf 14sx (get-4simplex 14id)) 
	       (setf 23sx (get-4simplex 23id)))
      (let* ((pts345 (4sx-hipts 23sx));;(3 4 5)
	     (pt6 (set-difference (4sx-hipts 14sx) pts345));;(6)
	     (pts12 (4sx-lopts 23sx));;(1 2)
	     (pt1 (4sx-lopts 14sx));;(1)
	     (pt2 (set-difference pts12 pt1));;(2)
	     (pts26 `(,@pt2 ,@pt6))
	     (tmlo (4sx-tmlo 14sx)) 
	     (tmhi (4sx-tmhi 23sx))
	     (nbors (set-difference (union (4sx-sx4ids 14sx) (4sx-sx4ids 23sx))
				    (list 14id 23id 0)))
	     (oldTL3sxs 
	      `((1 ,tmlo ,(intersection (4sx-points 14sx) (4sx-points 23sx)))))
	     (newsxdata nil))
	(unless (gethash `(1 ,tmlo ,pts26) *TL1SIMPLEX->ID*)
	  (setf newsxdata 
	      `((1 ,tmlo ,tmhi (,@pt2 ,@(4sx-hipts 14sx)))
		(2 ,tmlo ,tmhi (,@pts12 ,@pt6 ,(first pts345) ,(second pts345)))
		(2 ,tmlo ,tmhi (,@pts12 ,@pt6 ,(second pts345) ,(third pts345)))
		(2 ,tmlo ,tmhi (,@pts12 ,@pt6 ,(third pts345) ,(first pts345)))))
	  (unless (or (member 14id *FIXED-4SXID*) (member 23id *FIXED-4SXID*))
	    (return-from 2->4-v1-internal-12 
	      (list newsxdata nbors (list 14id 23id) 
		  oldTL3sxs nil
		  nil nil 
		  nil nil
		  DF24))))))))


(defun 2->4-v1-internal-43 (41id 32id)
;;  (setf CURRENT-MOVE-IDENTIFIER "2->4-v1-internal-43")
  (let ((41sx nil) (32sx nil))
    (when (and (setf 41sx (get-4simplex 41id)) 
	       (setf 32sx (get-4simplex 32id)))
      (let* ((pts345 (4sx-lopts 32sx));;(3 4 5)
	     (pt6 (set-difference (4sx-lopts 41sx) pts345));;(6)
	     (pts12 (4sx-hipts 32sx));;(1 2)
	     (pt1 (4sx-hipts 41sx));;(1)
	     (pt2 (set-difference pts12 pt1));;(2)
	     (pts62 `(,@pt6 ,@pt2))
	     (tmlo (4sx-tmlo 41sx)) 
	     (tmhi (4sx-tmhi 32sx))
	     (oldTL3sxs 
	      `((3 ,tmlo ,(intersection (4sx-points 41sx) (4sx-points 32sx)))))
	     (nbors (set-difference (union (4sx-sx4ids 41sx) (4sx-sx4ids 32sx))
				    (list 41id 32id 0)))
	     (newsxdata nil))
	(unless (gethash `(1 ,tmlo ,pts62) *TL1SIMPLEX->ID*)
	  (setf newsxdata 
	      `((4 ,tmlo ,tmhi (,@(4sx-lopts 41sx) ,@pt2))
		(3 ,tmlo ,tmhi (,(first pts345) ,(second pts345) ,@pt6 ,@pts12))
		(3 ,tmlo ,tmhi (,(second pts345) ,(third pts345) ,@pt6 ,@pts12))
		(3 ,tmlo ,tmhi (,(third pts345) ,(first pts345) ,@pt6 ,@pts12))))
	  (unless (or (member 41id *FIXED-4SXID*) (member 32id *FIXED-4SXID*))
	    (return-from 2->4-v1-internal-43 
	      (list newsxdata nbors (list 41id 32id) 
		  oldTL3sxs nil 
		  nil nil
		  nil nil
		  DF24))))))))


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


;;;--------------------------------------------------------------------------[7]
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
;;;
;;; or
;;;
;;; (3 4 5|1 2) + (3 4 6|1 2) + (3 5 6|1 2)
;;; +
;;; (3 4 5 6|2)
;;; ->
;;; (3 4 5 6|1)
;;; +
;;; (4 5 6|1 2)


;;;--------------------------------------------------------------------------[7]
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
				(when (and (setf sx1 (get-4simplex id1)) 
					   (setf sx2 (get-4simplex id2))
					   (setf sx3 (get-4simplex id3)) 
					   (4simplices-connected? id1 id2)
					   (4simplices-connected? id2 id3) 
					   (4simplices-connected? id3 id1))
				  (pushnew (list 1222 sxid id1 id2 id3) 
					   subcmplx :test #'set-equal?)))))))
	    ((= 2 (4sx-type sx))
	     (let ((23nbors (neighbors-of-type sx 2))
		   (14nbors (neighbors-of-type sx 1)))
	       (unless (or (< (length 23nbors) 2) (< (length 14nbors) 1))
		 (do-tuples/c (id1 id2) 23nbors
			      (let ((sx1 nil) (sx2 nil))
				(when (and (setf sx1 (get-4simplex id1)) 
					   (setf sx2 (get-4simplex id2))
					   (4simplices-connected? id1 id2))
				  (dolist (14nbor 14nbors)
				    (when (and 
					   (4simplices-connected? id1 14nbor)
					   (4simplices-connected? id2 14nbor))
				      (pushnew (list 1222 14nbor sxid id1 id2) 
					       subcmplx :test #'set-equal?)))))))))
	    ((= 3 (4sx-type sx))
	     (let ((32nbors (neighbors-of-type sx 3))
		   (41nbors (neighbors-of-type sx 4)))
	       (unless (or (< (length 32nbors) 2) (< (length 41nbors) 1))
		 (do-tuples/c (id1 id2) 32nbors
			      (let ((sx1 nil) (sx2 nil))
				(when (and (setf sx1 (get-4simplex id1)) 
					   (setf sx2 (get-4simplex id2))
					   (4simplices-connected? id1 id2))
				  (dolist (41nbor 41nbors)
				    (when (and 
					   (4simplices-connected? id1 41nbor)
					   (4simplices-connected? id2 41nbor))
				      (pushnew (list 4333 41nbor sxid id1 id2) 
					       subcmplx :test #'set-equal?)))))))))
	    ((= 4 (4sx-type sx))
	     (let ((32nbors (neighbors-of-type sx 3)))
	       (unless (< (length 32nbors) 3)
		 (do-tuples/c (id1 id2 id3) 32nbors
			      (let ((sx1 nil) (sx2 nil) (sx3 nil))
				(when (and (setf sx1 (get-4simplex id1)) 
					   (setf sx2 (get-4simplex id2))
					   (setf sx3 (get-4simplex id3)) 
					   (4simplices-connected? id1 id2)
					   (4simplices-connected? id2 id3) 
					   (4simplices-connected? id3 id1))
				  (pushnew (list 4333 sxid id1 id2 id3) 
					   subcmplx :test #'set-equal?)))))))))
    subcmplx))


(defun 4->2-v1-internal-1222 (14id 23id1 23id2 23id3)
;;  (setf CURRENT-MOVE-IDENTIFIER "4->2-v1-internal-1222")
  (let ((14sx nil) (23sx1 nil) (23sx2 nil) (23sx3 nil))
    (when (and (setf 14sx (get-4simplex 14id)) 
	       (setf 23sx1 (get-4simplex 23id1))
	       (setf 23sx2 (get-4simplex 23id2)) 
	       (setf 23sx3 (get-4simplex 23id3)))
      (let* ((14sxpts (4sx-points 14sx))
	     (23sx1pts (4sx-points 23sx1))
	     (23sx2pts (4sx-points 23sx2))
	     (23sx3pts (4sx-points 23sx3))
	     (pt3 (intersections 
		   (4sx-hipts 23sx1) (4sx-hipts 23sx2) (4sx-hipts 23sx3)));;(3)
	     (pts456 (set-difference (4sx-hipts 14sx) pt3));;(4 5 6)
	     (pts12 (4sx-lopts 23sx1));;(1 2)
	     (pt2 (4sx-lopts 14sx));(2)
	     (pt1 (set-difference pts12 pt2));;(1)
	     (tmlo (4sx-tmlo 14sx)) 
	     (tmhi (4sx-tmhi 23sx1))
	     (oldTL3sxs `((1 ,tmlo ,(intersection 14sxpts 23sx1pts))
			  (1 ,tmlo ,(intersection 14sxpts 23sx2pts))
			  (1 ,tmlo ,(intersection 14sxpts 23sx3pts))
			  (2 ,tmlo ,(intersection 23sx3pts 23sx1pts))
			  (2 ,tmlo ,(intersection 23sx1pts 23sx2pts))
			  (2 ,tmlo ,(intersection 23sx2pts 23sx3pts))))
	     (oldTL2sxs `((2 ,tmlo (,@pt1 ,@pt2 ,@pt3)) 
			  (1 ,tmlo (,@pt2 ,@pt3 ,(first pts456)))
			  (1 ,tmlo (,@pt2 ,@pt3 ,(second pts456)))
			  (1 ,tmlo (,@pt2 ,@pt3 ,(third pts456)))))
	     (oldTL1sxs `((1 ,tmlo (,@pt2 ,@pt3))));;(2 3)
	     (nbors (set-difference 
		     (unions (4sx-sx4ids 14sx) (4sx-sx4ids 23sx1) 
			     (4sx-sx4ids 23sx2) (4sx-sx4ids 23sx3))
		     (list 14id 23id1 23id2 23id3 0)))
	     (newsxdata nil))
	(unless (gethash `(1 ,tmlo (,@pt1 ,@pts456)) *TL3SIMPLEX->ID*)
	  (setf newsxdata 
		`((1 ,tmlo ,tmhi (,@pt1 ,@(4sx-hipts 14sx)))
		  (2 ,tmlo ,tmhi (,@pts12 ,@pts456))))
	  (unless (or (member 14id *FIXED-4SXID*) (member 23id1 *FIXED-4SXID*) (member 23id2 *FIXED-4SXID*) (member 23id3 *FIXED-4SXID*) )
	    (return-from 4->2-v1-internal-1222 
	      (list newsxdata nbors (list 14id 23id1 23id2 23id3)
		  oldTL3sxs nil
		  oldTL2sxs nil
		  oldTL1sxs nil
		  DF42))))))))


(defun 4->2-v1-internal-4333 (41id 32id1 32id2 32id3)
;;  (setf CURRENT-MOVE-IDENTIFIER "4->2-v1-internal-4333")
  (let ((41sx nil) (32sx1 nil) (32sx2 nil) (32sx3 nil))
    (when (and (setf 41sx (get-4simplex 41id)) 
	       (setf 32sx1 (get-4simplex 32id1))
	       (setf 32sx2 (get-4simplex 32id2)) 
	       (setf 32sx3 (get-4simplex 32id3)))
      (let* ((41sxpts (4sx-points 41sx))
	     (32sx1pts (4sx-points 32sx1))
	     (32sx2pts (4sx-points 32sx2))
	     (32sx3pts (4sx-points 32sx3))
	     (pt3 (intersections 
		   (4sx-lopts 32sx1) (4sx-lopts 32sx2) (4sx-lopts 32sx3)));;(3)
	     (pts456 (set-difference (4sx-lopts 41sx) pt3));;(4 5 6)
	     (pts12 (4sx-hipts 32sx1));;(1 2)
	     (pt2 (4sx-hipts 41sx));;(2)
	     (pt1 (set-difference pts12 pt2));;(1)
	     (tmlo (4sx-tmlo 41sx)) 
	     (tmhi (4sx-tmhi 32sx1))
	     (oldTL3sxs `((3 ,tmlo ,(intersection 41sxpts 32sx1pts))
			  (3 ,tmlo ,(intersection 41sxpts 32sx2pts))
			  (3 ,tmlo ,(intersection 41sxpts 32sx3pts))
			  (2 ,tmlo ,(intersection 32sx3pts 32sx1pts))
			  (2 ,tmlo ,(intersection 32sx1pts 32sx2pts))
			  (2 ,tmlo ,(intersection 32sx2pts 32sx3pts))))
	     (oldTL2sxs `((1 ,tmlo (,@pt3 ,@pt1 ,@pt2)) 
			  (2 ,tmlo (,@pt3 ,(first pts456) ,@pt2))
			  (2 ,tmlo (,@pt3 ,(second pts456) ,@pt2))
			  (2 ,tmlo (,@pt3 ,(third pts456) ,@pt2))))
	     (oldTL1sxs `((1 ,tmlo (,@pt3 ,@pt2))));;(3 2)
	     (nbors (set-difference 
		     (unions (4sx-sx4ids 41sx) (4sx-sx4ids 32sx1) 
			     (4sx-sx4ids 32sx2) (4sx-sx4ids 32sx3))
		     (list 41id 32id1 32id2 32id3 0)))
	     (newsxdata nil))
	(unless (gethash `(3 ,tmlo (,@pts456 ,@pt1)) *TL3SIMPLEX->ID*)
	  (setf newsxdata 
		`((4 ,tmlo ,tmhi (,@(4sx-lopts 41sx) ,@pt1))
		  (3 ,tmlo ,tmhi (,@pts456 ,@pts12))))
	  (unless (or (member 41id *FIXED-4SXID*) (member 32id1 *FIXED-4SXID*) (member 32id2 *FIXED-4SXID*) (member 32id3 *FIXED-4SXID*) )
	    (return-from 4->2-v1-internal-4333 
	      (list newsxdata nbors (list 41id 32id1 32id2 32id3)
		  oldTL3sxs nil
		  oldTL2sxs nil
		  oldTL1sxs nil
		  DF42))))))))


(defun try-4->2-v1 (sxid)
  (let ((subcmplx (4->2-v1-subcomplex sxid))
	(movedata nil))
    (unless (null subcmplx)
      (dolist (curr subcmplx)
	(cond ((= 1222 (first curr))
	       (setf movedata 
		     (4->2-v1-internal-1222 
		      (second curr) (third curr) (fourth curr) (fifth curr)))
	       (when movedata
		 (return-from try-4->2-v1 movedata)))
	      ((= 4333 (first curr))
	       (setf movedata 
		     (4->2-v1-internal-4333 
		      (second curr) (third curr) (fourth curr) (fifth curr)))
	       (when movedata
		 (return-from try-4->2-v1 movedata))))))))


;;;--------------------------------------------------------------------------[3]
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
;;;--------------------------------------------------------------------------[3]
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
;;  (setf CURRENT-MOVE-IDENTIFIER "2->4-v2")
  (dolist (subcx (2->4-v2-subcomplex sxid))
    (let* ((23id (first subcx)) 
	   (32id (second subcx))
	   (23sx (get-4simplex 23id)) 
	   (32sx (get-4simplex 32id))
	   (23sxpts (4sx-points 23sx))
	   (32sxpts (4sx-points 32sx))
	   (tmlo (4sx-tmlo 23sx)) 
	   (tmhi (4sx-tmhi 32sx))
	   (nbors (set-difference 
		   (union (4sx-sx4ids 23sx) (4sx-sx4ids 32sx)) 
		   (list 23id 32id 0)))
	   (oldTL3sxs `((2 ,tmlo ,(intersection 23sxpts 32sxpts))))
	   (pts12 (4sx-lopts 23sx)) 
	   (pts456 (4sx-hipts 23sx))
	   (pts123 (4sx-lopts 32sx)) 
	   (pts45 (4sx-hipts 32sx))
	   (pt3 (set-difference pts123 pts12)) 
	   (pt6 (set-difference pts456 pts45))
	   (newsxdata nil))
      (unless (gethash `(1 ,tmlo (,@pt3 ,@pt6)) *TL1SIMPLEX->ID*)
	(setf newsxdata
	      `((2 ,tmlo ,tmhi (,(first pts12) ,@pt3 ,@pts456))
		(2 ,tmlo ,tmhi (,(second pts12) ,@pt3 ,@pts456))
		(3 ,tmlo ,tmhi (,@pts123 ,(first pts45) ,@pt6))
		(3 ,tmlo ,tmhi (,@pts123 ,(second pts45) ,@pt6))))
	(unless (or (member 23id *FIXED-4SXID*) (member 32id *FIXED-4SXID*))
	  (return-from try-2->4-v2 
	    (list newsxdata nbors subcx 
		oldTL3sxs nil 
		nil nil
		nil nil
		DF24)))))))


;;;--------------------------------------------------------------------------[6]
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
;;;--------------------------------------------------------------------------[6]
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
				(when (and (setf sx1 (get-4simplex id1)) 
					   (setf sx2 (get-4simplex id2))
					   (4simplices-connected? id1 id2))
				  (dolist (23nbor 23nbors)
				    (when (and 
					   (4simplices-connected? id1 23nbor)
					   (4simplices-connected? id2 23nbor))
				      (pushnew (list sxid 23nbor id1 id2) 
					       subcmplx :test #'set-equal?)))))))))
	    ((= 3 (4sx-type sx))
	     (let ((23nbors (neighbors-of-type sx 2))
		   (32nbors (neighbors-of-type sx 3)))
	       (unless (or (< (length 23nbors) 2) (< (length 32nbors) 1))
		 (do-tuples/c (id1 id2) 23nbors
			      (let ((sx1 nil) (sx2 nil))
				(when (and (setf sx1 (get-4simplex id1)) 
					   (setf sx2 (get-4simplex id2))
					   (4simplices-connected? id1 id2))
				  (dolist (32nbor 32nbors)
				    (when (and 
					   (4simplices-connected? id1 32nbor)
					   (4simplices-connected? id2 32nbor))
				      (pushnew (list id1 id2 32nbor sxid) 
					       subcmplx :test #'set-equal?)))))))))))
    subcmplx))


(defun try-4->2-v2 (sxid)
;;  (setf CURRENT-MOVE-IDENTIFIER "4->2-v2")
  (dolist (subcx (4->2-v2-subcomplex sxid))
    (let* ((23id1 (first subcx)) 
	   (23id2 (second subcx))
	   (32id1 (third subcx)) 
	   (32id2 (fourth subcx))
	   (23sx1 (get-4simplex 23id1)) 
	   (23sx2 (get-4simplex 23id2))
	   (32sx1 (get-4simplex 32id1)) 
	   (32sx2 (get-4simplex 32id2))
	   (23sx1pts (4sx-points 23sx1))
	   (23sx2pts (4sx-points 23sx2))
	   (32sx1pts (4sx-points 32sx1))
	   (32sx2pts (4sx-points 32sx2))
	   (tmlo (4sx-tmlo 23sx1)) 
	   (tmhi (4sx-tmhi 32sx2))
	   (nbors (set-difference 
		   (unions (4sx-sx4ids 23sx1) (4sx-sx4ids 23sx2)
			   (4sx-sx4ids 32sx1) (4sx-sx4ids 32sx2))
		   (list 23id1 23id2 32id1 32id2 0)))
	   (oldTL3sxs `((1 ,tmlo ,(intersection 23sx1pts 23sx2pts))
			(3 ,tmlo ,(intersection 32sx1pts 32sx2pts))
			(2 ,tmlo ,(intersection 23sx1pts 32sx1pts))
			(2 ,tmlo ,(intersection 23sx1pts 32sx2pts))
			(2 ,tmlo ,(intersection 23sx2pts 32sx1pts))
			(2 ,tmlo ,(intersection 23sx2pts 32sx2pts))))
	   (pts13 (4sx-lopts 23sx1)) 
	   (pts23 (4sx-lopts 23sx2)) 
	   (pts456 (4sx-hipts 23sx1))
	   (pts123 (4sx-lopts 32sx1)) 
	   (pts46 (4sx-hipts 32sx1)) 
	   (pts56 (4sx-hipts 32sx2))
	   (pt6 (intersection pts46 pts56))
	   (pt3 (intersection pts13 pts23))
	   (oldTL2sxs `((1 ,tmlo (,@pt3 ,@pts46))
			(1 ,tmlo (,@pt3 ,@pts56)) 
			(2 ,tmlo (,@pts13 ,@pt6))
			(2 ,tmlo (,@pts23 ,@pt6))))
	   (oldTL1sxs `((1 ,tmlo (,@pt3 ,@pt6))))
	   (newsxdata nil))
      (unless (gethash `(2 ,tmlo (,@(set-exclusive-or pts13 pts23) 
				    ,@(set-exclusive-or pts46 pts56)))
		       *TL3SIMPLEX->ID*)
	(setf newsxdata 
	      `((2 ,tmlo ,tmhi (,@(set-exclusive-or pts13 pts23) ,@pts456))
		(3 ,tmlo ,tmhi (,@pts123 ,@(set-exclusive-or pts46 pts56)))))
      (unless (or (member 23id1 *FIXED-4SXID*) (member 23id2 *FIXED-4SXID*) (member 32id1 *FIXED-4SXID*) (member 32id2 *FIXED-4SXID*))
	(return-from try-4->2-v2 
	    (list newsxdata nbors subcx 
		oldTL3sxs nil
		oldTL2sxs nil
		oldTL1sxs nil
		DF42)))))))


;;;--------------------------------------------------------------------------[4]
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
;;; or
;;; 
;;; (1 2 3 4|5)
;;; +
;;; (1 2 3|5 6) + (1 2 4|5 6)
;;; ->
;;; (1 2 3 4|6)
;;; +
;;; (2 3 4|5 6) + (1 3 4|5 6)
;;;--------------------------------------------------------------------------[4]
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
				(when (and (setf sx1 (get-4simplex id1)) 
					   (setf sx2 (get-4simplex id2))
					   (4simplices-connected? id1 id2))
				  (pushnew (list 122 sxid id1 id2) 
					   subcmplx :test #'set-equal?)))))))
	    ((= 2 (4sx-type sx))
	     (let ((23nbors (neighbors-of-type sx 2))
		   (14nbors (neighbors-of-type sx 1)))
	       (dolist (23nbor 23nbors)
		 (dolist (14nbor 14nbors)
		   (let ((sx1 nil) (sx2 nil))
		     (when (and (setf sx1 (get-4simplex 23nbor)) 
				(setf sx2 (get-4simplex 14nbor))
				(4simplices-connected? 23nbor 14nbor))
		       (pushnew (list 122 14nbor sxid 23nbor) 
				subcmplx :test #'set-equal?)))))))
	    ((= 3 (4sx-type sx))
	     (let ((32nbors (neighbors-of-type sx 3))
		   (41nbors (neighbors-of-type sx 4)))
	       (dolist (32nbor 32nbors)
		 (dolist (41nbor 41nbors)
		   (let ((sx1 nil) (sx2 nil))
		     (when (and (setf sx1 (get-4simplex 32nbor)) 
				(setf sx2 (get-4simplex 41nbor))
				(4simplices-connected? 32nbor 41nbor))
		       (pushnew (list 433 41nbor sxid 32nbor) 
				subcmplx :test #'set-equal?)))))))
	    ((= 4 (4sx-type sx))
	     (let ((32nbors (neighbors-of-type sx 3)))
	       (unless (< (length 32nbors) 2)
		 (do-tuples/c (id1 id2) 32nbors
			      (let ((sx1 nil) (sx2 nil))
				(when (and (setf sx1 (get-4simplex id1)) 
					   (setf sx2 (get-4simplex id2))
					   (4simplices-connected? id1 id2))
				  (pushnew (list 433 sxid id1 id2) 
					   subcmplx :test #'set-equal?)))))))))
    subcmplx))


(defun 3->3-v1-internal-122 (14id 23id1 23id2)
;;  (setf CURRENT-MOVE-IDENTIFIER "3->3-v1-internal-122")
  (let ((14sx nil) (23sx1 nil) (23sx2 nil))
    (when (and (setf 14sx (get-4simplex 14id)) 
	       (setf 23sx1 (get-4simplex 23id1))
	       (setf 23sx2 (get-4simplex 23id2)))
      (let* ((14sxpts (4sx-points 14sx))
	     (23sx1pts (4sx-points 23sx1))
	     (23sx2pts (4sx-points 23sx2))
	     (pts12 (4sx-lopts 23sx1))
	     (pt1 (4sx-lopts 14sx))
	     (pt2 (set-difference pts12 pt1))
	     (pts34 (intersection (4sx-hipts 23sx1) (4sx-hipts 23sx2)))
	     (pts56 (set-exclusive-or (4sx-hipts 23sx1) (4sx-hipts 23sx2)))
	     (tmlo (4sx-tmlo 23sx1)) (tmhi (4sx-tmhi 23sx2))
	     (nbors (set-difference 
		     (unions (4sx-sx4ids 14sx) (4sx-sx4ids 23sx1)
			     (4sx-sx4ids 23sx2))
		     (list 14id 23id1 23id2 0)))
	     (oldTL3sxs `((1 ,tmlo ,(intersection 14sxpts 23sx1pts))
			  (1 ,tmlo ,(intersection 14sxpts 23sx2pts))
			  (2 ,tmlo ,(intersection 23sx1pts 23sx2pts))))
	     (oldTL2sxs `((1 ,tmlo (,@pt1 ,@pts34))))
	     (newsxdata nil))
	(unless (gethash `(1 ,tmlo (,@pt2 ,@pts56)) *TL2SIMPLEX->ID*)
	  (setf newsxdata 
		`((1 ,tmlo ,tmhi (,@pt2 ,@pts34 ,@pts56))
		  (2 ,tmlo ,tmhi (,@pts12 ,(first pts34) ,@pts56))
		  (2 ,tmlo ,tmhi (,@pts12 ,(second pts34) ,@pts56))))
	  (unless (or (member 14id *FIXED-4SXID*) (member 23id1 *FIXED-4SXID*) (member 23id2 *FIXED-4SXID*))
 	    (return-from 3->3-v1-internal-122 
	      (list newsxdata nbors (list 14id 23id1 23id2) 
		  oldTL3sxs nil
		  oldTL2sxs nil 
		  nil nil
		  DF33))))))))


(defun 3->3-v1-internal-433 (41id 32id1 32id2)
;;  (setf CURRENT-MOVE-IDENTIFIER "3->3-v1-internal-433")
  (let ((41sx nil) (32sx1 nil) (32sx2 nil))
    (when (and (setf 41sx (get-4simplex 41id)) 
	       (setf 32sx1 (get-4simplex 32id1))
	       (setf 32sx2 (get-4simplex 32id2)))
      (let* ((41sxpts (4sx-points 41sx))
	     (32sx1pts (4sx-points 32sx1))
	     (32sx2pts (4sx-points 32sx2))
	     (pts56 (4sx-hipts 32sx1))
	     (pt5 (4sx-hipts 41sx))
	     (pt6 (set-difference pts56 pt5))
	     (pts12 (intersection (4sx-lopts 32sx1) (4sx-lopts 32sx2)))
	     (pts34 (set-exclusive-or (4sx-lopts 32sx1) (4sx-lopts 32sx2)))
	     (tmlo (4sx-tmlo 32sx2)) 
	     (tmhi (4sx-tmhi 32sx1))
	     (nbors (set-difference (unions (4sx-sx4ids 41sx)(4sx-sx4ids 32sx1) 
					    (4sx-sx4ids 32sx2))
				    (list 41id 32id1 32id2 0)))
	     (oldTL3sxs `((3 ,tmlo ,(intersection 41sxpts 32sx1pts))
			  (3 ,tmlo ,(intersection 41sxpts 32sx2pts))
			  (2 ,tmlo ,(intersection 32sx1pts 32sx2pts))))
	     (oldTL2sxs `((2 ,tmlo (,@pts12 ,@pt5))))
	     (newsxdata nil))
	(unless (gethash `(2 ,tmlo (,@pts34 ,@pt6)) *TL2SIMPLEX->ID*)
	  (setf newsxdata 
		`((4 ,tmlo ,tmhi (,@pts34 ,@pts12 ,@pt6))
		  (3 ,tmlo ,tmhi (,@pts34 ,(first pts12) ,@pts56))
		  (3 ,tmlo ,tmhi (,@pts34 ,(second pts12) ,@pts56))))
	  (unless (or (member 41id *FIXED-4SXID*) (member 32id1 *FIXED-4SXID*) (member 32id2 *FIXED-4SXID*))
	  (return-from 3->3-v1-internal-433 
	    (list newsxdata nbors (list 41id 32id1 32id2) 
		  oldTL3sxs nil
		  oldTL2sxs nil 
		  nil nil
		  DF33))))))))


(defun try-3->3-v1 (sxid)
  (let ((subcmplx (3->3-v1-subcomplex sxid))
	(movedata nil))
    (unless (null subcmplx)
      (dolist (curr subcmplx)
	(cond ((= 122 (first curr))
	       (setf movedata (3->3-v1-internal-122 
			       (second curr) (third curr) (fourth curr)))
	       (when movedata
		 (return-from try-3->3-v1 movedata)))
	      ((= 433 (first curr))
	       (setf movedata (3->3-v1-internal-433 
			       (second curr) (third curr) (fourth curr)))
	       (when movedata
		 (return-from try-3->3-v1 movedata))))))))


;;;--------------------------------------------------------------------------[5]
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
;;;
;;; (1 2|4 5 6) + (1 3|4 5 6)
;;; +
;;; (1 2 3|4 5)
;;; ->
;;; (2 3|4 5 6)
;;; +
;;; (1 2 3|4 6) + (1 2 3|5 6)
;;;
;;;
;;; 2(3,2) + 1(2,3) -> 1(3,2) + 2(2,3)
;;;
;;; (1 2 3|4 6) + (1 2 3|5 6)
;;; +
;;; (2 3|4 5 6)
;;; ->
;;; (1 2 3|4 5)
;;; +
;;; (1 2|4 5 6) + (1 3|4 5 6)
;;;--------------------------------------------------------------------------[5]
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
		     (pushnew (list 223 sxid 23nbor 32nbor) 
			      subcmplx :test #'set-equal?))))
	       (when (>= (length 32nbors) 2)
		 (do-tuples/c (id1 id2) 32nbors
			      (when (4simplices-connected? id1 id2)
				(pushnew (list 332 id1 id2 sxid) 
					 subcmplx :test #'set-equal?)))))
	      ((= 3 (4sx-type sx))
	       (dolist (32nbor 32nbors)
		 (dolist (23nbor 23nbors)
		   (when (4simplices-connected? 32nbor 23nbor)
		     (pushnew (list 332 sxid 32nbor 23nbor) 
			      subcmplx :test #'set-equal?))))
	       (when (>= (length 23nbors) 2)
		 (do-tuples/c (id1 id2) 23nbors
			      (when (4simplices-connected? id1 id2)
				(pushnew (list 223 id1 id2 sxid) 
					 subcmplx :test #'set-equal?))))))))
    subcmplx))


(defun 3->3-v2-internal-223 (23id1 23id2 32id)
;;  (setf CURRENT-MOVE-IDENTIFIER "3->3-v2-internal-223")
  (let ((23sx1 nil) (23sx2 nil) (32sx nil))
    (when (and (setf 23sx1 (get-4simplex 23id1)) 
	       (setf 23sx2 (get-4simplex 23id2))
	       (setf 32sx (get-4simplex 32id)))
      (let* ((23sx1pts (4sx-points 23sx1))
	     (23sx2pts (4sx-points 23sx2))
	     (32sxpts (4sx-points 32sx))
	     (pts123 (4sx-lopts 32sx))
	     (pts23 (set-exclusive-or (4sx-lopts 23sx1) (4sx-lopts 23sx2)))
	     (pts456 (4sx-hipts 23sx1))
	     (pts45 (4sx-hipts 32sx))
	     (pt6 (set-difference pts456 pts45))
	     (pt1 (set-difference pts123 pts23))
	     (tmlo (4sx-tmlo 23sx1)) 
	     (tmhi (4sx-tmhi 23sx2))
	     (nbors (set-difference 
		     (unions (4sx-sx4ids 23sx1)
			     (4sx-sx4ids 23sx2) (4sx-sx4ids 32sx))
		     (list 23id1 23id2 32id 0)))
	     (oldTL3sxs `((1 ,tmlo ,(intersection 23sx1pts 23sx2pts))
			  (2 ,tmlo ,(intersection 23sx1pts 32sxpts))
			  (2 ,tmlo ,(intersection 23sx2pts 32sxpts))))
	     (oldTL2sxs `((1 ,tmlo (,@pt1 ,@pts45))))
	     (newsxdata nil))
	(unless (gethash `(2 ,tmlo (,@pts23 ,@pt6)) *TL2SIMPLEX->ID*)
	  (setf newsxdata 
		`((2 ,tmlo ,tmhi (,@pts23 ,@pts456))
		  (3 ,tmlo ,tmhi (,@pts123 ,(first pts45) ,@pt6))
		  (3 ,tmlo ,tmhi (,@pts123 ,(second pts45) ,@pt6))))
	  (unless (or (member 23id1 *FIXED-4SXID*) (member 23id2 *FIXED-4SXID*) (member 32id *FIXED-4SXID*))
	    (return-from 3->3-v2-internal-223 
	      (list newsxdata nbors (list 23id1 23id2 32id) 
		  oldTL3sxs nil
		  oldTL2sxs nil 
		  nil nil
		  DF33))))))))


(defun 3->3-v2-internal-332 (32id1 32id2 23id)
;;  (setf CURRENT-MOVE-IDENTIFIER "3->3-v2-internal-332")
  (let ((32sx1 nil) (32sx2 nil) (23sx nil))
    (when (and (setf 32sx1 (get-4simplex 32id1)) 
	       (setf 32sx2 (get-4simplex 32id2))
	       (setf 23sx (get-4simplex 23id)))
      (let* ((32sx1pts (4sx-points 32sx1))
	     (32sx2pts (4sx-points 32sx2))
	     (23sxpts (4sx-points 23sx))
	     (pts456 (4sx-hipts 23sx))
	     (pts45 (set-exclusive-or (4sx-hipts 32sx1) (4sx-hipts 32sx2)))
	     (pt6 (set-difference pts456 pts45))
	     (pts123 (4sx-lopts 32sx1))
	     (pts23 (4sx-lopts 23sx))
	     (pt1 (set-difference pts123 pts23))
	     (tmlo (4sx-tmlo 32sx1)) 
	     (tmhi (4sx-tmhi 32sx2))
	     (nbors (set-difference 
		     (unions (4sx-sx4ids 32sx1) (4sx-sx4ids 32sx2) 
			     (4sx-sx4ids 23sx))
		     (list 32id1 32id2 23id 0)))
	     (oldTL3sxs `((3 ,tmlo ,(intersection 32sx1pts 32sx2pts))
			  (2 ,tmlo ,(intersection 32sx1pts 23sxpts))
			  (2 ,tmlo ,(intersection 32sx2pts 23sxpts))))
	     (oldTL2sxs `((2 ,tmlo (,@pts23 ,@pt6))))
	     (newsxdata nil))
	(unless (gethash `(1 ,tmlo (,@pt1 ,@pts45)) *TL2SIMPLEX->ID*)
	  (setf newsxdata 
		`((3 ,tmlo ,tmhi (,@pts123 ,@pts45))
		  (2 ,tmlo ,tmhi (,@pt1 ,(first pts23) ,@pts456))
		  (2 ,tmlo ,tmhi (,@pt1 ,(second pts23) ,@pts456))))
	  (unless (or (member 32id1 *FIXED-4SXID*) (member 32id2 *FIXED-4SXID*) (member 23id *FIXED-4SXID*))
	    (return-from 3->3-v2-internal-332 
	      (list newsxdata nbors (list 32id1 32id2 23id) 
		  oldTL3sxs nil
		  oldTL2sxs nil
		  nil nil
		  DF33))))))))


(defun try-3->3-v2 (sxid)
  (let ((subcmplx (3->3-v2-subcomplex sxid))
	(movedata nil))
    (unless (null subcmplx)
      (dolist (curr subcmplx)
	(cond ((= 223 (first curr))
	       (setf movedata (3->3-v2-internal-223 
			       (second curr) (third curr) (fourth curr)))
	       (when movedata
		 (return-from try-3->3-v2 movedata)))
	      ((= 332 (first curr))
	       (setf movedata (3->3-v2-internal-332 
			       (second curr) (third curr) (fourth curr)))
	       (when movedata
		 (return-from try-3->3-v2 movedata))))))))

