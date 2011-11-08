;; try-a->b methods returns the following list, IFF the move can be 
;; successfully made (new3sxdata nbors old3sxids old2sxids fvector)
;; Otherwise it returns nil

;; The parameters are as follows
;; new3sxdata : The data necessary to construct the new simplices
;;            : ((type, tmlo, tmhi, (points)), ...)
;; nbors : The 3simplex ids that neighbor the subcomplex that was modified
;; old3sxids : The 3simplices in the modified subcomplex
;; old2sxids : The triangles internal to the modified subcomplex
;; fvector : the delta vector for our modification
;; point id : for the point, if any, created or deleted

(defmacro movedata-new3sxdata (movedata)
  `(nth 0 ,movedata))

(defmacro movedata-nbors (movedata)
  `(nth 1 ,movedata))

(defmacro movedata-old3sxids (movedata)
  `(nth 2 ,movedata))

(defmacro movedata-old2sxids (movedata)
  `(nth 3 ,movedata))

(defmacro movedata-fvector (movedata)
  `(nth 4 ,movedata))

(defmacro movedata-pointid (movedata)
  `(nth 5 ,movedata))

;; 2011-04-12 
;; sxdata should be changed to 
;; new3sxids 3sxnbors old3sxids old2sxids fvect news2sxids s2sxnbors olds2sxids
;; the last 3 elements will be null for 2->3 and 3->2 moves because new 
;; spatial triangles come into existence and old ones are deleted for the
;; 2->6, 6->2 and the 4->4 moves

(defun 2plus1move (mtype sxdata)
  ;; *Went here before
;;  (when (= mtype 26MTYPE)
;;      (setf (gethash (movedata-pointid sxdata) *0SIMPLEX->NS2SX*) 
;;	    (list (2sx-tmlo (gethash (car (movedata-old2sxids sxdata)) *ID->2SIMPLEX*)) '(0))))
;;  (when (= mtype 62MTYPE)
;;      (remhash (movedata-pointid sxdata) *0SIMPLEX->NS2SX*))
  (let ((new3sxids (make-3simplices-in-bulk (first sxdata))))
    (connect-3simplices-within-list new3sxids)
    (connect-3simplices-across-lists new3sxids (second sxdata))
    ;; * Next two lines moved from previous position
    (remove-3simplices (third sxdata))
    (remove-2simplices (fourth sxdata))
    (update-f-vector (fifth sxdata))))

;;
;; The code that suppports 2->6 moves
;;

;; Looks at the chosen simplex and identifies all the useful subcomplexes
;; In our case, there's only one relevant subcomplex that you can possibly
;; be part of

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

;; We'll try the operation on each of the supplied subcomplexes, hoping 
;; that it works on one. As soon as it works, we'll break and return
;; the results.
			  
(defun try-2->6 (sxid)
  (dolist (curr (2->6-subcomplex sxid))
    (let* ((sx13 (get-3simplex (first curr)))
	   (sx31 (get-3simplex (second curr)))
	   (old-internal-triangle (nth 0 (3sx-sx2ids sx13)))
	   (nbors (set-difference (union (3sx-sx3ids sx13) (3sx-sx3ids sx31)) 
				  curr))
	   (lopt (nth-point sx13 0))
	   (hipt (nth-point sx31 3))
	   (bigtr (2sx-points (get-2simplex old-internal-triangle)))
	   (13tmlo (3sx-tmlo sx13)) (13tmhi (3sx-tmhi sx13))
	   (31tmlo (3sx-tmlo sx31)) (31tmhi (3sx-tmhi sx31))
	   (newpt (next-pt))
	   (newsxdata 
	    (list (list 1 13tmlo 13tmhi 
			(list lopt newpt (first bigtr) (second bigtr)))
		  (list 1 13tmlo 13tmhi 
			(list lopt newpt (second bigtr) (third bigtr)))
		  (list 1 13tmlo 13tmhi 
			(list lopt newpt (third bigtr) (first bigtr)))
		  (list 3 31tmlo 31tmhi 
			(list (first bigtr) (second bigtr) newpt hipt))
		  (list 3 31tmlo 31tmhi 
			(list (second bigtr) (third bigtr) newpt hipt))
		  (list 3 31tmlo 31tmhi 
			(list (third bigtr) (first bigtr) newpt hipt)))))
      (return-from try-2->6 (list newsxdata nbors curr 
				  (list old-internal-triangle) DF26 newpt)))))

;;
;; Code that supports a 6->2 move
;;

(defun 6->2-subcomplex (sxid)
  "returns a list of the form ((13id1 13id2 13id3 31id3 31id2 31id1)...)"
  (let ((sx nil)
	(subcmplx nil))
    (when (setf sx (get-3simplex sxid))
      (cond ((= 3SX13 (3sx-type sx))
	     (when (= (elt *SVOLS* (3sx-tmhi sx)) 4)
	       (return-from 6->2-subcomplex))
	     (let ((31id (nth 0 (3sx-sx3ids sx)))
		   (31sx nil))
	       (when (setf 31sx (get-3simplex 31id))
		 (let ((13nbors (neighbors-of-type sx 1)))
		   (unless (< (length 13nbors) 2)
		     (do-tuples/c (currid nextid) 13nbors
		       (let ((curr nil) (next nil))
			 (when (and (setf curr (get-3simplex currid)) (setf next (get-3simplex nextid))
				    (3simplices-connected? currid nextid)
				    (3simplices-connected? (nth 0 (3sx-sx3ids curr))
							   (nth 0 (3sx-sx3ids next)))
				    (3simplices-connected? (nth 0 (3sx-sx3ids curr)) 31id)
				    (3simplices-connected? (nth 0 (3sx-sx3ids next)) 31id))
			   (pushnew (list sxid currid nextid 
					  (nth 0 (3sx-sx3ids next)) 
					  (nth 0 (3sx-sx3ids curr))
					  31id)
				    subcmplx :test #'set-equal?)))))))))
	    ((= 3SX31 (3sx-type sx))
	     (when (= (elt *SVOLS* (3sx-tmlo sx)) 4)
	       (return-from 6->2-subcomplex))
	     (let ((13id (nth 3 (3sx-sx3ids sx)))
		   (13sx nil))
	       (when (setf 13sx (get-3simplex 13id))
		 (let ((31nbors (neighbors-of-type sx 3)))
		   (unless (< (length 31nbors) 2)
		     (do-tuples/c (currid nextid) 31nbors
		       (let ((curr nil) (next nil))
			 (when (and (setf curr (get-3simplex currid)) (setf next (get-3simplex nextid))
				    (3simplices-connected? currid nextid)
				    (3simplices-connected? (nth 3 (3sx-sx3ids curr))
							   (nth 3 (3sx-sx3ids next)))
				    (3simplices-connected? (nth 3 (3sx-sx3ids curr)) 13id)
				    (3simplices-connected? (nth 3 (3sx-sx3ids next)) 13id))
			   (pushnew (list 13id (nth 3 (3sx-sx3ids next)) (nth 3 (3sx-sx3ids curr))
					  currid nextid sxid) 
				    subcmplx :test #'set-equal?)))))))))))
    subcmplx))
	    
(defun try-6->2 (sxid)
  (let ((subcmplx (6->2-subcomplex sxid))
	(old-internal-triangles nil)
	(new-internal-triangle nil)
	(old-pt nil)
	(nbors nil)
	(newsxdata nil)
	(lopt nil)
	(hipt nil)
	(13tmlo nil) (13tmhi nil) (31tmlo nil) (31tmhi nil))
    (unless (null subcmplx)
      (dolist (curr subcmplx)
	(setf old-internal-triangles (list (link-id (first curr) (sixth curr))
					   (link-id (first curr) (second curr))
					   (link-id (first curr) (third curr))
					   (link-id (second curr) (third curr))
					   (link-id (second curr) (fifth curr))
					   (link-id (third curr) (fourth curr))
					   (link-id (fourth curr) (fifth curr))
					   (link-id (fourth curr) (sixth curr))
					   (link-id (fifth curr) (sixth curr))))
	(setf nbors (set-difference (unions (3sx-sx3ids (get-3simplex (first curr)))
					    (3sx-sx3ids (get-3simplex (second curr)))
					    (3sx-sx3ids (get-3simplex (third curr)))
					    (3sx-sx3ids (get-3simplex (fourth curr)))
					    (3sx-sx3ids (get-3simplex (fifth curr)))
					    (3sx-sx3ids (get-3simplex (sixth curr))))
				    (list 0 (first curr) (second curr) (third curr) 
					  (fourth curr) (fifth curr) (sixth curr))))
	(setf 13tmlo (3sx-tmlo (get-3simplex (first curr))))
	(setf 13tmhi (3sx-tmhi (get-3simplex (first curr))))
	(setf 31tmlo (3sx-tmlo (get-3simplex (sixth curr))))
	(setf 31tmhi (3sx-tmhi (get-3simplex (sixth curr))))
	(setf lopt (nth-point (get-3simplex (first curr)) 0))
	(setf hipt (nth-point (get-3simplex (sixth curr)) 3))
	(setf old-pt (car (intersections 
			   (2sx-points (get-2simplex (first old-internal-triangles)))
			   (2sx-points (get-2simplex (fifth old-internal-triangles)))
			   (2sx-points (get-2simplex (sixth old-internal-triangles))))))
	(setf new-internal-triangle 
	      (list 0 13tmlo 13tmhi
		    (set-difference 
		     (unions (2sx-points (get-2simplex (first old-internal-triangles)))
			     (2sx-points (get-2simplex (fifth old-internal-triangles)))
			     (2sx-points (get-2simplex (sixth old-internal-triangles))))
		     (list old-pt))))
	(setf newsxdata (list (list 1 13tmlo 13tmhi (cons lopt (2sx-points new-internal-triangle))) 
			      (list 3 31tmlo 31tmhi (append (2sx-points new-internal-triangle) (list hipt)))))
	(return-from try-6->2 (list newsxdata nbors curr old-internal-triangles DF62 old-pt))))))
  
;;
;; Code that supports a 4->4 move
;;

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
		     (when (3simplices-connected? (nth 0 (3sx-sx3ids (get-3simplex 13nbor))) 31id)
		       (pushnew (list sxid 13nbor (nth 0 (3sx-sx3ids (get-3simplex 13nbor))) 31id)
				subcmplx :test #'set-equal?)))))))
	    ((= 3 (3sx-type sx))
	     (let ((13id (nth 3 (3sx-sx3ids sx))) (13sx nil))
	       (when (setf 13sx (get-3simplex 13id))
		 (let ((31nbors (neighbors-of-type sx 3)))
		   (dolist (31nbor 31nbors)
		     (when (3simplices-connected? (nth 3 (3sx-sx3ids (get-3simplex 31nbor))) 13id)
		       (pushnew (list 13id (nth 3 (3sx-sx3ids (get-3simplex 31nbor))) 31nbor sxid)
				subcmplx :test #'set-equal?)))))))))
    subcmplx))

(defun try-4->4 (sxid)
  (let ((subcmplx (4->4-subcomplex sxid))
	(old-internal-triangles nil)
	(new-internal-sl-triangle-1 nil) (new-internal-sl-triangle-2 nil)
	(new-internal-tl-triangle-1 nil) (new-internal-tl-triangle-2 nil)
	(shared nil) (unshared nil) (nbors nil) (newsxdata nil) (lopt nil) (hipt nil)
	(13tmlo nil) (13tmhi nil) (31tmlo nil) (31tmhi nil))
    (unless (null subcmplx)
      (dolist (curr subcmplx)
	(setf old-internal-triangles (list (link-id (first curr) (fourth curr)) ;; spacelike
					   (link-id (second curr) (third curr)) ;; spacelike
					   (link-id (first curr) (second curr)) ;; timelike
					   (link-id (third curr) (fourth curr)))) ;; timelike
	(setf 13tmlo (3sx-tmlo (get-3simplex (first curr))))
	(setf 13tmhi (3sx-tmhi (get-3simplex (first curr))))
	(setf 31tmlo (3sx-tmlo (get-3simplex (fourth curr))))
	(setf 31tmhi (3sx-tmhi (get-3simplex (fourth curr))))
	(setf lopt (nth-point (get-3simplex (first curr)) 0))
	(setf hipt (nth-point (get-3simplex (fourth curr)) 3))
	(setf nbors (set-difference (unions (3sx-sx3ids (get-3simplex (first curr)))
					    (3sx-sx3ids (get-3simplex (second curr)))
					    (3sx-sx3ids (get-3simplex (third curr)))
					    (3sx-sx3ids (get-3simplex (fourth curr))))
				    (list 0 (first curr) (second curr) (third curr) (fourth curr)))) 
	(setf shared (intersection (2sx-points (get-2simplex (first old-internal-triangles)))
				   (2sx-points (get-2simplex (second old-internal-triangles)))))
	(setf unshared (set-exclusive-or (2sx-points (get-2simplex (first old-internal-triangles)))
					 (2sx-points (get-2simplex (second old-internal-triangles)))))
	(setf new-internal-sl-triangle-1 (list 0 13tmlo 13tmhi (append unshared (butlast shared))))
	(setf new-internal-sl-triangle-2 (list 0 13tmlo 13tmhi (append unshared (last shared))))
	(setf new-internal-tl-triangle-1 (list 1 13tmlo 13tmhi (cons lopt unshared)))
	(setf new-internal-tl-triangle-2 (list 2 31tmlo 31tmhi (append unshared (cons hipt nil))))
	
	(unless (gethash unshared *1SIMPLEX?*)
	  (setf newsxdata (list (list 1 13tmlo 13tmhi (cons lopt (2sx-points new-internal-sl-triangle-1)))
				(list 1 13tmlo 13tmhi (cons lopt (2sx-points new-internal-sl-triangle-2)))
				(list 3 31tmlo 31tmhi (append (2sx-points new-internal-sl-triangle-1) 
							      (list hipt)))
				(list 3 31tmlo 31tmhi (append (2sx-points new-internal-sl-triangle-2) 
							      (list hipt)))))
	  (return-from try-4->4 (list newsxdata nbors curr old-internal-triangles DF44 nil)))))))

;;
;; Code that supports a 2->3 move
;;

(defun 2->3-subcomplex (sxid)
  "returns a list of the form ((1or3 13or31id 22id)...) where the first number 1 or 3 tells us about the
type of the simplex participating in the move"
  (let ((sx nil)
	(subcmplx nil))
    (when (setf sx (get-3simplex sxid))
      (cond ((or (= 1 (3sx-type sx)) (= 3 (3sx-type sx)))
	     (let ((22nbors (neighbors-of-type sx 2)))
	       (dolist (22nbor 22nbors)
		 (pushnew (list (3sx-type sx) sxid 22nbor) subcmplx :test #'set-equal?))))
	    ((= 2 (3sx-type sx))
	     (let ((13nbors (append (neighbors-of-type sx 1))))
	       (dolist (13nbor 13nbors)
		 (pushnew (list 1 13nbor sxid) subcmplx :test #'set-equal?)))
	     (let ((31nbors (append (neighbors-of-type sx 3))))
	       (dolist (31nbor 31nbors)
		 (pushnew (list 3 31nbor sxid) subcmplx :test #'set-equal?))))))
    subcmplx))

;; (1 | 2 3 4) (+) (1 5 | 3 4) --> (5 | 2 3 4) (+) (1 5 | 2 3) (+) (1 5 | 2 4)
(defun 2->3-move-internal-12 (13id 22id)
  "the 2,3 move performed on a (1,3) simplex attached to a (2,2) simplex"
  (let ((13sx nil) (22sx nil))
    (when (and (setf 13sx (get-3simplex 13id)) (setf 22sx (get-3simplex 22id)))
      (let* ((old2 (link-id 13id 22id))
	     (pts234 (3sx-hipts 13sx))    ;; hi points of the 1,3 simplex
	     (pts34 (3sx-hipts 22sx))  ;; hi points of the 2,2 simplex
	     (pts15 (3sx-lopts 22sx));; lo points of the 2,2 simplex
	     (pt5 (set-difference pts15 (3sx-lopts 13sx)))
	     (pt2 (set-difference pts234 pts34))
	     (new-internal-tlt-1 (list 1 (3sx-tmlo 13sx) (3sx-tmhi 13sx) 
				       (append pt5 pt2 (butlast pts34))))
	     (new-internal-tlt-2 (list 1 (3sx-tmlo 13sx) (3sx-tmhi 13sx)
				       (append pt5 pt2 (last pts34))))
	     (new-internal-tlt-3 (list 2 (3sx-tmlo 13sx) (3sx-tmhi 13sx) (append pts15 pt2)))
	     (new-13-pts (append pt5 pts234))
	     (new-22-pts-1 (append pts15 pt2 (butlast pts34)))
	     (new-22-pts-2 (append pts15 pt2 (last pts34)))
	     (nbors (set-difference (union (3sx-sx3ids 13sx) (3sx-sx3ids 22sx)) (list 0 13id 22id)))
	     (newsxdata nil))
	(unless (gethash (list (car pt2) (car pt5)) *1SIMPLEX?*)
	  (setf newsxdata (list (list 1 (3sx-tmlo 13sx) (3sx-tmhi 13sx) new-13-pts) 
				(list 2 (3sx-tmlo 22sx) (3sx-tmhi 22sx) new-22-pts-1) 
				(list 2 (3sx-tmlo 22sx) (3sx-tmhi 22sx) new-22-pts-2)))
	  (return-from 2->3-move-internal-12 (list newsxdata nbors (list 13id 22id) (list old2) DF23 nil)))))))

;; (2 3 4 | 1) (+) (3 4 | 1 5) --> (2 3 4 | 5) (+) (2 3 | 1 5) (+) (2 4 | 1 5)
(defun 2->3-move-internal-32 (31id 22id)
  "the 2,3 move performed on a (3,1) simplex attached to a (2,2) simplex"
    (let ((31sx nil) (22sx nil))
      (when (and (setf 31sx (get-3simplex 31id)) (setf 22sx (get-3simplex 22id)))
	(let* ((old2 (link-id 31id 22id))
	       (pts234 (3sx-lopts 31sx))    ;; lo points of the 3,1 simplex
	       (pts34 (3sx-lopts 22sx))  ;; lo points of the 2,2 simplex
	       (pts15 (3sx-hipts 22sx));; hi points of the 2,2 simplex
	       (pt5 (set-difference pts15 (3sx-hipts 31sx)))
	       (pt2 (set-difference pts234 pts34))
	       (new-internal-tlt-1 (list 2 (3sx-tmlo 31sx) (3sx-tmhi 31sx) 
					 (append pt2 (butlast pts34) pt5)))
	       (new-internal-tlt-2 (list 2 (3sx-tmlo 31sx) (3sx-tmhi 31sx)
					 (append pt2 (last pts34) pt5)))
	       (new-internal-tlt-3 (list 1 (3sx-tmlo 31sx) (3sx-tmhi 31sx) (append pt2 pts15)))
	       (new-31-pts (append pts234 pt5))
	       (new-22-pts-1 (append pt2 (butlast pts34) pts15))
	       (new-22-pts-2 (append pt2 (last pts34) pts15))
	       (nbors (set-difference (union (3sx-sx3ids 31sx) (3sx-sx3ids 22sx)) (list 0 31id 22id)))
	       (newsxdata nil))
	  (unless (gethash (list (car pt2) (car pt5)) *1SIMPLEX?*)
	    (setf newsxdata (list (list 3 (3sx-tmlo 31sx) (3sx-tmhi 31sx) new-31-pts) 
				  (list 2 (3sx-tmlo 22sx) (3sx-tmhi 22sx) new-22-pts-1) 
				  (list 2 (3sx-tmlo 22sx) (3sx-tmhi 22sx) new-22-pts-2)))
	    (return-from 2->3-move-internal-32 (list newsxdata nbors (list 31id 22id) (list old2) DF23 nil)))))))

(defun try-2->3 (sxid)
  (let ((subcmplx (2->3-subcomplex sxid))
	(movedata nil))
    (unless (null subcmplx)
      (dolist (curr subcmplx)
	(cond ((= 1 (first curr))
	       (setf movedata (2->3-move-internal-12 (second curr) (third curr)))
	       (when movedata
		 (return-from try-2->3 movedata)))
	      ((= 3 (first curr))
	       (setf movedata (2->3-move-internal-32 (second curr) (third curr)))
	       (when movedata
		 (return-from try-2->3 movedata))))))))

;;
;; Code that supports a 3->2 move
;;

(defun 3->2-subcomplex (sxid)
  "returns a list of the form ((1or3 13or31id 22id1 22id2)...) where the first number 1 or 3 tells us 
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
			   (pushnew (list (3sx-type sx) sxid 22nbor 22nborof22nbor) subcmplx 
				    :test #'set-equal?)))))))))
	    ((= 2 (3sx-type sx))
	     (let ((22nbors (neighbors-of-type sx 2)))
	       (dolist (22nbor 22nbors)
		 (let ((22sx (get-3simplex 22nbor)))
		   (when 22sx
		     (let ((13nborsof22nbor (neighbors-of-type 22sx 1)))
		       (dolist (13nborof22nbor 13nborsof22nbor)
			 (when (3simplices-connected? 13nborof22nbor sxid)
			   (pushnew (list 1 13nborof22nbor sxid 22nbor) subcmplx :test #'set-equal?))))
		     (let ((31nborsof22nbor (neighbors-of-type 22sx 3)))
		       (dolist (31nborof22nbor 31nborsof22nbor)
			 (when (3simplices-connected? 31nborof22nbor sxid)
			   (pushnew (list 3 31nborof22nbor sxid 22nbor) subcmplx :test #'set-equal?)))))))))))
    subcmplx))
;; (5 | 2 3 4) (+) (1 5 | 2 3) (+) (1 5 | 2 4) --> (1 | 2 3 4) (+) (1 5 | 3 4)
(defun 3->2-move-internal-122 (13id 22id1 22id2)
  "the (3,2) move performed on a (1,3) simplex attached to two (2,2) simplices"
  (let ((13sx nil) (22sx1 nil) (22sx2 nil))
    (when (and (setf 13sx (get-3simplex 13id)) 
	       (setf 22sx1 (get-3simplex 22id1))
	       (setf 22sx2 (get-3simplex 22id2)))
      (let* ((old2s (list (link-id 13id 22id1) (link-id 13id 22id2) (link-id 22id1 22id2)))
	     (pts234 (3sx-hipts 13sx))
	     (pts15 (3sx-lopts 22sx1))
	     (pt1 (set-difference pts15 (3sx-lopts 13sx)))
	     (pts34 (set-exclusive-or (3sx-hipts 22sx1) (3sx-hipts 22sx2)))
	     (new-internal-triangle (list 1 (3sx-tmlo 13sx) (3sx-tmhi 13sx) (append pt1 pts34)))
	     (new-13-pts (append pt1 pts234))
	     (new-22-pts (append pts15 pts34))
	     (nbors (set-difference (unions (3sx-sx3ids 13sx) (3sx-sx3ids 22sx1) (3sx-sx3ids 22sx2)) 
				    (list 0 13id 22id1 22id2)))
	     (newsxdata nil))
	(unless (gethash new-internal-triangle *2SIMPLEX->ID*)
	  (setf newsxdata (list (list 1 (3sx-tmlo 13sx) (3sx-tmhi 13sx) new-13-pts) 
				(list 2 (3sx-tmlo 22sx1) (3sx-tmhi 22sx2) new-22-pts)))
	  (return-from 3->2-move-internal-122 (list newsxdata nbors (list 13id 22id1 22id2) old2s DF32 nil)))))))
	
(defun 3->2-move-internal-322 (31id 22id1 22id2)
  "the (3,2) move performed on a (3,1) simplex attached to two (2,2) simplices"
  (let ((31sx nil) (22sx1 nil) (22sx2 nil))
    (when (and (setf 31sx (get-3simplex 31id)) 
	       (setf 22sx1 (get-3simplex 22id1))
	       (setf 22sx2 (get-3simplex 22id2)))
      (let* ((old2s (list (link-id 31id 22id1) (link-id 31id 22id2) (link-id 22id1 22id2)))
	     (pts234 (3sx-lopts 31sx)) ;; lo points of 3,1
	     (pts15 (3sx-hipts 22sx1))   ;; hi points of 2,2
	     (pt1 (set-difference pts15 (3sx-hipts 31sx))) ;; 2,2 hi - 3,1 hi
	     (pts34 (set-exclusive-or (3sx-lopts 22sx1) (3sx-lopts 22sx2)))
	     (new-internal-triangle (list 2 (3sx-tmlo 31sx) (3sx-tmhi 31sx) (append pts34 pt1)))
	     (new-31-pts (append pts234 pt1))
	     (new-22-pts (append pts34 pts15))
	     (nbors (set-difference (unions (3sx-sx3ids 31sx) (3sx-sx3ids 22sx1) (3sx-sx3ids 22sx2)) 
				    (list 0 31id 22id1 22id2)))
	     (newsxdata nil))
	(unless (gethash new-internal-triangle *2SIMPLEX->ID*)
	  (setf newsxdata (list (list 3 (3sx-tmlo 31sx) (3sx-tmhi 31sx) new-31-pts) 
				(list 2 (3sx-tmlo 22sx1) (3sx-tmhi 22sx2) new-22-pts))) 
	  (return-from 3->2-move-internal-322 (list newsxdata nbors (list 31id 22id1 22id2) old2s DF32 nil)))))))

(defun try-3->2 (sxid)
  (let ((subcmplx (3->2-subcomplex sxid))
	(movedata nil))
    (unless (null subcmplx)
      (dolist (curr subcmplx)
	(cond ((= 1 (first curr))
	       (setf movedata (3->2-move-internal-122 (second curr) (third curr) (fourth curr)))
	       (when movedata
		 (return-from try-3->2 movedata)))
	      ((= 3 (first curr))
	       (setf movedata (3->2-move-internal-322 (second curr) (third curr) (fourth curr)))
	       (when movedata 
		 (return-from try-3->2 movedata))))))))