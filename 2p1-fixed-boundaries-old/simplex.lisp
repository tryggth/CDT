;; 2simplex
;; (type tmlo tmhi (p0 p1 p2))
;; type = 0,1,2,3 tm[lo|hi] = [lo|hi] spatial time slice pj = points, nj = neighbor that does not have pj

;; spatial-2simplex
;; (time (p0 p1 p2) (n0 n1 n2))

;; 3simplex
;; (type tmlo tmhi (n0 n1 n2 n3) (t0 t1 t2 t3))
;; type = 1,2,3 tm[lo|hi] - [lo|hi] spatial time slice
;; tj = id of the 2sx that does not have pj
;; nj = id of the 3sx that does not have pj

;;macro to decide whether or not to apply a modulo division, depending on
;;whether periodic BC's are in use
(defmacro bc-mod (num)
  `(if (string= BCTYPE "PERIODIC")
       (mod ,num NUM-T)
       ,num))

(defun make-2simplex (type tmlo tmhi p0 p1 p2)
  "makes and returns the id of the 2-simplex with the specified data. If a 2-simplex with points 
p0 p1 p2 already exists, the id of that simplex is returned"
  (let* ((2sx (list type tmlo tmhi (list p0 p1 p2)))
	 (2sxid (gethash 2sx *2SIMPLEX->ID*)))
    (unless 2sxid
      (setf 2sxid (next-2simplex-id))
      (setf (gethash 2sx *2SIMPLEX->ID*) 2sxid)
      (setf (gethash 2sxid *ID->2SIMPLEX*) 2sx))
    2sxid))

(defmacro get-2simplex (sxid) `(gethash ,sxid *ID->2SIMPLEX*))

(defmacro 2sx-type (2sx) `(first ,2sx))
(defmacro 2sx-tmlo (2sx) `(second ,2sx))
(defmacro 2sx-tmhi (2sx) `(third ,2sx)) 
(defmacro 2sx-points (2sx) `(fourth ,2sx))

(defmacro s2sx-time (s2sx) `(first ,s2sx))
(defmacro s2sx-points (s2sx) `(second ,s2sx))
(defmacro s2sx-sx2ids (s2sx) `(third ,s2sx))

(defun make-s2simplex (triangle-time triangle-points)
  (let ((stid (next-spatial-2simplex-id)))
    (setf (gethash stid *ID->SPATIAL-2SIMPLEX*) 
	  (list triangle-time (copy-list triangle-points) (list 0 0 0)))
    stid))
;; this version is used for loading the simplex data from file
(defun make-s2simplex-v2 (s2sxid s2sxtm p0 p1 p2 n0 n1 n2)
  (setf (gethash s2sxid *ID->SPATIAL-2SIMPLEX*)
	(list s2sxtm (list p0 p1 p2) (list n0 n1 n2))))
(defmacro get-s2simplex (stid) `(gethash ,stid *ID->SPATIAL-2SIMPLEX*))

(defmacro remove-2simplex (2sxid)
  `(let ((2sx (gethash ,2sxid *ID->2SIMPLEX*)))
     (remhash ,2sxid *ID->2SIMPLEX*)
     (remhash 2sx *2SIMPLEX->ID*)))

(defun remove-2simplices (2sxids)
  (dolist (2sxid 2sxids)
    (let ((2sx (gethash 2sxid *ID->2SIMPLEX*)))
      (remhash 2sxid *ID->2SIMPLEX*)
      (remhash 2sx *2SIMPLEX->ID*))))

(defun show-2simplex->id-store ()
  (maphash #'(lambda (2sx 2sxid) (format t "~A [~A]~%" 2sx 2sxid)) *2SIMPLEX->ID*))

(defun show-id->2simplex-store ()
  (maphash #'(lambda (2sxid 2sx) (format t "[~A] ~A~%" 2sxid 2sx)) *ID->2SIMPLEX*))

(defun connect-spatial-2simplices (st1id st2id)
  (let ((st1 nil) (st2 nil))
    (when (and (setf st1 (get-s2simplex st1id)) (setf st2 (get-s2simplex st2id)))
      (let* ((points1 (s2sx-points st1))
	     (points2 (s2sx-points st2))
	     (line (intersection points1 points2)))
	(when (= 2 (length line))
	  (let ((pos1 (position (first (set-difference points1 line)) points1))
		(pos2 (position (first (set-difference points2 line)) points2)))
	    (setf (nth pos1 (s2sx-sx2ids st1)) st2id)(setf (nth pos2 (s2sx-sx2ids st2)) st1id)))))))

(defun connect-spatial-2simplices-within-list (sx1ids)
  (for (n 0 (1- (length sx1ids)))
       (for (m (1+ n) (1- (length sx1ids)))
	    (connect-spatial-2simplices (nth n sx1ids) (nth m sx1ids)))))

(defun make-3simplex (type tmlo tmhi p0 p1 p2 p3)
  (let ((t0type nil) (t1type nil) (t2type nil) (t3type nil) (sx3id (next-3simplex-id)))
    (ecase type
      (1 (setf t0type 0) (setf t1type 1) (setf t2type 1)(setf t3type 1))
      (2 (setf t0type 1) (setf t1type 1) (setf t2type 2)(setf t3type 2))
      (3 (setf t0type 2) (setf t1type 2) (setf t2type 2)(setf t3type 3)))
    (setf (gethash sx3id *ID->3SIMPLEX*)
	  (list type tmlo tmhi
		(list p0 p1 p2 p3)
		(list 0 0 0 0)
		(list (make-2simplex t0type tmlo tmhi p1 p2 p3)
		      (make-2simplex t1type tmlo tmhi p0 p2 p3)
		      (make-2simplex t2type tmlo tmhi p0 p1 p3)
		      (make-2simplex t3type tmlo tmhi p0 p1 p2))))
    sx3id))

;; same as above except the points are "packed" into a list; useful during the moves, when the points
;; of the simplex are computed via append / unions / intersections in the form of a list
(defun make-3simplex-v2 (type tmlo tmhi pts)
  (let ((t0type nil) (t1type nil) (t2type nil) (t3type nil) (sx3id (next-3simplex-id)))
    (ecase type
      (1 (setf t0type 0) (setf t1type 1) (setf t2type 1)(setf t3type 1))
      (2 (setf t0type 1) (setf t1type 1) (setf t2type 2)(setf t3type 2))
      (3 (setf t0type 2) (setf t1type 2) (setf t2type 2)(setf t3type 3)))
    (setf (gethash sx3id *ID->3SIMPLEX*)
	  (list type tmlo tmhi
		(copy-list pts)
		(list 0 0 0 0)
		(list (make-2simplex t0type tmlo tmhi (nth 1 pts) (nth 2 pts) (nth 3 pts))
		      (make-2simplex t1type tmlo tmhi (nth 0 pts) (nth 2 pts) (nth 3 pts))
		      (make-2simplex t2type tmlo tmhi (nth 0 pts) (nth 1 pts) (nth 3 pts))
		      (make-2simplex t3type tmlo tmhi (nth 0 pts) (nth 1 pts) (nth 2 pts)))))
    sx3id))

;; this version is used only during initialization. If periodic b.c. are specified, it adjusts the
;; points on the final time slice, since the t=T slice is identified with t=0 slice.
(defun make-3simplex-v3 (type tmlo tmhitmp p0tmp p1tmp p2tmp p3tmp)
  (let ((t0type nil) (t1type nil) (t2type nil) (t3type nil) (sx3id (next-3simplex-id))
	(p0 p0tmp) (p1 p1tmp) (p2 p2tmp) (p3 p3tmp) (tmhi tmhitmp))
    (when (and (string= BCTYPE "PERIODIC") (= NUM-T tmhitmp))
      (setf tmhi 0)
      (cond ((= 1 type)
	     (decf p1 (* N0-PER-SLICE NUM-T)) 
	     (decf p2 (* N0-PER-SLICE NUM-T)) 
	     (decf p3 (* N0-PER-SLICE NUM-T)))
	    ((= 2 type)
	     (decf p2 (* N0-PER-SLICE NUM-T)) 
	     (decf p3 (* N0-PER-SLICE NUM-T)))
	    ((= 3 type)
	     (decf p3 (* N0-PER-SLICE NUM-T)))))
    (ecase type
      (1 (setf t0type 0) (setf t1type 1) (setf t2type 1)(setf t3type 1))
      (2 (setf t0type 1) (setf t1type 1) (setf t2type 2)(setf t3type 2))
      (3 (setf t0type 2) (setf t1type 2) (setf t2type 2)(setf t3type 3)))

    ;;set last used point if a new point is used
    (setf *LAST-USED-POINT* (max *LAST-USED-POINT* p0 p1 p2 p3))

    (setf (gethash sx3id *ID->3SIMPLEX*)
	  (list type tmlo tmhi
		(list p0 p1 p2 p3)
		(list 0 0 0 0)
		(list (make-2simplex t0type tmlo tmhi p1 p2 p3)
		      (make-2simplex t1type tmlo tmhi p0 p2 p3)
		      (make-2simplex t2type tmlo tmhi p0 p1 p3)
		      (make-2simplex t3type tmlo tmhi p0 p1 p2))))))

;; this version is used for loading the simplex data from file
(defun make-3simplex-v4 (type tmlo tmhi p0 p1 p2 p3 n0 n1 n2 n3 sx3id)
  (let ((t0type nil) (t1type nil) (t2type nil) (t3type nil))
    (ecase type
      (1 (setf t0type 0) (setf t1type 1) (setf t2type 1)(setf t3type 1))
      (2 (setf t0type 1) (setf t1type 1) (setf t2type 2)(setf t3type 2))
      (3 (setf t0type 2) (setf t1type 2) (setf t2type 2)(setf t3type 3)))
    (setf (gethash sx3id *ID->3SIMPLEX*)
	  (list type tmlo tmhi
		(list p0 p1 p2 p3)
		(list n0 n1 n2 n3)
		(list (make-2simplex t0type tmlo tmhi p1 p2 p3)
		      (make-2simplex t1type tmlo tmhi p0 p2 p3)
		      (make-2simplex t2type tmlo tmhi p0 p1 p3)
		      (make-2simplex t3type tmlo tmhi p0 p1 p2))))))

;; a replacement for v2
(defun make-3simplex-v5 (simplex-data)
  (let ((type (first simplex-data))
	(tmlo (second simplex-data))
	(tmhi (third simplex-data))
	(pts (fourth simplex-data))
	(t0type nil) (t1type nil) (t2type nil) (t3type nil) (sx3id (next-3simplex-id)))
    (ecase type
      (1 (setf t0type 0) (setf t1type 1) (setf t2type 1)(setf t3type 1))
      (2 (setf t0type 1) (setf t1type 1) (setf t2type 2)(setf t3type 2))
      (3 (setf t0type 2) (setf t1type 2) (setf t2type 2)(setf t3type 3)))
    (setf (gethash sx3id *ID->3SIMPLEX*)
	  (list type tmlo tmhi
		(copy-list pts)
		(list 0 0 0 0)
		(list (make-2simplex t0type tmlo tmhi (nth 1 pts) (nth 2 pts) (nth 3 pts))
		      (make-2simplex t1type tmlo tmhi (nth 0 pts) (nth 2 pts) (nth 3 pts))
		      (make-2simplex t2type tmlo tmhi (nth 0 pts) (nth 1 pts) (nth 3 pts))
		      (make-2simplex t3type tmlo tmhi (nth 0 pts) (nth 1 pts) (nth 2 pts)))))
    sx3id))

;; simplex-data-list = ((typ tmlo tmhi (p0 p1 p2 p3)) (typ tmlo tmhi (p0 p1 p2 p3))...)
;; the ids of the simplices are returned
(defun make-3simplices-in-bulk (simplex-data-list)
  (let ((3sxids nil))
    (dolist (simplex-data simplex-data-list)
      (push (make-3simplex-v5 simplex-data) 3sxids))
    3sxids))

(defmacro 3sx-type (sx) `(nth 0 ,sx))
(defmacro 3sx-tmlo (sx) `(nth 1 ,sx))
(defmacro 3sx-tmhi (sx) `(nth 2 ,sx))
(defmacro 3sx-points (sx) `(nth 3 ,sx))
(defmacro 3sx-sx3ids (sx) `(nth 4 ,sx))
(defmacro 3sx-sx2ids (sx) `(nth 5 ,sx))
(defmacro 3sx-lopts (sx) `(subseq (3sx-points ,sx) 0 (3sx-type ,sx)))
(defmacro 3sx-hipts (sx) `(subseq (3sx-points ,sx) (3sx-type ,sx)))

(defmacro get-3simplex (sxid)
  `(gethash ,sxid *ID->3SIMPLEX*))
(defmacro remove-3simplex (3sxid)
  `(progn
     (remhash ,3sxid *ID->3SIMPLEX*)
     (recycle-3simplex-id ,3sxid)))
(defmacro remove-3simplices (3sxids)
  `(dolist (3sxid ,3sxids)
     (remhash 3sxid *ID->3SIMPLEX*)
     (recycle-3simplex-id 3sxid)))
(defun show-id->3simplex-store ()
  (maphash #'(lambda (3sxid 3sx) (format t "[~A] ~A~%" 3sxid 3sx)) *ID->3SIMPLEX*))

(defmacro nth-point (sx n)
  `(nth ,n (3sx-points ,sx)))

(defun connect-3simplices (sx1id sx2id)
  (let ((sx1 nil) (sx2 nil))
    (when (and (setf sx1 (get-3simplex sx1id)) (setf sx2 (get-3simplex sx2id)))
      (let ((2sxlinkid (intersection (3sx-sx2ids sx1) (3sx-sx2ids sx2))))
	(when (= 1 (length 2sxlinkid))
	  (let ((pos1 (position (first 2sxlinkid) (3sx-sx2ids sx1)))
		(pos2 (position (first 2sxlinkid) (3sx-sx2ids sx2))))
	    (setf (nth pos1 (3sx-sx3ids sx1)) sx2id) (setf (nth pos2 (3sx-sx3ids sx2)) sx1id)))))))

(defun connect-3simplices-within-list (sx1ids)
  (for (n 0 (1- (length sx1ids)))
       (for (m (1+ n) (1- (length sx1ids)))
	    (connect-3simplices (nth n sx1ids) (nth m sx1ids)))))

(defun connect-3simplices-across-lists (sx1ids sx2ids)
  (dolist (sx1id sx1ids)
    (dolist (sx2id sx2ids)
      (connect-3simplices sx1id sx2id))))

(defmacro 3simplices-connected? (sxid1 sxid2)
  `(let ((sx1 nil) (sx2 nil))
     (and (setf sx1 (get-3simplex ,sxid1)) (setf sx2 (get-3simplex ,sxid2))
	  (find ,sxid1 (3sx-sx3ids sx2)) (find ,sxid2 (3sx-sx3ids sx1)))))

(defun get-simplices-in-sandwich (tlo thi)
  (let ((sxids '()))
    (maphash #'(lambda (id sx)
		 (when (and (= (3sx-tmlo sx) (bc-mod tlo)) 
			    (= (3sx-tmhi sx) (bc-mod thi)))
		   (push id sxids)))
	     *ID->3SIMPLEX*)
    sxids))

(defun get-simplices-in-sandwich-of-type (tlo thi typ)
  (let ((sxids '()))
    (maphash #'(lambda (id sx)
		 (when (and (= (3sx-tmlo sx) (bc-mod tlo)) 
			    (= (3sx-tmhi sx) (bc-mod thi))
			    (= (3sx-type sx) typ))
		   (push id sxids)))
	     *ID->3SIMPLEX*)
    sxids))

(defun get-2simplices-in-sandwich-of-type (tlo thi typ)
  (let ((sxids '()))
    (maphash #'(lambda (id sx)
		 (when (and (= (2sx-tmlo sx) (bc-mod tlo)) 
			    (= (2sx-tmhi sx) (bc-mod thi))
			    (= (2sx-type sx) typ))
		   (push id sxids)))
	     *ID->2SIMPLEX*)
    sxids))

;; returns (1ids 2ids 3ids) where 1ids is a list of (1,3) ids in the sandwich etc.
(defun get-simplices-in-sandwich-ordered-by-type (tlo thi)
  (let ((1ids nil) (2ids nil) (3ids nil))
    (maphash #'(lambda (id sx)
		 (when (and (= (3sx-tmlo sx) (bc-mod tlo)) 
			    (= (3sx-tmhi sx) (bc-mod thi)))
		   (ecase (3sx-type sx)
		     (1 (push id 1ids))
		     (2 (push id 2ids))
		     (3 (push id 3ids)))))
	     *ID->3SIMPLEX*)
    (values 1ids 2ids 3ids)))

;;define macro to help weed out fictitious simplices added for open boundaries
(defmacro is-real-simplex (sx)
  `(and (>= (3sx-tmlo ,sx) 0) (<= (3sx-tmhi ,sx) NUM-T)))

(defun get-simplices-of-type (typ)
  (let ((sxids '()))
    (maphash #'(lambda (id sx)
		 (if (and (is-real-simplex sx) (= (3sx-type sx) typ))
		     (push id sxids)))
	     *ID->3SIMPLEX*)
    sxids))

;;function to grab even the fake simplices beyond the boundaries
(defun get-real-and-fake-simplices-of-type (typ)
  (let ((sxids '()))
    (maphash #'(lambda (id sx)
		 (if (= (3sx-type sx) typ)
		     (push id sxids)))
	     *ID->3SIMPLEX*)
    sxids))

(defun count-simplices-of-type (typ)
  (let ((count 0))
    (maphash #'(lambda (id sx)
		 (declare (ignore id))
		 (if (and (is-real-simplex sx) (= (3sx-type sx) typ))
		     (incf count)))
	     *ID->3SIMPLEX*)
    count))

(defun count-simplices-in-sandwich (tlo thi)
  (let ((count 0))
    (maphash #'(lambda (id sx)
		 (declare (ignore id))
		 (when (and (= (3sx-tmlo sx) (bc-mod tlo)) 
			    (= (3sx-tmhi sx) (bc-mod thi)))
		     (incf count)))
	     *ID->3SIMPLEX*)
    count))

(defun count-simplices-of-all-types ()
  (let ((1count 0) (2count 0) (3count 0))
    (maphash #'(lambda (id sx)
		 (declare (ignore id))
		 (when (is-real-simplex sx)
		   (ecase (3sx-type sx)
		     (1 (incf 1count))
		     (2 (incf 2count))
		     (3 (incf 3count)))))
	     *ID->3SIMPLEX*)
    (list 1count 2count 3count (+ 1count 2count 3count))))

(defun count-boundary-vs-bulk ()
  (let ((first-13 0) (last-13 0) (bulk-13 0)
	(first-22 0) (last-22 0) (bulk-22 0)
	(first-31 0) (last-31 0) (bulk-31 0))
    (maphash #'(lambda (id sx)
		 (declare (ignore id))
		 (let* ((ty (3sx-type sx))
			(tl (3sx-tmlo sx))
			(th (3sx-tmhi sx))
			(in-first (= tl 0))
			(in-last  (= th NUM-T))
			(outside  (or (= tl -1) (= tl NUM-T)))
			(in-bulk  (not (or in-first in-last outside))))
		   (cond
		     ((and in-first (= ty 1)) (incf first-13))
		     ((and in-first (= ty 2)) (incf first-22))
		     ((and in-first (= ty 3)) (incf first-31))
		     ((and in-bulk  (= ty 1)) (incf bulk-13))
		     ((and in-bulk  (= ty 2)) (incf bulk-22))
		     ((and in-bulk  (= ty 3)) (incf bulk-31))
		     ((and in-last  (= ty 1)) (incf last-13))
		     ((and in-last  (= ty 2)) (incf last-22))
		     ((and in-last  (= ty 3)) (incf last-31)))))
	     *ID->3SIMPLEX*)
    (list first-13 first-22 first-31
	  bulk-13  bulk-22  bulk-31
	  last-13  last-22  last-31
	  (+ first-13 first-22 first-31
	     bulk-13 bulk-22 bulk-31
	     last-13 last-22 last-31))))

(defun count-simplices-in-sandwich-of-type (tlo thi typ)
  (let ((count 0))
    (maphash #'(lambda (id sx)
		 (declare (ignore id))
		 (when (and (= (3sx-tmlo sx) (bc-mod tlo)) 
			    (= (3sx-tmhi sx) (bc-mod thi))
			    (= (3sx-type sx) typ))
		   (incf count)))
	     *ID->3SIMPLEX*)
    count))

;;function to count the number of points in a particular time slice
(defun count-points-at-time (t0)

  ;;for all but the first time slice, we can get the upper points from all
  ;;simplicies in the (1- t0)/t0 time sandwich.  otherwise, we need to get 
  ;;the lower points of the first sandwich.

  (let ((list-of-points-with-duplicates () )
	(previous-element 0)
	(count 0))
    (if (= t0 0)
	(dolist (sxid (get-simplices-in-sandwich 0 1))
	  (dolist (i (3sx-lopts (get-3simplex sxid)))
	    (push i list-of-points-with-duplicates)))
	(dolist (sxid (get-simplices-in-sandwich (1- t0) t0))
	  (dolist (i (3sx-hipts (get-3simplex sxid)))
	    (push i list-of-points-with-duplicates))))

    ;;sorting the list allows me to count the number of unique points faster
    (dolist (i (sort list-of-points-with-duplicates #'<))
      (when (not (= previous-element i))
	(incf count)
	(setf previous-element i)))
    count))

;;function to count the number of timelike links in a sandwich
(defun count-timelike-links-in-sandwich (t0 t1)
 
  ;;approach: get all of the simplices in the sandwich, then repeatedly
  ;;perform set unions with the links from each simplex
  
  (let ((list-of-links () ))
    (dolist (sxid (get-simplices-in-sandwich t0 t1))
      (let* ((sx (get-3simplex sxid))
	     (lopts (3sx-lopts sx))
	     (hipts (3sx-hipts sx))
	     (sx-tl-links (let ((retval ()))
			    (dolist (lopt lopts)
			      (dolist (hipt hipts)
				(push (list lopt hipt) retval)))
			    retval)))
	(setf list-of-links (union list-of-links sx-tl-links 
				   :test #'(lambda (x y) (not (set-difference x y)))))))
    (length list-of-links)))

;;count the number of spacelike triangles at a certain time
(defun count-spacelike-triangles-at-time (t0)

  ;;must treat one end as a special case; choosing to treat
  ;;time 0 specially. number of triangles is just related to the
  ;;number of type 1/3 simplices in nieghboring sandwiches

  (if (= t0 0)
      (count-simplices-in-sandwich-of-type 0 1 3)
      (count-simplices-in-sandwich-of-type (1- t0) t0 1)))

;;count the number of spacelike links at a particular time, using the euler characteristic
;;and the numbers of triangles and points
(defun count-spacelike-links-at-time (t0)

  (let* ((chi (euler-char))
	 (n0  (count-points-at-time t0))
	 (n2  (count-spacelike-triangles-at-time t0)))

    ;;  (chi = n0 - n1 + n2)  =>  (n1 = n0 + n2 - chi)
    (- (+ n0 n2) chi)))
				    
;;count the number of timelike triangles in a particular sandwich
(defun count-timelike-triangles-in-sandwich (t0 t1)
  
  ;;similar approach to counting timelike links, except sets of three points
  ;;are collected and counted

  (let ((list-of-triangles () ))
    (dolist (sxid (get-simplices-in-sandwich t0 t1))
      (let* ((sx (get-3simplex sxid))
	     (ty (3sx-type sx))
	     (pts (3sx-points sx))

	     (sx-tl-triangles
	      (case ty
		(1 (list (list (first pts)  (second pts) (third pts))
			 (list (first pts)  (second pts) (fourth pts))
			 (list (first pts)  (third pts)  (fourth pts))))
		(2 (list (list (first pts)  (third pts)  (fourth pts))
			 (list (second pts) (third pts)  (fourth pts))
			 (list (first pts)  (second pts) (third pts))
			 (list (first pts)  (second pts) (fourth pts))))
		(3 (list (list (fourth pts) (first pts)  (second pts))
			 (list (fourth pts) (first pts)  (third pts))
			 (list (fourth pts) (second pts) (third pts)))))))

	(setf list-of-triangles (union list-of-triangles sx-tl-triangles
				       :test #'(lambda (x y) (not (set-difference x y)))))))
    (length list-of-triangles)))


;;(defun connect-simplices-in-sandwich (tlo thi)
;;(connect-3simplices-within-list (get-simplices-in-sandwich tlo thi)))

;; in a given sandwich
;; (1,3) can be connected to a (2,2) and cannot be connected to a (3,1)
;; a (2,2) can be connected to a (1,3) and a (3,1)
;; a (3,2) can be connected to a (2,2) and cannot be connected to a (1,3)
(defun connect-simplices-in-sandwich (tlo thi)
  (multiple-value-bind (1ids 2ids 3ids) (get-simplices-in-sandwich-ordered-by-type tlo thi)
    (connect-3simplices-across-lists 1ids 2ids)
    (connect-3simplices-across-lists 2ids 3ids)
    (connect-3simplices-within-list 1ids)
    (connect-3simplices-within-list 2ids)
    (connect-3simplices-within-list 3ids)))

(defun connect-simplices-in-adjacent-sandwiches (tl tm th)
  (connect-3simplices-across-lists 
   (get-simplices-in-sandwich-of-type tl tm 1)
   (get-simplices-in-sandwich-of-type tm th 3)))

(defun check-13-and-31 (tlo thi)
  (let ((13ids (get-simplices-in-sandwich-of-type tlo thi 1))
	(31ids (get-simplices-in-sandwich-of-type tlo thi 3))
	(problem-ids '()))
    (dolist (s13 13ids)
      (dolist (d13 13ids)
	(when (and (/= s13 d13) (set-equal? (subseq (3sx-points (get-3simplex s13)) 1)
					    (subseq (3sx-points (get-3simplex d13)) 1)))
	  (push (list s13 d13) problem-ids))))
    (dolist (s31 31ids)
      (dolist (d31 31ids)
	(when (and (/= s31 d31) (set-equal? (subseq (3sx-points (get-3simplex s31)) 0 3)
					    (subseq (3sx-points (get-3simplex d31)) 0 3)))
	  (push (list s31 d31) problem-ids))))
    problem-ids))

(defun check-all-slices-for-problem-simplices ()
  (for (ts 0 (- NUM-T 1))
       (format t "slice ~A has ~A problem simplices~%" ts (check-13-and-31 ts (1+ ts)))
       (finish-output)))

(defun check-all-slices-for-simplices-with-missing-neighbors ()
  (let ((problem-ids '()))
    (maphash #'(lambda (id sx)
		 (for (n 0 3)
		   (if (= 0 (nth n (3sx-sx3ids sx)))
		       (push id problem-ids))))
	     *ID->3SIMPLEX*)
    problem-ids))

;; if the 3-simplices are connected, returns the id of the linking 2-simplex else returns 0
(defmacro link-id (sxid1 sxid2)
  `(let ((sx1 nil) (sx2 nil) (link nil))
     (if (and (setf sx1 (get-3simplex ,sxid1)) 
	      (setf sx2 (get-3simplex ,sxid2))
	      (setf link (intersection (3sx-sx2ids sx2) (3sx-sx2ids sx1))))
	 (first link)
	 0)))

(defun neighbors-of-type (sx type)
  (let ((nbors nil)
	(nsx nil)
	(nids (3sx-sx3ids sx)))
    (for (n 0 3)
      (when (and (setf nsx (get-3simplex (nth n nids))) (= type (3sx-type nsx)))
	(pushnew (nth n nids) nbors)))
    nbors))

(defun save-spacetime-to-file (outfile)
  (format outfile "~A ~A ~A ~A ~A ~A ~A ~A ~A ~A ~A ~A ~A ~A ~A ~A ~A~%" 
	  BCTYPE STOPOLOGY NUM-T N-INIT *LAST-USED-POINT* *LAST-USED-3SXID* 
	  N0 N1-SL N1-TL N2-SL N2-TL N3-TL-31 N3-TL-22 eps k0 k3 *alpha*)
  (maphash #'(lambda (k v)
	       (format outfile "~A ~A ~A ~A ~A ~A ~A ~A ~A ~A ~A ~A~%" 
		       (3sx-type v) (3sx-tmlo v) (3sx-tmhi v)
		       (nth-point v 0) (nth-point v 1) (nth-point v 2) (nth-point v 3)
		       (nth 0 (3sx-sx3ids v)) (nth 1 (3sx-sx3ids v)) 
		       (nth 2 (3sx-sx3ids v)) (nth 3 (3sx-sx3ids v)) k))
	   *ID->3SIMPLEX*))

(defun save-s2simplex-data-to-file (outfile)
  (maphash #'(lambda (k v)
	       (let ((pts (s2sx-points v))
		     (nbors (s2sx-sx2ids v)))
		 (format outfile "~A ~A ~A ~A ~A ~A ~A ~A~%"
			 k (s2sx-time v) (nth 0 pts) (nth 1 pts) (nth 2 pts) (nth 0 nbors) (nth 1 nbors)
			 (nth 2 nbors))))
	   *ID->SPATIAL-2SIMPLEX*))

(defun parse-parameters-line (line)
  (with-input-from-string (s line)
    (let ((data (loop
		   :for num := (read s nil nil)
		   :while num
		   :collect num)))
      (setf BCTYPE (nth 0 data)) (setf STOPOLOGY (nth 1 data))
      (setf NUM-T (nth 2 data)) (setf N-INIT (nth 3 data))
      (setf *LAST-USED-POINT* (nth 4 data)) (setf *LAST-USED-3SXID* (nth 5 data))
      (setf N0 (nth 6 data)) (setf N1-SL (nth 7 data)) (setf N1-TL (nth 8 data))
      (setf N2-SL (nth 9 data)) (setf N2-TL (nth 10 data)) (setf N3-TL-31 (nth 11 data))
      (setf N3-TL-22 (nth 12 data)) (setf eps (nth 13 data)) 
      (set-k0-k3-alpha (nth 14 data) (nth 15 data) (nth 16 data)))))

(defun parse-simplex-data-line (line)
  (with-input-from-string (s line)
    (let ((data (loop
		   :for num := (read s nil nil)
		   :while num
		   :collect num)))
      (make-3simplex-v4 (nth 0 data) (nth 1 data) (nth 2 data) (nth 3 data) 
			(nth 4 data) (nth 5 data) (nth 6 data) (nth 7 data)
			(nth 8 data) (nth 9 data) (nth 10 data) (nth 11 data)))))

(defun load-spacetime-from-file (infile)
  (parse-parameters-line (read-line infile nil))
  (loop for line = (read-line infile nil)
     while line do (parse-simplex-data-line line)))

(defun parse-s2simplex-data-line (line)
  (with-input-from-string (s line)
    (let ((data (loop
		   :for num := (read s nil nil)
		   :while num
		   :collect num)))
      (make-s2simplex-v2 (nth 0 data) (nth 1 data) (nth 2 data) (nth 3 data) (nth 4 data)
			 (nth 5 data) (nth 6 data) (nth 7 data)))))

(defun load-s2simplex-data-from-file (infile)
  (clrhash *ID->SPATIAL-2SIMPLEX*)
  (loop for line = (read-line infile nil)
     while line do (parse-s2simplex-data-line line)))
