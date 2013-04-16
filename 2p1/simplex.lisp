(declaim (optimize (speed 3)
		   (compilation-speed 0)
		   (debug 0)
		   (safety 0)))
;; sl 1-simplex is (tslice (p0 p1))
;; tl 1-simplex is (type tmlo (p0 p1)) where type = 1
;; sl 2-simplex is (tslice (p0 p1 p2))
;; tl 2-simplex is (type tmlo (p0 p1 p2)) where type = 1,2
;; tl 3-simplex is (type tmlo tmhi (p0 p1 p2 p3) (n0 n1 n2 n3))
;; where type = 1,2,3,4 
;; nj = id of the 4sx that does not have pj

(defun make-2simplices (sxtype tmlo tmhi p0 p1 p2 p3)
  "makes all 2-simplices with the specified data iff they don't already exist"
  (cond ((= 1 sxtype) ; (p0 | p1 p2 p3)
	 (setf (gethash `(1 ,tmlo (,p0 ,p1 ,p2)) *TL2SIMPLEX->ID*) 0)
	 (setf (gethash `(1 ,tmlo (,p0 ,p1 ,p3)) *TL2SIMPLEX->ID*) 0)
	 (setf (gethash `(1 ,tmlo (,p0 ,p2 ,p3)) *TL2SIMPLEX->ID*) 0)
	 (setf (gethash `(,tmhi (,p1 ,p2 ,p3)) *SL2SIMPLEX->ID*) 0))
	((= 2 sxtype) ; (p0 p1 | p2 p3)
	 (setf (gethash `(2 ,tmlo (,p0 ,p1 ,p2)) *TL2SIMPLEX->ID*) 0)
	 (setf (gethash `(2 ,tmlo (,p0 ,p1 ,p3)) *TL2SIMPLEX->ID*) 0)
	 (setf (gethash `(1 ,tmlo (,p0 ,p2 ,p3)) *TL2SIMPLEX->ID*) 0)
	 (setf (gethash `(1 ,tmlo (,p1 ,p2 ,p3)) *TL2SIMPLEX->ID*) 0))
	((= 3 sxtype) ; (p0 p1 p2 | p3)
	 (setf (gethash `(,tmlo (,p0 ,p1 ,p2)) *SL2SIMPLEX->ID*) 0)
	 (setf (gethash `(2 ,tmlo (,p0 ,p1 ,p3)) *TL2SIMPLEX->ID*) 0)
	 (setf (gethash `(2 ,tmlo (,p0 ,p2 ,p3)) *TL2SIMPLEX->ID*) 0)
	 (setf (gethash `(2 ,tmlo (,p1 ,p2 ,p3)) *TL2SIMPLEX->ID*) 0))))

(defun make-1simplices (sxtype tmlo tmhi p0 p1 p2 p3)
  "makes all 1-simplices with the specified data iff they don't already exist"
  (cond ((= 1 sxtype) ;; p0 | p1 p2 p3
	 (setf (gethash `(1 ,tmlo (,p0 ,p1)) *TL1SIMPLEX->ID*) 0)
	 (setf (gethash `(1 ,tmlo (,p0 ,p2)) *TL1SIMPLEX->ID*) 0)
	 (setf (gethash `(1 ,tmlo (,p0 ,p3)) *TL1SIMPLEX->ID*) 0)
	 (setf (gethash `(,tmhi (,p1 ,p2)) *SL1SIMPLEX->ID*) 0)
	 (setf (gethash `(,tmhi (,p1 ,p3)) *SL1SIMPLEX->ID*) 0)
	 (setf (gethash `(,tmhi (,p2 ,p3)) *SL1SIMPLEX->ID*) 0))
	((= 2 sxtype) ;; p0 p1 | p2 p3
	 (setf (gethash `(,tmlo (,p0 ,p1)) *SL1SIMPLEX->ID*) 0)
	 (setf (gethash `(1 ,tmlo (,p0 ,p2)) *TL1SIMPLEX->ID*) 0)
	 (setf (gethash `(1 ,tmlo (,p0 ,p3)) *TL1SIMPLEX->ID*) 0)
	 (setf (gethash `(1 ,tmlo (,p1 ,p2)) *TL1SIMPLEX->ID*) 0)
	 (setf (gethash `(1 ,tmlo (,p1 ,p3)) *TL1SIMPLEX->ID*) 0)
	 (setf (gethash `(,tmhi (,p2 ,p3)) *SL1SIMPLEX->ID*) 0))
	((= 3 sxtype) ;; p0 p1 p2 | p3
	 (setf (gethash `(,tmlo (,p0 ,p1)) *SL1SIMPLEX->ID*) 0)
	 (setf (gethash `(,tmlo (,p0 ,p2)) *SL1SIMPLEX->ID*) 0)
	 (setf (gethash `(1 ,tmlo (,p0 ,p3)) *TL1SIMPLEX->ID*) 0)
	 (setf (gethash `(,tmlo (,p1 ,p2)) *SL1SIMPLEX->ID*) 0)
	 (setf (gethash `(1 ,tmlo (,p1 ,p3)) *TL1SIMPLEX->ID*) 0)
	 (setf (gethash `(1 ,tmlo (,p2 ,p3)) *TL1SIMPLEX->ID*) 0))))

(defun remove-tl2simplex (tl2sx)
  (remhash tl2sx *TL2SIMPLEX->ID*))
(defun remove-sl2simplex (sl2sx)
  (remhash sl2sx *SL2SIMPLEX->ID*))
(defun remove-tl2simplices (tl2sxs)
  (dolist (tl2sx tl2sxs)
    (remove-tl2simplex tl2sx)))
(defun remove-sl2simplices (sl2sxs)
  (dolist (sl2sx sl2sxs)
    (remove-sl2simplex sl2sx)))

(defun remove-tl1simplex (tl1sx)
  (remhash tl1sx *TL1SIMPLEX->ID*))
(defun remove-sl1simplex (sl1sx)
  (remhash sl1sx *SL1SIMPLEX->ID*))
(defun remove-tl1simplices (tl1sxs)
  (dolist (tl1sx tl1sxs)
    (remove-tl1simplex tl1sx)))
(defun remove-sl1simplices (sl1sxs)
  (dolist (sl1sx sl1sxs)
    (remove-sl1simplex sl1sx)))

(defun make-3simplex (sxtype tmlo tmhi p0 p1 p2 p3)
  (let ((sx3id (next-3simplex-id)))
    (setf (gethash sx3id *ID->3SIMPLEX*)
	  (list sxtype tmlo tmhi (list p0 p1 p2 p3) (list 0 0 0 0)))
    (make-2simplices sxtype tmlo tmhi p0 p1 p2 p3)
    (make-1simplices sxtype tmlo tmhi p0 p1 p2 p3)
    sx3id))

;; same as above except the points are "packed" into a list
(defun make-3simplex-v2 (sxtype tmlo tmhi pts)
    (let ((p0 (nth 0 pts))
	(p1 (nth 1 pts))
	(p2 (nth 2 pts))
	(p3 (nth 3 pts)))
    (make-3simplex sxtype tmlo tmhi p0 p1 p2 p3)))

;; this version is used only during initialization. If periodic b.c. are 
;; specified, it adjusts the points on the final time slice, since the t=T 
;; slice is identified with t=0 slice.
(defun make-3simplex-v3 (sxtype tmlo tmhitmp p0tmp p1tmp p2tmp p3tmp)
  (let ((p0 p0tmp) (p1 p1tmp) (p2 p2tmp) (p3 p3tmp) (tmhi tmhitmp))
    (when (and (string= BCTYPE "PERIODIC") (= NUM-T tmhitmp))
      (setf tmhi 0)
      (cond ((= 1 sxtype)
	     (decf p1 (* N0-PER-SLICE NUM-T)) 
	     (decf p2 (* N0-PER-SLICE NUM-T)) 
	     (decf p3 (* N0-PER-SLICE NUM-T)))
	    ((= 2 sxtype)
	     (decf p2 (* N0-PER-SLICE NUM-T)) 
	     (decf p3 (* N0-PER-SLICE NUM-T)))
	    ((= 3 sxtype)
	     (decf p3 (* N0-PER-SLICE NUM-T)))))
    (make-3simplex sxtype tmlo tmhi p0 p1 p2 p3)))

;; this version is used for loading the simplex data from file
(defun make-3simplex-v4 (sxtype tmlo tmhi p0 p1 p2 p3 n0 n1 n2 n3 sx3id)
  (setf (gethash sx3id *ID->3SIMPLEX*)
	(list sxtype tmlo tmhi (list p0 p1 p2 p3) (list n0 n1 n2 n3)))
  (make-2simplices sxtype tmlo tmhi p0 p1 p2 p3)
  (make-1simplices sxtype tmlo tmhi p0 p1 p2 p3))

;; all the simplex data, not just the points, is packed into a list
(defun make-3simplex-v5 (simplex-data)
  (let ((sxtype (first simplex-data))
	(tmlo (second simplex-data))
	(tmhi (third simplex-data))
	(pts (fourth simplex-data)))
    (make-3simplex-v2 sxtype tmlo tmhi pts)))

;; simplex-data-list = ((typ tmlo tmhi (p0 p1 p2 p3))...)
;; the ids of the simplices are returned
(defun make-3simplices-in-bulk (simplex-data-list)
  (let ((3sxids nil))
    (dolist (simplex-data simplex-data-list)
      (push (make-3simplex-v5 simplex-data) 3sxids))
    3sxids))

(defmacro 3sx-type (sx) `(first ,sx))
(defmacro 3sx-tmlo (sx) `(second ,sx))
(defmacro 3sx-tmhi (sx) `(third ,sx))
(defmacro 3sx-points (sx) `(fourth ,sx))
(defmacro 3sx-sx3ids (sx) `(fifth ,sx))
(defmacro 3sx-lopts (sx) `(subseq (3sx-points ,sx) 0 (3sx-type ,sx)))
(defmacro 3sx-hipts (sx) `(subseq (3sx-points ,sx) (3sx-type ,sx)))
(defmacro nth-point (sx n) `(nth ,n (3sx-points ,sx)))
(defmacro nth-neighbor (sx n) `(nth ,n (3sx-sx3ids ,sx)))

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
  (maphash #'(lambda (3sxid 3sx) 
	       (cond ((= 1 (3sx-type 3sx))
		      (format t "[~A] (~A ~A ~A (~A|~A ~A ~A) (~A ~A ~A ~A))~%"
			      3sxid (3sx-type 3sx) (3sx-tmlo 3sx) (3sx-tmhi 3sx)
			      (nth-point 3sx 0) (nth-point 3sx 1) 
			      (nth-point 3sx 2) (nth-point 3sx 3)
			      (nth-neighbor 3sx 0) (nth-neighbor 3sx 1) 
			      (nth-neighbor 3sx 2) (nth-neighbor 3sx 3)))
		     ((= 2 (3sx-type 3sx))
		      (format t "[~A] (~A ~A ~A (~A ~A|~A ~A) (~A ~A ~A ~A))~%"
			      3sxid (3sx-type 3sx) (3sx-tmlo 3sx) (3sx-tmhi 3sx)
			      (nth-point 3sx 0) (nth-point 3sx 1) 
			      (nth-point 3sx 2) (nth-point 3sx 3)
			      (nth-neighbor 3sx 0) (nth-neighbor 3sx 1) 
			      (nth-neighbor 3sx 2) (nth-neighbor 3sx 3)))
		     ((= 3 (3sx-type 3sx))
		      (format t "[~A] (~A ~A ~A (~A ~A ~A|~A) (~A ~A ~A ~A))~%"
			      3sxid (3sx-type 3sx) (3sx-tmlo 3sx) (3sx-tmhi 3sx)
			      (nth-point 3sx 0) (nth-point 3sx 1) 
			      (nth-point 3sx 2) (nth-point 3sx 3)
			      (nth-neighbor 3sx 0) (nth-neighbor 3sx 1) 
			      (nth-neighbor 3sx 2) (nth-neighbor 3sx 3)))))
	   *ID->3SIMPLEX*))

(defun show-tl2simplex-store ()
  (let ((count 1))
    (maphash #'(lambda (tl2sx id)
		 (format t "[~A] ~A ~A~%" count tl2sx id)
		 (incf count))
	     *TL2SIMPLEX->ID*)))
(defun show-sl2simplex-store ()
  (let ((count 1))
    (maphash #'(lambda (sl2sx id)
		 (format t "[~A] ~A ~A~%" count sl2sx id)
		 (incf count))
	     *SL2SIMPLEX->ID*)))

(defun show-tl1simplex-store ()
  (let ((count 1))
    (maphash #'(lambda (tl1sx id)
		 (format t "[~A] ~A ~A~%" count tl1sx id)
		 (incf count))
	     *TL1SIMPLEX->ID*)))
(defun show-sl1simplex-store ()
  (let ((count 1))
    (maphash #'(lambda (sl1sx id)
		 (format t "[~A] ~A ~A~%" count sl1sx id)
		 (incf count))
	     *SL1SIMPLEX->ID*)))

(defun connect-3simplices (sx1id sx2id)
  (let ((sx1 nil) (sx2 nil))
    (when (and (setf sx1 (get-3simplex sx1id))
	       (setf sx2 (get-3simplex sx2id)))
      (let ((tri (intersection (3sx-points sx1) (3sx-points sx2))))
	(when (= 3 (length tri))
	  (let* ((pts1 (3sx-points sx1))
		 (pts2 (3sx-points sx2))
		 (pos1 (position (first (set-difference pts1 tri)) pts1))
		 (pos2 (position (first (set-difference pts2 tri)) pts2)))
	    (setf (nth pos1 (3sx-sx3ids sx1)) sx2id 
		  (nth pos2 (3sx-sx3ids sx2)) sx1id)))))))

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
     (and (setf sx1 (get-3simplex ,sxid1)) 
	  (setf sx2 (get-3simplex ,sxid2))
	  (find ,sxid1 (3sx-sx3ids sx2)) 
	  (find ,sxid2 (3sx-sx3ids sx1)))))

(defun get-simplices-in-sandwich (tlo thi)
  (let ((sxids '()))
    (maphash #'(lambda (id sx)
		 (when (and (= (3sx-tmlo sx) (mod tlo NUM-T)) 
			    (= (3sx-tmhi sx) (mod thi NUM-T)))
		   (push id sxids)))
	     *ID->3SIMPLEX*)
    sxids))

(defun get-simplices-in-sandwich-of-type (tlo thi typ)
  (let ((sxids '()))
    (maphash #'(lambda (id sx)
		 (when (and (= (3sx-tmlo sx) (mod tlo NUM-T)) 
			    (= (3sx-tmhi sx) (mod thi NUM-T))
			    (= (3sx-type sx) typ))
		   (push id sxids)))
	     *ID->3SIMPLEX*)
    sxids))

;; returns (1ids 2ids 3ids) where 1ids is a list of (1,3) ids in the sandwich
;;  etc.
(defun get-simplices-in-sandwich-ordered-by-type (tlo thi)
  (let ((1ids nil) (2ids nil) (3ids nil))
    (maphash #'(lambda (id sx)
		 (when (and (= (3sx-tmlo sx) (mod tlo NUM-T)) 
			    (= (3sx-tmhi sx) (mod thi NUM-T)))
		   (ecase (3sx-type sx)
		     (1 (push id 1ids))
		     (2 (push id 2ids))
		     (3 (push id 3ids)))))
	     *ID->3SIMPLEX*)
    (values 1ids 2ids 3ids)))

(defun get-simplices-of-type (typ)
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
		 (if (= (3sx-type sx) typ)
		     (incf count)))
	     *ID->3SIMPLEX*)
    count))

(defun count-simplices-in-sandwich (tlo thi)
  (let ((count 0))
    (maphash #'(lambda (id sx)
		 (declare (ignore id))
		 (when (and (= (3sx-tmlo sx) (mod tlo NUM-T)) 
			    (= (3sx-tmhi sx) (mod thi NUM-T)))
		     (incf count)))
	     *ID->3SIMPLEX*)
    count))

(defun count-simplices-of-all-types ()
  (let ((1count 0) (2count 0) (3count 0))
    (maphash #'(lambda (id sx)
		 (declare (ignore id))
		 (ecase (3sx-type sx)
		   (1 (incf 1count))
		   (2 (incf 2count))
		   (3 (incf 3count))))
	     *ID->3SIMPLEX*)
    (list 1count 2count 3count (+ 1count 2count 3count))))

;; in a given sandwich
;; (1,3) can be connected to a (2,2) and cannot be connected to a (3,1)
;; a (2,2) can be connected to a (1,3) and a (3,1)
;; a (3,2) can be connected to a (2,2) and cannot be connected to a (1,3)
(defun connect-simplices-in-sandwich (tlo thi)
  (multiple-value-bind (1ids 2ids 3ids) 
      (get-simplices-in-sandwich-ordered-by-type tlo thi)
    (connect-3simplices-across-lists 1ids 2ids)
    (connect-3simplices-across-lists 2ids 3ids)
    (connect-3simplices-within-list 1ids)
    (connect-3simplices-within-list 2ids)
    (connect-3simplices-within-list 3ids)))

(defun connect-simplices-in-adjacent-sandwiches (tl tm th)
  (connect-3simplices-across-lists 
   (get-simplices-in-sandwich-of-type tl tm 1)
   (get-simplices-in-sandwich-of-type tm th 3)))

(defun neighbors-of-type (sx sxtype)
  (let ((nbors nil)
	(nsx nil)
	(nids (3sx-sx3ids sx)))
    (for (n 0 3)
      (when (and (setf nsx (get-3simplex (nth n nids))) 
		 (= sxtype (3sx-type nsx)))
	(pushnew (nth n nids) nbors)))
    nbors))

(defun save-spacetime-to-file (outfile)
  (format outfile "~A ~A ~A ~A ~A ~A ~A ~A ~A ~A ~A ~A ~A ~A ~A ~A ~A~%" 
	  BCTYPE STOPOLOGY NUM-T N-INIT *LAST-USED-POINT* *LAST-USED-3SXID* 
	  N0 N1-SL N1-TL N2-SL N2-TL N3-TL-31 N3-TL-22 *eps* *k0* *k3* *alpha*)
  (maphash #'(lambda (k v)
	       (format outfile "~A ~A ~A ~A ~A ~A ~A ~A ~A ~A ~A ~A~%" 
		       (3sx-type v) (3sx-tmlo v) (3sx-tmhi v)
		       (nth-point v 0) (nth-point v 1) 
		       (nth-point v 2) (nth-point v 3)
		       (nth-neighbor v 0) (nth-neighbor v 1)
		       (nth-neighbor v 2) (nth-neighbor v 3) k))
	   *ID->3SIMPLEX*))

(defun parse-parameters-line (line)
  (with-input-from-string (s line)
    (let ((data (loop
		   :for num := (read s nil nil)
		   :while num
		   :collect num)))
      (setf BCTYPE (nth 0 data) STOPOLOGY (nth 1 data) NUM-T (nth 2 data) 
	    N-INIT (nth 3 data) *LAST-USED-POINT* (nth 4 data) 
	    *LAST-USED-3SXID* (nth 5 data) N0 (nth 6 data) N1-SL (nth 7 data) 
	    N1-TL (nth 8 data) N2-SL (nth 9 data) N2-TL (nth 10 data)
	    N3-TL-31 (nth 11 data) N3-TL-22 (nth 12 data) *eps* (nth 13 data))
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

;; spatial 2-simplex is a spatial triangle; this information is needed for 
;; computing the spectral and hausdorff dimensions of the spatial slice, among 
;; other things.
(defun make-s2simplex (triangle-time triangle-points)
  (let ((stid (next-s2simplex-id)))
    (setf (gethash stid *ID->SPATIAL-2SIMPLEX*) 
	  (list triangle-time (copy-list triangle-points) (list 0 0 0)))
    stid))

;; this version is used for loading the simplex data from file
(defun make-s2simplex-v2 (s2sxid s2sxtm p0 p1 p2 n0 n1 n2)
  (setf (gethash s2sxid *ID->SPATIAL-2SIMPLEX*)
	(list s2sxtm (list p0 p1 p2) (list n0 n1 n2))))

(defmacro get-s2simplex (stid) `(gethash ,stid *ID->SPATIAL-2SIMPLEX*))
(defmacro s2sx-time (s2sx) `(first ,s2sx))
(defmacro s2sx-points (s2sx) `(second ,s2sx))
(defmacro s2sx-sx2ids (s2sx) `(third ,s2sx))


(defun connect-s2simplices (st1id st2id)
  (let ((st1 nil) (st2 nil))
    (when (and (setf st1 (get-s2simplex st1id)) 
	       (setf st2 (get-s2simplex st2id)))
      (let* ((points1 (s2sx-points st1))
	     (points2 (s2sx-points st2))
	     (line (intersection points1 points2)))
	(when (= 2 (length line))
	  (let ((pos1 (position (first (set-difference points1 line)) points1))
		(pos2 (position (first (set-difference points2 line)) points2)))
	    (setf (nth pos1 (s2sx-sx2ids st1)) st2id 
		  (nth pos2 (s2sx-sx2ids st2)) st1id)))))))

(defun connect-s2simplices-within-list (sx1ids)
  (for (n 0 (1- (length sx1ids)))
       (for (m (1+ n) (1- (length sx1ids)))
	    (connect-s2simplices (nth n sx1ids) (nth m sx1ids)))))

(defun count-s2simplices-in-slice (ts)
  "counts the number of spaeclike 2-simplices (i.e. spatial triangles) in spatial slice ts"
  (let ((count 0))
    (maphash #'(lambda (id sx)
		 (declare (ignore id))
		 (when (= (s2sx-time sx) (mod ts NUM-T)) 
		   (incf count)))
	     *ID->SPATIAL-2SIMPLEX*)
    count))

(defun parse-s2simplex-parameters-line (line)
  "parses line for parameters but does nothing with them. The parameters are
stored in the 2-simplex data file for purposes of identification only"
  (with-input-from-string (s line)))

(defun parse-s2simplex-data-line (line)
  "parses line for spatial 2-simplex data"
  (with-input-from-string (s line)
    (let ((data (loop
		   :for num := (read s nil nil)
		   :while num
		   :collect num)))
      (make-s2simplex-v2 (nth 0 data) (nth 1 data) (nth 2 data) (nth 3 data) 
			 (nth 4 data) (nth 5 data) (nth 6 data) (nth 7 data)))))

(defun load-s2simplex-data-from-file (infile)
  "loads spatial 2-simplex data from infile"
  (parse-s2simplex-parameters-line (read-line infile nil))
  (clrhash *ID->SPATIAL-2SIMPLEX*)
  (loop for line = (read-line infile nil)
     while line do (parse-s2simplex-data-line line)))

(defun save-s2simplex-data-to-file (outfile)
  "saves the spatial 2-simplex data to outfile"
  (format outfile "~A ~A ~A ~A ~A ~A ~A ~A~%" BCTYPE STOPOLOGY NUM-T N-INIT 
	  N2-SL (hash-table-count *ID->SPATIAL-2SIMPLEX*) *k0* *k3*)
  (maphash #'(lambda (k v)
	       (let ((pts (s2sx-points v))
		     (nbors (s2sx-sx2ids v)))
		 (format outfile "~A ~A ~A ~A ~A ~A ~A ~A~%"
			 k (s2sx-time v) (nth 0 pts) (nth 1 pts) (nth 2 pts) 
			 (nth 0 nbors) (nth 1 nbors) (nth 2 nbors))))
	   *ID->SPATIAL-2SIMPLEX*))

(defun 3sx2p1->s2sx2p1 ()
  "3sx2p1->s2sx2p1 generates the spatial 2-simplex information for each spatial slice 
from the 3-simplex data for the entire spacetime."
  (clrhash *ID->SPATIAL-2SIMPLEX*)
  (setf *LAST-USED-S2SXID* 0)
  (for (ts 0 (1- NUM-T))
       (let ((31simplices (get-simplices-in-sandwich-of-type ts (1+ ts) 3)) ;; list of ids
	     (spatial-triangles '()))
	 (dolist (31simplex 31simplices)
	   (push (make-s2simplex ts (3sx-lopts (get-3simplex 31simplex))) 
		 spatial-triangles))
	 (connect-s2simplices-within-list spatial-triangles))))

(defun generate-s2sx2p1-files (3sx2p1files)
  "3sx2p1files is a list of .3sx2p1 files. For each file in this list, this
function generates a .s2sx2p1 file. The prefix for the .3sx2p1 and the .s2sx2p1
file are identical"
  (loop for line = (read-line 3sx2p1files nil)
       while line do
       (let ((outfilename (change-file-suffix line "s2sx2p1")))
	 (with-open-file (indatafile line :direction :input)
	   (load-spacetime-from-file indatafile))
	 (3sx2p1->s2sx2p1)
	 (clrhash *ID->3SIMPLEX*)
	 (sb-ext:gc :gen 1000)
	 (with-open-file (outdatafile outfilename :direction :output)
	   (save-s2simplex-data-to-file outdatafile))
	 (format t "finished 3sx2p1->s2sx2p1 for ~A at ~A~%"
		 outfilename (cdt-now-str)))))

	     
