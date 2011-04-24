;; 3simplex
;; (type tmlo tmhi (p0 p1 p2 p3))
;; type = 0,1,2,3,4 tm[lo|hi] = [lo|hi] spatial time slice pj = points

;; 4simplex
;; (type tmlo tmhi (p0 p1 p2 p3 p4) (n0 n1 n2 n3 n4) (t0 t1 t2 t3 t4))
;; type = 1,2,3,4 tm[lo|hi] - [lo|hi] spatial time slice
;; tj = id of the 2sx that does not have pj
;; nj = id of the 3sx that does not have pj

(defun make-3simplex (type tmlo tmhi p0 p1 p2 p3)
  "makes and returns the id of the 3-simplex with the specified data. If a 3-simplex with points 
p0 p1 p2 p3 already exists, the id of that simplex is returned"
  (let* ((3sx (list type tmlo tmhi (list p0 p1 p2 p3)))
	 (3sxid (gethash 3sx *3SIMPLEX->ID*)))
    (unless 3sxid
      (setf 3sxid (next-3simplex-id))
      (setf (gethash 3sx *3SIMPLEX->ID*) 3sxid)
      (setf (gethash 3sxid *ID->3SIMPLEX*) 3sx))
    3sxid))

(defmacro get-3simplex (sxid) `(gethash ,sxid *ID->3SIMPLEX*))

(defmacro 3sx-type (sx) `(first ,sx))
(defmacro 3sx-tmlo (sx) `(second ,sx))
(defmacro 3sx-tmhi (sx) `(third ,sx))
(defmacro 3sx-points (sx) `(fourth ,sx))
(defmacro 3sx-lopts (sx) `(subseq (3sx-points ,sx) 0 (3sx-type ,sx)))
(defmacro 3sx-hipts (sx) `(subseq (3sx-points ,sx) (3sx-type ,sx)))

(defmacro remove-3simplex (3sxid)
  `(let ((3sx (gethash ,3sxid *ID->3SIMPLEX*)))
     (remhash ,3sxid *ID->3SIMPLEX*)
     (remhash 3sx *3SIMPLEX->ID*)))

(defun remove-3simplices (3sxids)
  (dolist (3sxid 3sxids)
    (let ((3sx (gethash 3sxid *ID->3SIMPLEX*)))
      (remhash 3sxid *ID->3SIMPLEX*)
      (remhash 3sx *3SIMPLEX->ID*))))

(defun show-3simplex->id-store ()
  (maphash #'(lambda (3sx 3sxid) (format t "~A [~A]~%" 3sx 3sxid)) *3SIMPLEX->ID*))

(defun show-id->3simplex-store ()
  (maphash #'(lambda (3sxid 3sx) (format t "[~A] ~A~%" 3sxid 3sx)) *ID->3SIMPLEX*))

(defun make-4simplex (type tmlo tmhi p0 p1 p2 p3 p4)
  (let ((t0type nil) (t1type nil) (t2type nil) (t3type nil) (t4type nil) (sx4id (next-4simplex-id)))
    (ecase type
      (1 (setf t0type 0 t1type 1 t2type 1 t3type 1 t4type 1))
      (2 (setf t0type 1 t1type 1 t2type 2 t3type 2 t4type 2))
      (3 (setf t0type 2 t1type 2 t2type 2 t3type 3 t4type 3))
      (4 (setf t0type 3 t1type 3 t2type 3 t3type 3 t4type 4)))
    (setf (gethash sx4id *ID->4SIMPLEX*)
	  (list type tmlo tmhi
		(list p0 p1 p2 p3 p4)
		(list 0 0 0 0 0)
		(list (make-3simplex t0type tmlo tmhi p1 p2 p3 p4)
		      (make-3simplex t1type tmlo tmhi p0 p2 p3 p4)
		      (make-3simplex t2type tmlo tmhi p0 p1 p3 p4)
		      (make-3simplex t3type tmlo tmhi p0 p1 p2 p4)
		      (make-3simplex t4type tmlo tmhi p0 p1 p2 p3))))
    sx4id))

;; same as make-4simplex except the points are packed in a list called pts
(defun make-4simplex-v2 (type tmlo tmhi pts)
  (let ((t0type nil) (t1type nil) (t2type nil) (t3type nil) (t4type nil) 
	(sx4id (next-4simplex-id)))
    (ecase type
      (1 (setf t0type 0 t1type 1 t2type 1 t3type 1 t4type 1))
      (2 (setf t0type 1 t1type 1 t2type 2 t3type 2 t4type 2))
      (3 (setf t0type 2 t1type 2 t2type 2 t3type 3 t4type 3))
      (4 (setf t0type 3 t1type 3 t2type 3 t3type 3 t4type 4)))
    (setf (gethash sx4id *ID->4SIMPLEX*)
	  (list type tmlo tmhi
		(copy-list pts)
		(list 0 0 0 0 0)
		(list (make-3simplex 
		       t0type tmlo tmhi 
		       (nth 1 pts) (nth 2 pts) (nth 3 pts) (nth 4 pts))
		      (make-3simplex 
		       t1type tmlo tmhi 
		       (nth 0 pts) (nth 2 pts) (nth 3 pts) (nth 4 pts))
		      (make-3simplex 
		       t2type tmlo tmhi 
		       (nth 0 pts) (nth 1 pts) (nth 3 pts) (nth 4 pts))
		      (make-3simplex 
		       t3type tmlo tmhi 
		       (nth 0 pts) (nth 1 pts) (nth 2 pts) (nth 4 pts))
		      (make-3simplex 
		       t4type tmlo tmhi 
		       (nth 0 pts) (nth 1 pts) (nth 2 pts) (nth 3 pts)))))
    sx4id))

;; this version is used only during initialization. If periodic b.c. are 
;; specified, it adjusts the points on the final time slice, since the t=T 
;; slice is identified with t=0 slice.
(defun make-4simplex-v3 (type tmlo tmhitmp p0tmp p1tmp p2tmp p3tmp p4tmp)
  (let ((t0type nil) (t1type nil) (t2type nil) (t3type nil) (t4type nil) 
	(sx4id (next-4simplex-id))
	(p0 p0tmp) (p1 p1tmp) (p2 p2tmp) (p3 p3tmp) (p4 p4tmp) (tmhi tmhitmp))
    (when (and (string= BCTYPE "PERIODIC") (= NUM-T tmhitmp))
      (setf tmhi 0)
      (cond ((= 1 type)
	     (decf p1 (* 5 NUM-T)) (decf p2 (* 5 NUM-T)) (decf p3 (* 5 NUM-T)) 
	     (decf p4 (* 5 NUM-T)))
	    ((= 2 type)
	     (decf p2 (* 5 NUM-T)) (decf p3 (* 5 NUM-T)) (decf p4 (* 5 NUM-T)))
	    ((= 3 type)
	     (decf p3 (* 5 NUM-T)) (decf p4 (* 5 NUM-T)))
	    ((= 4 type)
	     (decf p4 (* 5 NUM-T)))))
    (ecase type
      (1 (setf t0type 0 t1type 1 t2type 1 t3type 1 t4type 1))
      (2 (setf t0type 1 t1type 1 t2type 2 t3type 2 t4type 2))
      (3 (setf t0type 2 t1type 2 t2type 2 t3type 3 t4type 3))
      (4 (setf t0type 3 t1type 3 t2type 3 t3type 3 t4type 4)))
    (setf (gethash sx4id *ID->4SIMPLEX*)
	  (list type tmlo tmhi
		(list p0 p1 p2 p3 p4)
		(list 0 0 0 0 0)
		(list (make-3simplex t0type tmlo tmhi p1 p2 p3 p4)
		      (make-3simplex t1type tmlo tmhi p0 p2 p3 p4)
		      (make-3simplex t2type tmlo tmhi p0 p1 p3 p4)
		      (make-3simplex t3type tmlo tmhi p0 p1 p2 p4)
		      (make-3simplex t4type tmlo tmhi p0 p1 p2 p3))))
    sx4id))

;; this version is used for loading the simplex data from file
(defun make-4simplex-v4 (type tmlo tmhi p0 p1 p2 p3 p4 n0 n1 n2 n3 n4 sx4id)
  (let ((t0type nil) (t1type nil) (t2type nil) (t3type nil) (t4type nil))
    (ecase type
      (1 (setf t0type 0 t1type 1 t2type 1 t3type 1 t4type 1))
      (2 (setf t0type 1 t1type 1 t2type 2 t3type 2 t4type 2))
      (3 (setf t0type 2 t1type 2 t2type 2 t3type 3 t4type 3))
      (4 (setf t0type 3 t1type 3 t2type 3 t3type 3 t4type 4)))
    (setf (gethash sx4id *ID->4SIMPLEX*)
	  (list type tmlo tmhi
		(list p0 p1 p2 p3 p4)
		(list n0 n1 n2 n3 n4)
		(list (make-3simplex t0type tmlo tmhi p1 p2 p3 p4)
		      (make-3simplex t1type tmlo tmhi p0 p2 p3 p4)
		      (make-3simplex t2type tmlo tmhi p0 p1 p3 p4)
		      (make-3simplex t3type tmlo tmhi p0 p1 p2 p4)
		      (make-3simplex t4type tmlo tmhi p0 p1 p2 p3))))))

;; a replacement for v2
(defun make-4simplex-v5 (simplex-data)
  (let ((type (first simplex-data))
	(tmlo (second simplex-data))
	(tmhi (third simplex-data))
	(pts (fourth simplex-data))
	(t0type nil) (t1type nil) (t2type nil) (t3type nil) (t4type nil) 
	(sx4id (next-4simplex-id)))
    (ecase type
      (1 (setf t0type 0 t1type 1 t2type 1 t3type 1 t4type 1))
      (2 (setf t0type 1 t1type 1 t2type 2 t3type 2 t4type 2))
      (3 (setf t0type 2 t1type 2 t2type 2 t3type 3 t4type 3))
      (4 (setf t0type 3 t1type 3 t2type 3 t3type 3 t4type 4)))
    (setf (gethash sx4id *ID->4SIMPLEX*)
	  (list type tmlo tmhi
		(copy-list pts)
		(list 0 0 0 0 0)
		(list (make-3simplex 
		       t0type tmlo tmhi 
		       (nth 1 pts) (nth 2 pts) (nth 3 pts) (nth 4 pts))
		      (make-3simplex 
		       t1type tmlo tmhi 
		       (nth 0 pts) (nth 2 pts) (nth 3 pts) (nth 4 pts))
		      (make-3simplex 
		       t2type tmlo tmhi 
		       (nth 0 pts) (nth 1 pts) (nth 3 pts) (nth 4 pts))
		      (make-3simplex 
		       t3type tmlo tmhi 
		       (nth 0 pts) (nth 1 pts) (nth 2 pts) (nth 4 pts))
		      (make-3simplex 
		       t4type tmlo tmhi 
		       (nth 0 pts) (nth 1 pts) (nth 2 pts) (nth 3 pts)))))
    sx4id))


;; simplex-data-list = ((typ tmlo tmhi (p0 p1 p2 p3 p4))...)
;; the ids of the simplices are returned
(defun make-4simplices-in-bulk (simplex-data-list)
  (let ((4sxids nil))
    (dolist (simplex-data simplex-data-list)
      (push (make-4simplex-v5 simplex-data) 4sxids))
    4sxids))

(defmacro 4sx-type (sx) `(first ,sx))
(defmacro 4sx-tmlo (sx) `(second ,sx))
(defmacro 4sx-tmhi (sx) `(third ,sx))
(defmacro 4sx-points (sx) `(fourth ,sx))
(defmacro 4sx-sx4ids (sx) `(fifth ,sx))
(defmacro 4sx-sx3ids (sx) `(sixth ,sx))
(defmacro 4sx-lopts (sx) `(subseq (4sx-points ,sx) 0 (4sx-type ,sx)))
(defmacro 4sx-hipts (sx) `(subseq (4sx-points ,sx) (4sx-type ,sx)))

(defmacro get-4simplex (sxid) `(gethash ,sxid *ID->4SIMPLEX*))

(defmacro remove-4simplex (4sxid)
  `(progn
     (remhash ,4sxid *ID->4SIMPLEX*)
     (recycle-4simplex-id ,4sxid)))
(defmacro remove-4simplices (4sxids)
  `(dolist (4sxid ,4sxids)
     (remhash 4sxid *ID->4SIMPLEX*)
     (recycle-4simplex-id 4sxid)))

(defun show-id->4simplex-store ()
  (maphash #'(lambda (4sxid 4sx) (format t "[~A] ~A~%" 4sxid 4sx)) 
	   *ID->4SIMPLEX*))

(defmacro nth-point (sx n)
  `(nth ,n (4sx-points ,sx)))

(defun connect-4simplices (sx1id sx2id)
  (let ((sx1 nil) (sx2 nil))
    (when (and (setf sx1 (get-4simplex sx1id)) (setf sx2 (get-4simplex sx2id)))
      (let ((3sxlinkid (intersection (4sx-sx3ids sx1) (4sx-sx3ids sx2))))
	(when (= 1 (length 3sxlinkid))
	  (let ((pos1 (position (first 3sxlinkid) (4sx-sx3ids sx1)))
		(pos2 (position (first 3sxlinkid) (4sx-sx3ids sx2))))
	    (setf (nth pos1 (4sx-sx4ids sx1)) sx2id 
		  (nth pos2 (4sx-sx4ids sx2)) sx1id)))))))

(defun connect-4simplices-within-list (sx1ids)
  (for (n 0 (1- (length sx1ids)))
       (for (m (1+ n) (1- (length sx1ids)))
	    (connect-4simplices (nth n sx1ids) (nth m sx1ids)))))

(defun connect-4simplices-across-lists (sx1ids sx2ids)
  (dolist (sx1id sx1ids)
    (dolist (sx2id sx2ids)
      (connect-4simplices sx1id sx2id))))

(defmacro 4simplices-connected? (sxid1 sxid2)
  `(let ((foosx1 nil) (foosx2 nil))
     (and (setf foosx1 (get-4simplex ,sxid1)) 
	  (setf foosx2 (get-4simplex ,sxid2))
	  (find ,sxid1 (4sx-sx4ids foosx2)) 
	  (find ,sxid2 (4sx-sx4ids foosx1)))))

(defun get-simplices-in-sandwich (tlo thi)
  (let ((sxids '()))
    (maphash #'(lambda (id sx)
		 (when (and (= (4sx-tmlo sx) (mod tlo NUM-T)) 
			    (= (4sx-tmhi sx) (mod thi NUM-T)))
		   (push id sxids)))
	     *ID->4SIMPLEX*)
    sxids))

(defun get-simplices-in-sandwich-of-type (tlo thi typ)
  (let ((sxids '()))
    (maphash #'(lambda (id sx)
		 (when (and (= (4sx-tmlo sx) (mod tlo NUM-T)) 
			    (= (4sx-tmhi sx) (mod thi NUM-T))
			    (= (4sx-type sx) typ))
		   (push id sxids)))
	     *ID->4SIMPLEX*)
    sxids))

(defun get-3simplices-in-sandwich-of-type (tlo thi typ)
  (let ((sxids '()))
    (maphash #'(lambda (id sx)
		 (when (and (= (3sx-tmlo sx) (mod tlo NUM-T)) 
			    (= (3sx-tmhi sx) (mod thi NUM-T))
			    (= (3sx-type sx) typ))
		   (push id sxids)))
	     *ID->3SIMPLEX*)
    sxids))

;; returns (1ids 2ids 3ids 4ids) where 1ids is a list of (1,4) ids in the 
;; sandwich etc.
(defun get-simplices-in-sandwich-ordered-by-type (tlo thi)
  (let ((1ids nil) (2ids nil) (3ids nil) (4ids nil))
    (maphash #'(lambda (id sx)
		 (when (and (= (4sx-tmlo sx) (mod tlo NUM-T)) 
			    (= (4sx-tmhi sx) (mod thi NUM-T)))
		   (ecase (4sx-type sx)
		     (1 (push id 1ids))
		     (2 (push id 2ids))
		     (3 (push id 3ids))
		     (4 (push id 4ids)))))
	     *ID->4SIMPLEX*)
    (values 1ids 2ids 3ids 4ids)))


(defun get-simplices-of-type (typ)
  (let ((sxids '()))
    (maphash #'(lambda (id sx)
		 (if (= (4sx-type sx) typ)
		     (push id sxids)))
	     *ID->4SIMPLEX*)
    sxids))

(defun count-simplices-of-type (typ)
  (let ((count 0))
    (maphash #'(lambda (id sx)
		 (declare (ignore id))
		 (if (= (4sx-type sx) typ)
		     (incf count)))
	     *ID->4SIMPLEX*)
    count))

(defun count-simplices-in-sandwich (tlo thi)
  (let ((count 0))
    (maphash #'(lambda (id sx)
		 (declare (ignore id))
		 (when (and (= (4sx-tmlo sx) (mod tlo NUM-T)) 
			    (= (4sx-tmhi sx) (mod thi NUM-T)))
		     (incf count)))
	     *ID->4SIMPLEX*)
    count))

(defun count-simplices-of-all-types ()
  (let ((1count 0) (2count 0) (3count 0) (4count 0))
    (maphash #'(lambda (id sx)
		 (declare (ignore id))
		 (ecase (4sx-type sx)
		   (1 (incf 1count))
		   (2 (incf 2count))
		   (3 (incf 3count))
		   (4 (incf 4count))))
	     *ID->4SIMPLEX*)
    (list 1count 2count 3count 4count (+ 1count 2count 3count 4count))))

;;(defun connect-simplices-in-sandwich (tlo thi)
;;  (connect-4simplices-within-list (get-simplices-in-sandwich tlo thi)))

;; in a given sandwich,
;; (1,4) can be connected to (2,3) and cannot be connected to (3,2) or (4,1)
;; (2,3) can be connected to (3,2) and (1,4) and cannot be connected to (4,1)
;; (3,2) can be connected to (2,3) and (4,1) and cannot be connected to (1,4)
;; (4,1) can be connected to (3,2) and cannot be connected to (2,3) or (1,4)
(defun connect-simplices-in-sandwich (tlo thi)
  (multiple-value-bind (1ids 2ids 3ids 4ids) 
      (get-simplices-in-sandwich-ordered-by-type tlo thi)
    (connect-4simplices-across-lists 1ids 2ids)
    (connect-4simplices-across-lists 2ids 3ids)
    (connect-4simplices-across-lists 3ids 4ids)
    (connect-4simplices-within-list 1ids)
    (connect-4simplices-within-list 2ids)
    (connect-4simplices-within-list 3ids)
    (connect-4simplices-within-list 4ids)))

(defun connect-simplices-in-adjacent-sandwiches (tl tm th)
  (connect-4simplices-across-lists 
   (get-simplices-in-sandwich-of-type tl tm 1)
   (get-simplices-in-sandwich-of-type tm th 4)))

(defun check-14-and-41 (tlo thi)
  (let ((14ids (get-simplices-in-sandwich-of-type tlo thi 1))
	(41ids (get-simplices-in-sandwich-of-type tlo thi 4))
	(problem-ids '()))
    (dolist (s14 14ids)
      (dolist (d14 14ids)
	(when (and (/= s14 d14) (set-equal? 
				 (subseq (4sx-points (get-4simplex s14)) 1)
				 (subseq (4sx-points (get-4simplex d14)) 1)))
	  (push (list s14 d14) problem-ids))))
    (dolist (s41 41ids)
      (dolist (d41 41ids)
	(when (and (/= s41 d41) (set-equal? 
				 (subseq (4sx-points (get-4simplex s41)) 0 4)
				 (subseq (4sx-points (get-4simplex d41)) 0 4)))
	  (push (list s41 d41) problem-ids))))
    problem-ids))

(defun check-all-slices-for-problem-simplices ()
  (for (ts 0 (- NUM-T 1))
       (format t "slice ~A has ~A problem simplices~%" ts 
	       (check-14-and-41 ts (1+ ts)))
       (finish-output)))

(defun check-all-slices-for-simplices-with-missing-neighbors ()
  (let ((problem-ids '()))
    (maphash #'(lambda (id sx)
		 (for (n 0 4)
		   (if (= 0 (nth n (4sx-sx4ids sx)))
		       (push id problem-ids))))
	     *ID->4SIMPLEX*)
    problem-ids))

(defmacro link-id (sxid1 sxid2)
  `(let ((sx1 nil) (sx2 nil) (link nil))
     (if (and (setf sx1 (get-4simplex ,sxid1)) 
	      (setf sx2 (get-4simplex ,sxid2))
	      (setf link (intersection (4sx-sx3ids sx2) (4sx-sx3ids sx1))))
	 (first link)
	 0)))

(defun neighbors-of-type (sx type)
  (let ((nbors nil)
	(nsx nil)
	(nids (4sx-sx4ids sx)))
    (for (n 0 4)
      (when (and (setf nsx (get-4simplex (nth n nids))) (= type (4sx-type nsx)))
	(pushnew (nth n nids) nbors)))
    nbors))

(defun save-spacetime-to-file (outfile)
  ;; the first line with simulation parameters
  (format outfile 
	  "~A ~A ~A ~A ~A ~A ~A ~A ~A ~A ~A ~A ~A ~A ~A ~A ~A ~A ~A ~A~%" 
	  BCTYPE STOPOLOGY NUM-T N-INIT *LAST-USED-POINT* *LAST-USED-4SXID* 
	  N0 N1-SL N1-TL N2-SL N2-TL N3-TL-31 N3-TL-22 N3-SL N4-TL-41 N4-TL-32 
	  EPS KAPPA-0 DELTA KAPPA-4)
  ;; each subsequent line has the information for a single 4-simplex in the 
  ;; following format
  ;; type tmlo tmhi p0 p1 p2 p3 p4 n0 n1 n2 n3 n4 id
  (maphash #'(lambda (k v)
	       (format outfile "~A ~A ~A ~A ~A ~A ~A ~A ~A ~A ~A ~A ~A ~A~%" 
		       (3sx-type v) (3sx-tmlo v) (3sx-tmhi v)
		       (nth-point v 0) (nth-point v 1) (nth-point v 2) 
		       (nth-point v 3) (nth-point v 4)
		       (nth 0 (4sx-sx4ids v)) (nth 1 (4sx-sx4ids v)) 
		       (nth 2 (4sx-sx4ids v)) 
		       (nth 3 (4sx-sx4ids v)) (nth 4 (4sx-sx4ids v)) k))
	   *ID->4SIMPLEX*))

(defun parse-parameters-line (line)
  (with-input-from-string (s line)
    (let ((data (loop
		   :for num := (read s nil nil)
		   :while num
		   :collect num)))
      (setf BCTYPE (nth 0 data))(setf STOPOLOGY (nth 1 data))
      (setf NUM-T (nth 2 data)) (setf N-INIT (nth 3 data))
      (setf *LAST-USED-POINT* (nth 4 data)) (setf *LAST-USED-4SXID* 
						  (nth 5 data))
      (setf N0 (nth 6 data)) (setf N1-SL (nth 7 data))
      (setf N1-TL (nth 8 data)) (setf N2-SL (nth 9 data))
      (setf N2-TL (nth 10 data)) (setf N3-TL-31 (nth 11 data))
      (setf N3-TL-22 (nth 12 data)) (setf N3-SL (nth 13 data))
      (setf N4-TL-41 (nth 14 data)) (setf N4-TL-32 (nth 15 data))
      (setf EPS (nth 16 data)) 
      (set-kappa0-delta-kappa4 (nth 17 data) (nth 18 data) (nth 19 data)))))

(defun parse-simplex-data-line (line)
  (with-input-from-string (s line)
    (let ((data (loop
		   :for num := (read s nil nil)
		   :while num
		   :collect num)))
      (make-4simplex-v4 (nth 0 data) (nth 1 data) (nth 2 data)
			(nth 3 data) (nth 4 data) (nth 5 data) (nth 6 data) 
			(nth 7 data)
			(nth 8 data) (nth 9 data) (nth 10 data) (nth 11 data) (nth 12 data)
			(nth 13 data)))))

(defun load-spacetime-from-file (infile)
  (parse-parameters-line (read-line infile nil))
  (loop for line = (read-line infile nil)
     while line do (parse-simplex-data-line line)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; spatial 3-simplex is a spatial tetrahedron; this information is needed for 
;; computing the spectral and hausdorff dimensions of the spatial slices, among 
;; other things.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
(defun make-s3simplex (tetrahedron-time tetrahedron-points)
  (let ((stid (next-s3simplex-id)))
    (setf (gethash stid *ID->SPATIAL-3SIMPLEX*) 
	  (list tetrahedron-time (copy-list tetrahedron-points) (list 0 0 0 0)))
    stid))

;; this version is used for loading the simplex data from file
(defun make-s3simplex-v2 (s3sxid s3sxtm p0 p1 p2 p3 n0 n1 n2 n3)
  (setf (gethash s3sxid *ID->SPATIAL-3SIMPLEX*)
	(list s3sxtm (list p0 p1 p2 p3) (list n0 n1 n2 n3))))

(defmacro get-s3simplex (stid) `(gethash ,stid *ID->SPATIAL-3SIMPLEX*))
(defmacro s3sx-time (s3sx) `(first ,s3sx))
(defmacro s3sx-points (s3sx) `(second ,s3sx))
(defmacro s3sx-sx3ids (s3sx) `(third ,s3sx))

(defun connect-s3simplices (st1id st2id)
  (let ((st1 nil) (st2 nil))
    (when (and (setf st1 (get-s3simplex st1id)) (setf st2 (get-s3simplex st2id)))
      (let* ((points1 (s3sx-points st1))
	     (points2 (s3sx-points st2))
	     (triangle (intersection points1 points2)))
	(when (= 3 (length triangle))
	  (let ((pos1 (position (first (set-difference points1 triangle)) points1))
		(pos2 (position (first (set-difference points2 triangle)) points2)))
	    (setf (nth pos1 (s3sx-sx3ids st1)) st2id)(setf (nth pos2 (s3sx-sx3ids st2)) st1id)))))))

(defun connect-s3simplices-within-list (sx1ids)
  (for (n 0 (1- (length sx1ids)))
       (for (m (1+ n) (1- (length sx1ids)))
	    (connect-s3simplices (nth n sx1ids) (nth m sx1ids)))))

(defun count-s3simplices-in-slice (ts)
  "counts the number of spacelike 3-simplices (i.e. spatial tetrahedra) in spatial slice ts"
  (let ((count 0))
    (maphash #'(lambda (id sx)
		 (declare (ignore id))
		 (when (= (s3sx-time sx) (mod ts NUM-T)) 
		   (incf count)))
	     *ID->SPATIAL-3SIMPLEX*)
    count))

(defun parse-s3simplex-parameters-line (line)
  "parses line for parameters but does nothing with them. The parameters are
stored in the 3-simplex data file for purposes of identification only"
  (with-input-from-string (s line)))

(defun parse-s3simplex-data-line (line)
  "parses line for spatial 3-simplex data"
  (with-input-from-string (s line)
    (let ((data (loop
		   :for num := (read s nil nil)
		   :while num
		   :collect num)))
      (make-s3simplex-v2 (nth 0 data) (nth 1 data) (nth 2 data) (nth 3 data) 
			 (nth 4 data) (nth 5 data) (nth 6 data) (nth 7 data)
			 (nth 8 data) (nth 9 data)))))

(defun load-s3simplex-data-from-file (infile)
  "loads spatial 3-simplex data from infile"
  (parse-s3simplex-parameters-line (read-line infile nil))
  (clrhash *ID->SPATIAL-3SIMPLEX*)
  (loop for line = (read-line infile nil)
     while line do (parse-s3simplex-data-line line)))

(defun save-s3simplex-data-to-file (outfile)
  "saves the spatial 3-simplex data to outfile"
  (format outfile "~A ~A ~A ~A ~A ~A ~A ~A ~A~%" BCTYPE STOPOLOGY NUM-T N-INIT 
	  N3-SL (hash-table-count *ID->SPATIAL-3SIMPLEX*) KAPPA-0 DELTA KAPPA-4)
  (maphash #'(lambda (k v)
	       (let ((pts (s3sx-points v))
		     (nbors (s3sx-sx3ids v)))
		 (format outfile "~A ~A ~A ~A ~A ~A ~A ~A ~A ~A~%"
			 k (s3sx-time v) (nth 0 pts) (nth 1 pts) (nth 2 pts) (nth 3 pts) 
			 (nth 0 nbors) (nth 1 nbors) (nth 2 nbors) (nth 3 nbors))))
	   *ID->SPATIAL-3SIMPLEX*))

(defun 4sx3p1->s3sx3p1 ()
  "4sx3p1->s3sx3p1 generates the spatial 3-simplex information for each 
spatial slice from the 4-simplex data for the entire spacetime."
  (clrhash *ID->SPATIAL-3SIMPLEX*)
  (setf *LAST-USED-S3SXID* 0)
  (for (ts 0 (1- NUM-T))
       (let ((41simplices (get-simplices-in-sandwich-of-type ts (1+ ts) 4)) ;; list of ids
	     (spatial-tetrahedra '()))
	 (dolist (41simplex 41simplices)
	   (push (make-s3simplex ts (4sx-lopts (get-4simplex 41simplex))) 
		 spatial-tetrahedra))
	 (connect-s3simplices-within-list spatial-tetrahedra))))

(defun generate-s3sx3p1-files (4sx3p1files)
  "4sx3p1files is a list of .4sx3p1 files. For each file in this list, this
function generates a .s3sx3p1 file. The prefix for the .4sx3p1 and the .s3sx3p1
file are identical"
  (loop for line = (read-line 4sx3p1files nil)
     while line do
       (let ((outfilename (change-file-suffix line "s3sx3p1")))
	 (with-open-file (indatafile line :direction :input)
	   (load-spacetime-from-file indatafile))
	 (4sx3p1->s3sx3p1)
	 (clrhash *ID->4SIMPLEX*)
	 (sb-ext:gc)
	 (with-open-file (outdatafile outfilename :direction :output)
	   (save-s3simplex-data-to-file outdatafile)))))

(defun add-subsimplices-to-stores (4sx)
  (let ((sxtype (4sx-type 4sx))
	(p0 (nth 0 (4sx-points 4sx)))
	(p1 (nth 1 (4sx-points 4sx)))
	(p2 (nth 2 (4sx-points 4sx)))
	(p3 (nth 3 (4sx-points 4sx)))
	(p4 (nth 4 (4sx-points 4sx))))
    ;; N0
    (setf (gethash p0 *N0STORE*) "")
    (setf (gethash p1 *N0STORE*) "")
    (setf (gethash p2 *N0STORE*) "")
    (setf (gethash p3 *N0STORE*) "")
    (setf (gethash p4 *N0STORE*) "")
    (cond ((= 1 sxtype) ;; (p0 | p1 p2 p3 p4)
	   ;; N1-SL
	   (setf (gethash (list p1 p2) *N1SLSTORE*) "")
	   (setf (gethash (list p1 p3) *N1SLSTORE*) "")
	   (setf (gethash (list p1 p4) *N1SLSTORE*) "")
	   (setf (gethash (list p2 p3) *N1SLSTORE*) "")
	   (setf (gethash (list p2 p4) *N1SLSTORE*) "")
	   (setf (gethash (list p3 p4) *N1SLSTORE*) "")
	   ;; N1-TL
	   (setf (gethash (list p0 p1) *N1TLSTORE*) "")
	   (setf (gethash (list p0 p2) *N1TLSTORE*) "")
	   (setf (gethash (list p0 p3) *N1TLSTORE*) "")
	   (setf (gethash (list p0 p4) *N1TLSTORE*) "")
	   ;; N2-SL
	   (setf (gethash (list p1 p2 p3) *N2SLSTORE*) "")
	   (setf (gethash (list p1 p2 p4) *N2SLSTORE*) "")
	   (setf (gethash (list p1 p3 p4) *N2SLSTORE*) "")
	   (setf (gethash (list p2 p3 p4) *N2SLSTORE*) "")
	   ;; N2-TL
	   (setf (gethash (list p0 p1 p2) *N2TLSTORE*) "")
	   (setf (gethash (list p0 p1 p3) *N2TLSTORE*) "")
	   (setf (gethash (list p0 p1 p4) *N2TLSTORE*) "")
	   (setf (gethash (list p0 p2 p3) *N2TLSTORE*) "")
	   (setf (gethash (list p0 p2 p4) *N2TLSTORE*) "")
	   (setf (gethash (list p0 p3 p4) *N2TLSTORE*) "")
	   ;; N3-SL
	   (setf (gethash (list p1 p2 p3 p4) *N3SLSTORE*) "")
	   ;; N3-TL-13
	   (setf (gethash (list p0 p1 p2 p3) *N3TL13STORE*) "")
	   (setf (gethash (list p0 p1 p2 p4) *N3TL13STORE*) "")
	   (setf (gethash (list p0 p1 p3 p4) *N3TL13STORE*) "")
	   (setf (gethash (list p0 p2 p3 p4) *N3TL13STORE*) "")
	   ;; N4-TL-14
	   (setf (gethash (list p0 p1 p2 p3 p4) *N4TL14STORE*) ""))
	  ((= 2 sxtype) ;; (p0 p1 | p2 p3 p4)
	   ;; N1-SL
	   (setf (gethash (list p0 p1) *N1SLSTORE*) "")
	   (setf (gethash (list p2 p3) *N1SLSTORE*) "")
	   (setf (gethash (list p2 p4) *N1SLSTORE*) "")
	   (setf (gethash (list p3 p4) *N1SLSTORE*) "")
	   ;; N1-TL
	   (setf (gethash (list p0 p2) *N1TLSTORE*) "")
	   (setf (gethash (list p0 p3) *N1TLSTORE*) "")
	   (setf (gethash (list p0 p4) *N1TLSTORE*) "")
	   (setf (gethash (list p1 p2) *N1TLSTORE*) "")
	   (setf (gethash (list p1 p3) *N1TLSTORE*) "")
	   (setf (gethash (list p1 p4) *N1TLSTORE*) "")
	   ;; N2-SL
	   (setf (gethash (list p2 p3 p4) *N2SLSTORE*) "")
	   ;; N2-TL
	   (setf (gethash (list p0 p1 p2) *N2TLSTORE*) "")
	   (setf (gethash (list p0 p1 p3) *N2TLSTORE*) "")
	   (setf (gethash (list p0 p1 p4) *N2TLSTORE*) "")
	   (setf (gethash (list p0 p2 p3) *N2TLSTORE*) "")
	   (setf (gethash (list p0 p2 p4) *N2TLSTORE*) "")
	   (setf (gethash (list p0 p3 p4) *N2TLSTORE*) "")
	   (setf (gethash (list p1 p2 p3) *N2TLSTORE*) "")
	   (setf (gethash (list p1 p2 p4) *N2TLSTORE*) "")
	   (setf (gethash (list p1 p3 p4) *N2TLSTORE*) "")
	   ;; N3-TL-13
	   (setf (gethash (list p0 p2 p3 p4) *N3TL13STORE*) "")
	   (setf (gethash (list p1 p2 p3 p4) *N3TL13STORE*) "")
	   ;; M3-TL-22
	   (setf (gethash (list p0 p1 p2 p3) *N3TL22STORE*) "")
	   (setf (gethash (list p0 p1 p2 p4) *N3TL22STORE*) "")
	   (setf (gethash (list p0 p1 p3 p4) *N3TL22STORE*) "")
	   ;; N4-TL-23
	   (setf (gethash (list p0 p1 p2 p3 p4) *N4TL23STORE*) ""))
	  ((= 3 sxtype) ;; (p0 p1 p2 | p3 p4)
	   ;; N1-SL
	   (setf (gethash (list p0 p1) *N1SLSTORE*) "")
	   (setf (gethash (list p0 p2) *N1SLSTORE*) "")
	   (setf (gethash (list p1 p2) *N1SLSTORE*) "")
	   (setf (gethash (list p3 p4) *N1SLSTORE*) "")
	   ;; N1-TL
	   (setf (gethash (list p0 p3) *N1TLSTORE*) "")
	   (setf (gethash (list p0 p4) *N1TLSTORE*) "")
	   (setf (gethash (list p1 p3) *N1TLSTORE*) "")
	   (setf (gethash (list p1 p4) *N1TLSTORE*) "")
	   (setf (gethash (list p2 p3) *N1TLSTORE*) "")
	   (setf (gethash (list p2 p4) *N1TLSTORE*) "")
	   ;; N2-SL
	   (setf (gethash (list p0 p1 p2) *N2SLSTORE*) "")
	   ;; N2-TL
	   (setf (gethash (list p0 p1 p3) *N2TLSTORE*) "")
	   (setf (gethash (list p0 p2 p3) *N2TLSTORE*) "")
	   (setf (gethash (list p1 p2 p3) *N2TLSTORE*) "")
	   (setf (gethash (list p0 p1 p4) *N2TLSTORE*) "")
	   (setf (gethash (list p0 p2 p4) *N2TLSTORE*) "")
	   (setf (gethash (list p1 p2 p4) *N2TLSTORE*) "")
	   (setf (gethash (list p0 p3 p4) *N2TLSTORE*) "")
	   (setf (gethash (list p1 p3 p4) *N2TLSTORE*) "")
	   (setf (gethash (list p2 p3 p4) *N2TLSTORE*) "")
	    ;; M3-TL-22
	   (setf (gethash (list p0 p1 p3 p4) *N3TL22STORE*) "")
	   (setf (gethash (list p0 p2 p3 p4) *N3TL22STORE*) "")
	   (setf (gethash (list p1 p2 p3 p4) *N3TL22STORE*) "")
	   ;; N3-TL-31
	   (setf (gethash (list p0 p1 p2 p3) *N3TL31STORE*) "")
	   (setf (gethash (list p0 p1 p2 p4) *N3TL31STORE*) "")
	   ;; N4-TL-32
	   (setf (gethash (list p0 p1 p2 p3 p4) *N4TL32STORE*) ""))
	  ((= 4 sxtype) ;; (p0 p1 p2 p3 | p4)
	   ;; N1-SL
	   (setf (gethash (list p0 p1) *N1SLSTORE*) "")
	   (setf (gethash (list p0 p2) *N1SLSTORE*) "")
	   (setf (gethash (list p0 p3) *N1SLSTORE*) "")
	   (setf (gethash (list p1 p2) *N1SLSTORE*) "")
	   (setf (gethash (list p1 p3) *N1SLSTORE*) "")
	   (setf (gethash (list p2 p3) *N1SLSTORE*) "")
	   ;; N1-TL
	   (setf (gethash (list p0 p4) *N1TLSTORE*) "")
	   (setf (gethash (list p1 p4) *N1TLSTORE*) "")
	   (setf (gethash (list p2 p4) *N1TLSTORE*) "")
	   (setf (gethash (list p3 p4) *N1TLSTORE*) "")
	   ;; N2-SL
	   (setf (gethash (list p0 p1 p2) *N2SLSTORE*) "")
	   (setf (gethash (list p0 p1 p3) *N2SLSTORE*) "")
	   (setf (gethash (list p0 p2 p3) *N2SLSTORE*) "")
	   (setf (gethash (list p1 p2 p3) *N2SLSTORE*) "")
	   ;; N2-TL
	   (setf (gethash (list p0 p1 p4) *N2TLSTORE*) "")
	   (setf (gethash (list p0 p2 p4) *N2TLSTORE*) "")
	   (setf (gethash (list p0 p3 p4) *N2TLSTORE*) "")
	   (setf (gethash (list p1 p2 p4) *N2TLSTORE*) "")
	   (setf (gethash (list p1 p3 p4) *N2TLSTORE*) "")
	   (setf (gethash (list p2 p3 p4) *N2TLSTORE*) "")
	   ;; N3-SL
	   (setf (gethash (list p0 p1 p2 p3) *N3SLSTORE*) "")
	   ;; N3-TL-31
	   (setf (gethash (list p0 p1 p2 p4) *N3TL31STORE*) "")
	   (setf (gethash (list p0 p1 p3 p4) *N3TL31STORE*) "")
	   (setf (gethash (list p0 p2 p3 p4) *N3TL31STORE*) "")
	   (setf (gethash (list p1 p2 p3 p4) *N3TL31STORE*) "")
	   ;; N4-TL-41
	   (setf (gethash (list p0 p1 p2 p3 p4) *N4TL41STORE*) "")))))

(defun count-and-display-bulk-variables ()
  (clrhash *N0STORE*)
  (clrhash *N1SLSTORE*)
  (clrhash *N1TLSTORE*)
  (clrhash *N2SLSTORE*)
  (clrhash *N2TLSTORE*)
  (clrhash *N3SLSTORE*)
  (clrhash *N3TL31STORE*)
  (clrhash *N3TL13STORE*)
  (clrhash *N3TL22STORE*)
  (clrhash *N4TL14STORE*)
  (clrhash *N4TL23STORE*)
  (clrhash *N4TL32STORE*)
  (clrhash *N4TL41STORE*)

  (maphash #'(lambda (id sx)
	       (declare (ignore id))
	       (add-subsimplices-to-stores sx))
	   *ID->4SIMPLEX*)
  (format t "N0 from store = ~A and bulk variable = ~A~%"
	  (hash-table-count *N0STORE*) N0)
  (format t "N1-SL from store = ~A and bulk variable = ~A~%"
	  (hash-table-count *N1SLSTORE*) N1-SL)
  (format t "N1-TL from store = ~A and bulk variable = ~A~%"
	  (hash-table-count *N1TLSTORE*) N1-TL)
  (format t "N2-SL from store = ~A and bulk variable = ~A~%"
	  (hash-table-count *N2SLSTORE*) N2-SL)
  (format t "N2-TL from store = ~A and bulk variable = ~A~%"
	  (hash-table-count *N2TLSTORE*) N2-TL)
  (format t "N3-SL from store = ~A and bulk variable = ~A~%"
	  (hash-table-count *N3SLSTORE*) N3-SL)
  (format t "N3-TL-13 + N3-TL-31 from store = ~A bulk variable = ~A~%"
	  (+ (hash-table-count *N3TL13STORE*) (hash-table-count *N3TL31STORE*))
	  N3-TL-31)
  (format t "N3-TL-22 from store = ~A and bulk variable = ~A~%"
	  (hash-table-count *N3TL22STORE*) N3-TL-22)
  (format t "N4-TL-14 + N4-TL-41 from store = ~A bulk variable = ~A~%"
	  (+ (hash-table-count *N4TL14STORE*) (hash-table-count *N4TL41STORE*))
	  N4-TL-41)
  (format t "N4-TL-23 + N4-TL-32 from store = ~A bulk variable = ~A~%"
	  (+ (hash-table-count *N4TL23STORE*) (hash-table-count *N4TL32STORE*))
	  N4-TL-32))