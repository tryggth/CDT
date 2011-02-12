;; 2simplex
;; (type tmlo tmhi (p0 p1 p2))
;; type = 0,1,2,3 tm[lo|hi] = [lo|hi] spatial time slice pj = points, nj = neighbor that does not have pj

;; 3simplex
;; (type tmlo tmhi (t0 t1 t2 t3) (n0 n1 n2 n3))
;; type = 1,2,3 tm[lo|hi] - [lo|hi] spatial time slice
;; tj = id of the 2sx that does not have pj
;; nj = id of the 3sx that does not have pj

;; for example given the following points
;; 
;; 9,10,11,12
;;------------- t = 2
;;
;; 5,6,7,8
;;------------- t = 1
;;
;; 1,2,3,4
;;------------- t = 0
;;
;;
;; the following 3 simplices 
;;
;; (1 | 5 6 7) (1 2 | 5 6) (1 2 | 6 7)
;;
;; will result in the following 2 simplex database
;;
;; [1] (1 0 1 (1 5 6))
;; [2] (1 0 1 (1 6 7))
;; [3] (1 0 1 (1 5 7))
;; [4] (0 0 1 (5 6 7))
;; [5] (2 0 1 (1 2 5))
;; [6] (2 0 1 (1 2 6))
;; [7] (1 0 1 (2 5 6))
;; [8] (2 0 1 (1 2 7))
;; [9] (1 0 1 (2 6 7))
;;
;; and the following 3 simplex database
;;
;; [1] (1 0 1 (4 2 3 1) (? 3 ? 2))
;; [2] (2 0 1 (7 1 6 5) (? 1 3 ?))
;; [3] (2 0 1 (9 2 8 6) (? 1 ? 2))


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
(defmacro get-2simplex (sxid)
  `(gethash ,sxid *ID->2SIMPLEX*))
(defmacro 2sx-points (2sx)
  `(fourth ,2sx))
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
	     (decf p1 (* 22 NUM-T)) (decf p2 (* 22 NUM-T)) (decf p3 (* 22 NUM-T)))
	    ((= 2 type)
	     (decf p2 (* 22 NUM-T)) (decf p3 (* 22 NUM-T)))
	    ((= 3 type)
	     (decf p3 (* 22 NUM-T)))))
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

;; simplex-data = ((typ tmlo tmhi (p0 p1 p2 p3)) (typ tmlo tmhi (p0 p1 p2 p3))...)
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
(defmacro 3sx-sx2ids (sx) `(sixth ,sx))
(defmacro get-3simplex (sxid)
  `(gethash ,sxid *ID->3SIMPLEX*))
(defun remove-3simplices (3sxids)
  (dolist (3sxid 3sxids)
    (let ((3sx (gethash 3sxid *ID->3SIMPLEX*)))
      (dolist (nid (3sx-sx3ids 3sx))
	(unless (= 0 nid)
	  (substitute 0 3sxid (3sx-sx3ids (get-3simplex nid))))))
    (remhash 3sxid *ID->3SIMPLEX*)))
(defun show-id->3simplex-store ()
  (maphash #'(lambda (3sxid 3sx) (format t "[~A] ~A~%" 3sxid 3sx)) *ID->3SIMPLEX*))

#|
(defun nth-point (3sx posn)
  (let ((others nil))
    (for (n (+ posn 1) (+ posn 3))
      (setf others (union others (fourth (gethash (nth (mod n 4) (3sx-sx2ids 3sx)) *ID->2SIMPLEX*)))))
    (first (set-difference others (fourth (gethash (nth posn (3sx-sx2ids 3sx)) *ID->2SIMPLEX*))))))
|#

(defmacro nth-point (sx n)
  `(nth ,n (3sx-points ,sx)))
(defmacro lo-points (sx)
  `(subseq (3sx-points ,sx) 0 (3sx-type ,sx)))
(defmacro hi-points (sx)
  `(subseq (3sx-points ,sx) (3sx-type ,sx)))

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

;; returns (1ids 2ids 3ids) where 1ids is a list of (1,3) ids in the sandwich etc.
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

;;(defun connect-simplices-in-sandwich (tlo thi)
;;(connect-3simplices-within-list (get-simplices-in-sandwich tlo thi)))

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


;;(defun save-spacetime-to-file (outfile)
;;  (format outfile "~A ~A ~A ~A ~A ~A ~A ~A ~A ~A ~A ~A ~A ~A ~A~%" 
;;	  BCTYPE STOPOLOGY NUM-T N-INIT *LAST-USED-POINT* 
;;	  N0 N1-SL N1-TL N2-SL N2-TL N3-TL-31 N3-TL-22 eps k0 k3)
;;  (maphash #'(lambda (k v)
;;	       (declare (ignore k))
;;	       (format outfile "~A ~A ~A ~A ~A ~A ~A~%" 
;;		       (3sx-type v) (3sx-tmlo v) (3sx-tmhi v)
;;		       (nth-point v 0) (nth-point v 1) (nth-point v 2) (nth-point v 3)))
;;	   *ID->3SIMPLEX*))

;; same as above, except also saves the neighboring 3simplex ids to files. This avoids reconnection
;; after loading the data from file.
(defun save-spacetime-to-file (outfile)
  (format outfile "~A ~A ~A ~A ~A ~A ~A ~A ~A ~A ~A ~A ~A ~A ~A ~A~%" 
	  BCTYPE STOPOLOGY NUM-T N-INIT *LAST-USED-POINT* *LAST-USED-3SXID* 
	  N0 N1-SL N1-TL N2-SL N2-TL N3-TL-31 N3-TL-22 eps k0 k3)
  (maphash #'(lambda (k v)
	       (format outfile "~A ~A ~A ~A ~A ~A ~A ~A ~A ~A ~A ~A~%" 
		       (3sx-type v) (3sx-tmlo v) (3sx-tmhi v)
		       (nth-point v 0) (nth-point v 1) (nth-point v 2) (nth-point v 3)
		       (nth 0 (3sx-sx3ids v)) (nth 1 (3sx-sx3ids v)) 
		       (nth 2 (3sx-sx3ids v)) (nth 3 (3sx-sx3ids v)) k))
	   *ID->3SIMPLEX*))

;;(defun parse-parameters-line (line)
;;  (with-input-from-string (s line)
;;  (let ((data (loop
;;		   :for num := (read s nil nil)
;;		   :while num
;;		   :collect num)))
;;    (setf BCTYPE (nth 0 data))
;;    (setf STOPOLOGY (nth 1 data))
;;    (setf NUM-T (nth 2 data))
;;    (setf N-INIT (nth 3 data))
;;    (setf *LAST-USED-POINT* (nth 4 data))
;;    (setf N0 (nth 5 data))
;;    (setf N1-SL (nth 6 data))
;;    (setf N1-TL (nth 7 data))
;;    (setf N2-SL (nth 8 data))
;;    (setf N2-TL (nth 9 data))
;;    (setf N3-TL-31 (nth 10 data))
;;    (setf N3-TL-22 (nth 11 data))
;;    (setf eps (nth 12 data))
;;    (setf k0 (nth 13 data))
;;    (setf k3 (nth 14 data)))))

(defun parse-parameters-line (line)
  (with-input-from-string (s line)
    (let ((data (loop
		   :for num := (read s nil nil)
		   :while num
		   :collect num)))
      (setf BCTYPE (nth 0 data))
      (setf STOPOLOGY (nth 1 data))
      (setf NUM-T (nth 2 data))
      (setf N-INIT (nth 3 data))
      (setf *LAST-USED-POINT* (nth 4 data))
      (setf *LAST-USED-3SXID* (nth 5 data))
      (setf N0 (nth 6 data))
      (setf N1-SL (nth 7 data))
      (setf N1-TL (nth 8 data))
      (setf N2-SL (nth 9 data))
      (setf N2-TL (nth 10 data))
      (setf N3-TL-31 (nth 11 data))
      (setf N3-TL-22 (nth 12 data))
      (setf eps (nth 13 data))
      (setf k0 (nth 14 data))
      (setf k3 (nth 15 data)))))

;;(defun parse-simplex-data-line (line)
;;  (with-input-from-string (s line)
;;  (let ((data (loop
;;		   :for num := (read s nil nil)
;;		   :while num
;;		   :collect num)))
;;    (make-3simplex (first data) (second data) (third data) 
;;		    (fourth data) (fifth data) (sixth data) (seventh data)))))

(defun parse-simplex-data-line (line)
  (with-input-from-string (s line)
    (let ((data (loop
		   :for num := (read s nil nil)
		   :while num
		   :collect num)))
      (make-3simplex-v4 (nth 0 data) (nth 1 data) (nth 2 data) (nth 3 data) 
			(nth 4 data) (nth 5 data) (nth 6 data) (nth 7 data)
			(nth 8 data) (nth 9 data) (nth 10 data) (nth 11 data)))))

;;(defun load-spacetime-from-file (infile)
;;  (parse-parameters-line (read-line infile nil))
;;  (loop for line = (read-line infile nil)
;;     while line do (parse-simplex-data-line line))
;;  (for (ts 0 (- NUM-T 1))
;;       (connect-simplices-in-sandwich ts (1+ ts) )
;;       (connect-simplices-in-adjacent-sandwiches ts (+ ts 1) (+ ts 2)))
;;  (when (and (string= STOPOLOGY "S2") (string= BCTYPE "PERIODIC"))
;;    (setf (symbol-function 'action) #'action-S1xS2)))

(defun load-spacetime-from-file (infile)
  (parse-parameters-line (read-line infile nil))
  (loop for line = (read-line infile nil)
     while line do (parse-simplex-data-line line))
  (when (and (string= STOPOLOGY "S2") (string= BCTYPE "PERIODIC"))
    (setf (symbol-function 'action) #'action-S1xS2)))

;; if the 3-simplices are connected, returns the id of the linking 2-simplex else returns 0
(defmacro link-id (sxid1 sxid2)
  `(let ((sx1 nil) (sx2 nil) (link nil))
     (if (and (setf sx1 (get-3simplex ,sxid1)) 
	      (setf sx2 (get-3simplex ,sxid2))
	      (setf link (intersection (3sx-sx2ids sx2) (3sx-sx2ids sx1))))
	 (first link)
	 0)))

;; if the 3simplices are connected returns (id pos1 pos2) where id is the id of the linking 2-simplex
;; pos1 is the position of this linkid in the sx2ids list of the first simplex
;; pos2 is the position of this linkid in the sx2ids list of the second simplex
;; if the two simplices are not connected, returns (0 -1 -1)
;;(defmacro link-id-and-positions (sxid1 sxid2)
;;  `(let ((sx1 nil) (sx2 nil) (link nil))
;;   (if (and (setf sx1 (get-3simplex ,sxid1))
;;	      (setf sx2 (get-3simplex ,sxid2))
;;	      (setf link (intersection (3sx-sx2ids sx1) (3sx-sx2ids sx2))))
;;	 (list (first link) (position (first link) (3sx-sx2ids sx1)) (position (first link) (3sx-sx2ids sx2)))
;;	 (list 0 -1 -1))))
  
(defun neighbors-of-type (sx type)
  (let ((nbors nil)
	(nsx nil)
	(nids (3sx-sx3ids sx)))
    (for (n 0 3)
      (when (and (setf nsx (get-3simplex (nth n nids))) (= type (3sx-type nsx)))
	(pushnew (nth n nids) nbors)))
    nbors))
