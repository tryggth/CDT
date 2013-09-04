(declaim (optimize (speed 3)
		   (compilation-speed 0)
		   (debug 0)
		   (safety 0)))
;; sl 1-simplex is (tslice (p0 p1))
;; tl 1-simplex is (type tmlo (p0 p1)) where type = 1
;; sl 2-simplex is (tslice (p0 p1 p2))
;; tl 2-simplex is (type tmlo (p0 p1 p2)) where type = 1,2
;; sl 3-simplex is (tslice (p0 p1 p2 p3))
;; tl 3-simplex is (type tmlo (p0 p1 p2 p3)) where type = 1,2,3
;; tl 4-simplex is (type tmlo tmhi (p0 p1 p2 p3 p4) (n0 n1 n2 n3 n4))
;; where type = 1,2,3,4 
;; nj = id of the 4sx that does not have pj


(defun make-3simplices (type tmlo tmhi p0 p1 p2 p3 p4)
  "makes all 3-simplices with the specified data iff they don't already exist"
  (cond ((= 1 type)
	 (setf (gethash `(1 ,tmlo (,p0 ,p1 ,p2 ,p3)) *TL3SIMPLEX->ID*) 0)
	 (setf (gethash `(1 ,tmlo (,p0 ,p1 ,p2 ,p4)) *TL3SIMPLEX->ID*) 0)
	 (setf (gethash `(1 ,tmlo (,p0 ,p1 ,p3 ,p4)) *TL3SIMPLEX->ID*) 0)
	 (setf (gethash `(1 ,tmlo (,p0 ,p2 ,p3 ,p4)) *TL3SIMPLEX->ID*) 0)
	 (setf (gethash `(,tmhi (,p1 ,p2 ,p3 ,p4)) *SL3SIMPLEX->ID*) 0))
	((= 2 type)
	 (setf (gethash `(2 ,tmlo (,p0 ,p1 ,p2 ,p3)) *TL3SIMPLEX->ID*) 0)
	 (setf (gethash `(2 ,tmlo (,p0 ,p1 ,p2 ,p4)) *TL3SIMPLEX->ID*) 0)
	 (setf (gethash `(2 ,tmlo (,p0 ,p1 ,p3 ,p4)) *TL3SIMPLEX->ID*) 0)
	 (setf (gethash `(1 ,tmlo (,p0 ,p2 ,p3 ,p4)) *TL3SIMPLEX->ID*) 0)
	 (setf (gethash `(1 ,tmlo (,p1 ,p2 ,p3 ,p4)) *TL3SIMPLEX->ID*) 0))
	((= 3 type)
	 (setf (gethash `(3 ,tmlo (,p0 ,p1 ,p2 ,p3)) *TL3SIMPLEX->ID*) 0)
	 (setf (gethash `(3 ,tmlo (,p0 ,p1 ,p2 ,p4)) *TL3SIMPLEX->ID*) 0)
	 (setf (gethash `(2 ,tmlo (,p0 ,p1 ,p3 ,p4)) *TL3SIMPLEX->ID*) 0)
	 (setf (gethash `(2 ,tmlo (,p0 ,p2 ,p3 ,p4)) *TL3SIMPLEX->ID*) 0)
	 (setf (gethash `(2 ,tmlo (,p1 ,p2 ,p3 ,p4)) *TL3SIMPLEX->ID*) 0))
	((= 4 type)
	 (setf (gethash `(,tmlo (,p0 ,p1 ,p2 ,p3)) *SL3SIMPLEX->ID*) 0)
	 (setf (gethash `(3 ,tmlo (,p0 ,p1 ,p2 ,p4)) *TL3SIMPLEX->ID*) 0)
	 (setf (gethash `(3 ,tmlo (,p0 ,p1 ,p3 ,p4)) *TL3SIMPLEX->ID*) 0)
	 (setf (gethash `(3 ,tmlo (,p0 ,p2 ,p3 ,p4)) *TL3SIMPLEX->ID*) 0)
	 (setf (gethash `(3 ,tmlo (,p1 ,p2 ,p3 ,p4)) *TL3SIMPLEX->ID*) 0))))


(defun make-2simplices (type tmlo tmhi p0 p1 p2 p3 p4)
  "makes all 2-simplices with the specified data iff they don't already exist"
  (cond ((= 1 type)
	 (setf (gethash `(1 ,tmlo (,p0 ,p1 ,p2)) *TL2SIMPLEX->ID*) 0)
	 (setf (gethash `(1 ,tmlo (,p0 ,p1 ,p3)) *TL2SIMPLEX->ID*) 0)
	 (setf (gethash `(1 ,tmlo (,p0 ,p1 ,p4)) *TL2SIMPLEX->ID*) 0)
	 (setf (gethash `(1 ,tmlo (,p0 ,p2 ,p3)) *TL2SIMPLEX->ID*) 0)
	 (setf (gethash `(1 ,tmlo (,p0 ,p2 ,p4)) *TL2SIMPLEX->ID*) 0)
	 (setf (gethash `(1 ,tmlo (,p0 ,p3 ,p4)) *TL2SIMPLEX->ID*) 0)
	 (setf (gethash `(,tmhi (,p1 ,p2 ,p3)) *SL2SIMPLEX->ID*) 0)
	 (setf (gethash `(,tmhi (,p1 ,p2 ,p4)) *SL2SIMPLEX->ID*) 0)
	 (setf (gethash `(,tmhi (,p1 ,p3 ,p4)) *SL2SIMPLEX->ID*) 0)
	 (setf (gethash `(,tmhi (,p2 ,p3 ,p4)) *SL2SIMPLEX->ID*) 0))
	((= 2 type)
	 (setf (gethash `(2 ,tmlo (,p0 ,p1 ,p2)) *TL2SIMPLEX->ID*) 0)
	 (setf (gethash `(2 ,tmlo (,p0 ,p1 ,p3)) *TL2SIMPLEX->ID*) 0)
	 (setf (gethash `(2 ,tmlo (,p0 ,p1 ,p4)) *TL2SIMPLEX->ID*) 0)
	 (setf (gethash `(1 ,tmlo (,p0 ,p2 ,p3)) *TL2SIMPLEX->ID*) 0)
	 (setf (gethash `(1 ,tmlo (,p0 ,p2 ,p4)) *TL2SIMPLEX->ID*) 0)
	 (setf (gethash `(1 ,tmlo (,p0 ,p3 ,p4)) *TL2SIMPLEX->ID*) 0)
	 (setf (gethash `(1 ,tmlo (,p1 ,p2 ,p3)) *TL2SIMPLEX->ID*) 0)
	 (setf (gethash `(1 ,tmlo (,p1 ,p2 ,p4)) *TL2SIMPLEX->ID*) 0)
	 (setf (gethash `(1 ,tmlo (,p1 ,p3 ,p4)) *TL2SIMPLEX->ID*) 0)
	 (setf (gethash `(,tmhi (,p2 ,p3 ,p4)) *SL2SIMPLEX->ID*) 0))
	((= 3 type)
	 (setf (gethash `(,tmlo (,p0 ,p1 ,p2)) *SL2SIMPLEX->ID*) 0)
	 (setf (gethash `(2 ,tmlo (,p0 ,p1 ,p3)) *TL2SIMPLEX->ID*) 0)
	 (setf (gethash `(2 ,tmlo (,p0 ,p1 ,p4)) *TL2SIMPLEX->ID*) 0)
	 (setf (gethash `(2 ,tmlo (,p0 ,p2 ,p3)) *TL2SIMPLEX->ID*) 0)
	 (setf (gethash `(2 ,tmlo (,p0 ,p2 ,p4)) *TL2SIMPLEX->ID*) 0)
	 (setf (gethash `(1 ,tmlo (,p0 ,p3 ,p4)) *TL2SIMPLEX->ID*) 0)
	 (setf (gethash `(2 ,tmlo (,p1 ,p2 ,p3)) *TL2SIMPLEX->ID*) 0)
	 (setf (gethash `(2 ,tmlo (,p1 ,p2 ,p4)) *TL2SIMPLEX->ID*) 0)
	 (setf (gethash `(1 ,tmlo (,p1 ,p3 ,p4)) *TL2SIMPLEX->ID*) 0)
	 (setf (gethash `(1 ,tmlo (,p2 ,p3 ,p4)) *TL2SIMPLEX->ID*) 0))
	((= 4 type)
	 (setf (gethash `(,tmlo (,p0 ,p1 ,p2)) *SL2SIMPLEX->ID*) 0)
	 (setf (gethash `(,tmlo (,p0 ,p1 ,p3)) *SL2SIMPLEX->ID*) 0)
	 (setf (gethash `(2 ,tmlo (,p0 ,p1 ,p4)) *TL2SIMPLEX->ID*) 0)
	 (setf (gethash `(,tmlo (,p0 ,p2 ,p3)) *SL2SIMPLEX->ID*) 0)
	 (setf (gethash `(2 ,tmlo (,p0 ,p2 ,p4)) *TL2SIMPLEX->ID*) 0)
	 (setf (gethash `(2 ,tmlo (,p0 ,p3 ,p4)) *TL2SIMPLEX->ID*) 0)
	 (setf (gethash `(,tmlo (,p1 ,p2 ,p3)) *SL2SIMPLEX->ID*) 0)
	 (setf (gethash `(2 ,tmlo (,p1 ,p2 ,p4)) *TL2SIMPLEX->ID*) 0)
	 (setf (gethash `(2 ,tmlo (,p1 ,p3 ,p4)) *TL2SIMPLEX->ID*) 0)
	 (setf (gethash `(2 ,tmlo (,p2 ,p3 ,p4)) *TL2SIMPLEX->ID*) 0))))


(defun make-1simplices (type tmlo tmhi p0 p1 p2 p3 p4)
  "makes all 1-simplices with the specified data iff they don't already exist"
  (cond ((= 1 type) ;; p0 | p1 p2 p3 p4
	 (setf (gethash `(1 ,tmlo (,p0 ,p1)) *TL1SIMPLEX->ID*) 0)
	 (setf (gethash `(1 ,tmlo (,p0 ,p2)) *TL1SIMPLEX->ID*) 0)
	 (setf (gethash `(1 ,tmlo (,p0 ,p3)) *TL1SIMPLEX->ID*) 0)
	 (setf (gethash `(1 ,tmlo (,p0 ,p4)) *TL1SIMPLEX->ID*) 0)
	 (setf (gethash `(,tmhi (,p1 ,p2)) *SL1SIMPLEX->ID*) 0)
	 (setf (gethash `(,tmhi (,p1 ,p3)) *SL1SIMPLEX->ID*) 0)
	 (setf (gethash `(,tmhi (,p1 ,p4)) *SL1SIMPLEX->ID*) 0)
	 (setf (gethash `(,tmhi (,p2 ,p3)) *SL1SIMPLEX->ID*) 0)
	 (setf (gethash `(,tmhi (,p2 ,p4)) *SL1SIMPLEX->ID*) 0)
	 (setf (gethash `(,tmhi (,p3 ,p4)) *SL1SIMPLEX->ID*) 0))
	((= 2 type) ;; p0 p1 | p2 p3 p4
	 (setf (gethash `(,tmlo (,p0 ,p1)) *SL1SIMPLEX->ID*) 0)
	 (setf (gethash `(1 ,tmlo (,p0 ,p2)) *TL1SIMPLEX->ID*) 0)
	 (setf (gethash `(1 ,tmlo (,p0 ,p3)) *TL1SIMPLEX->ID*) 0)
	 (setf (gethash `(1 ,tmlo (,p0 ,p4)) *TL1SIMPLEX->ID*) 0)
	 (setf (gethash `(1 ,tmlo (,p1 ,p2)) *TL1SIMPLEX->ID*) 0)
	 (setf (gethash `(1 ,tmlo (,p1 ,p3)) *TL1SIMPLEX->ID*) 0)
	 (setf (gethash `(1 ,tmlo (,p1 ,p4)) *TL1SIMPLEX->ID*) 0)
	 (setf (gethash `(,tmhi (,p2 ,p3)) *SL1SIMPLEX->ID*) 0)
	 (setf (gethash `(,tmhi (,p2 ,p4)) *SL1SIMPLEX->ID*) 0)
	 (setf (gethash `(,tmhi (,p3 ,p4)) *SL1SIMPLEX->ID*) 0))
	((= 3 type) ;; p0 p1 p2 | p3 p4
	 (setf (gethash `(,tmlo (,p0 ,p1)) *SL1SIMPLEX->ID*) 0)
	 (setf (gethash `(,tmlo (,p0 ,p2)) *SL1SIMPLEX->ID*) 0)
	 (setf (gethash `(1 ,tmlo (,p0 ,p3)) *TL1SIMPLEX->ID*) 0)
	 (setf (gethash `(1 ,tmlo (,p0 ,p4)) *TL1SIMPLEX->ID*) 0)
	 (setf (gethash `(,tmlo (,p1 ,p2)) *SL1SIMPLEX->ID*) 0)
	 (setf (gethash `(1 ,tmlo (,p1 ,p3)) *TL1SIMPLEX->ID*) 0)
	 (setf (gethash `(1 ,tmlo (,p1 ,p4)) *TL1SIMPLEX->ID*) 0)
	 (setf (gethash `(1 ,tmlo (,p2 ,p3)) *TL1SIMPLEX->ID*) 0)
	 (setf (gethash `(1 ,tmlo (,p2 ,p4)) *TL1SIMPLEX->ID*) 0)
	 (setf (gethash `(,tmhi (,p3 ,p4)) *SL1SIMPLEX->ID*) 0))
	((= 4 type) ;; p0 p1 p2 p3 | p4
	 (setf (gethash `(,tmlo (,p0 ,p1)) *SL1SIMPLEX->ID*) 0)
	 (setf (gethash `(,tmlo (,p0 ,p2)) *SL1SIMPLEX->ID*) 0)
	 (setf (gethash `(,tmlo (,p0 ,p3)) *SL1SIMPLEX->ID*) 0)
	 (setf (gethash `(1 ,tmlo (,p0 ,p4)) *TL1SIMPLEX->ID*) 0)
	 (setf (gethash `(,tmlo (,p1 ,p2)) *SL1SIMPLEX->ID*) 0)
	 (setf (gethash `(,tmlo (,p1 ,p3)) *SL1SIMPLEX->ID*) 0)
	 (setf (gethash `(1 ,tmlo (,p1 ,p4)) *TL1SIMPLEX->ID*) 0)
	 (setf (gethash `(,tmlo (,p2 ,p3)) *SL1SIMPLEX->ID*) 0)
	 (setf (gethash `(1 ,tmlo (,p2 ,p4)) *TL1SIMPLEX->ID*) 0)
	 (setf (gethash `(1 ,tmlo (,p3 ,p4)) *TL1SIMPLEX->ID*) 0))))


(defun remove-tl3simplex (tl3sx)
  (remhash tl3sx *TL3SIMPLEX->ID*))
(defun remove-sl3simplex (sl3sx)
  (remhash sl3sx *SL3SIMPLEX->ID*))
(defun remove-tl3simplices (tl3sxs)
  (dolist (tl3sx tl3sxs)
    (remove-tl3simplex tl3sx)))
(defun remove-sl3simplices (sl3sxs)
  (dolist (sl3sx sl3sxs)
    (remove-sl3simplex sl3sx)))


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


(defun make-4simplex (type tmlo tmhi p0 p1 p2 p3 p4)
  (let ((sx4id (next-4simplex-id)))
    (setf (gethash sx4id *ID->4SIMPLEX*)
	  (list type tmlo tmhi (list p0 p1 p2 p3 p4) (list 0 0 0 0 0)))
    (make-3simplices type tmlo tmhi p0 p1 p2 p3 p4)
    (make-2simplices type tmlo tmhi p0 p1 p2 p3 p4)
    (make-1simplices type tmlo tmhi p0 p1 p2 p3 p4)
    sx4id))

;;make fixed simplex
(defun make-4simplex-f (type tmlo tmhi p0 p1 p2 p3 p4)
  (let ((sx4id (next-4simplex-id)))
    (push sx4id *FIXED-4SXID*)
    (setf (gethash sx4id *ID->4SIMPLEX*)
	  (list type tmlo tmhi (list p0 p1 p2 p3 p4) (list 0 0 0 0 0)))
    (make-3simplices type tmlo tmhi p0 p1 p2 p3 p4)
    (make-2simplices type tmlo tmhi p0 p1 p2 p3 p4)
    (make-1simplices type tmlo tmhi p0 p1 p2 p3 p4)
    sx4id))

;; same as make-4simplex except the points are packed in a list called pts
(defun make-4simplex-v2 (type tmlo tmhi pts)
  (let ((p0 (nth 0 pts))
	(p1 (nth 1 pts))
	(p2 (nth 2 pts))
	(p3 (nth 3 pts))
	(p4 (nth 4 pts)))
    (make-4simplex type tmlo tmhi p0 p1 p2 p3 p4)))


;; this version is used only during initialization. If periodic b.c. are 
;; specified, it adjusts the points on the final time slice, since the t=T 
;; slice is identified with t=0 slice.
(defun make-4simplex-v3 (type tmlo tmhitmp p0tmp p1tmp p2tmp p3tmp p4tmp)
  (let ((p0 p0tmp) (p1 p1tmp) (p2 p2tmp) (p3 p3tmp) (p4 p4tmp) (tmhi tmhitmp))
    (when (and (string= BCTYPE "PERIODIC") (= NUM-T tmhitmp))
      (setf tmhi 0)
      (cond ((= 1 type)
	     (decf p1 (* 25 NUM-T)) (decf p2 (* 25 NUM-T)) (decf p3 (* 
25 
NUM-T)) 
	     (decf p4 (* 25 NUM-T)))
	    ((= 2 type)
	     (decf p2 (* 25 NUM-T)) (decf p3 (* 25 NUM-T)) (decf p4 (* 
25 
NUM-T)))
	    ((= 3 type)
	     (decf p3 (* 25 NUM-T)) (decf p4 (* 25 NUM-T)))
	    ((= 4 type)
	     (decf p4 (* 25 NUM-T)))))
    (make-4simplex type tmlo tmhi p0 p1 p2 p3 p4)))

;;make fixed 4simplex.
(defun make-4simplex-v3f (type tmlo tmhitmp p0tmp p1tmp p2tmp p3tmp p4tmp)
  (let ((p0 p0tmp) (p1 p1tmp) (p2 p2tmp) (p3 p3tmp) (p4 p4tmp) (tmhi tmhitmp))
    (when (and (string= BCTYPE "PERIODIC") (= NUM-T tmhitmp))
      (setf tmhi 0)
      (cond ((= 1 type)
	     (decf p1 (* 25 NUM-T)) (decf p2 (* 25 NUM-T)) (decf p3 (* 
25 
NUM-T)) 
	     (decf p4 (* 25 NUM-T)))
	    ((= 2 type)
	     (decf p2 (* 25 NUM-T)) (decf p3 (* 25 NUM-T)) (decf p4 (* 
25 
NUM-T)))
	    ((= 3 type)
	     (decf p3 (* 25 NUM-T)) (decf p4 (* 25 NUM-T)))
	    ((= 4 type)
	     (decf p4 (* 25 NUM-T)))))
    (make-4simplex-f type tmlo tmhi p0 p1 p2 p3 p4)))


;; this version is used for loading the simplex data from file
(defun make-4simplex-v4 (type tmlo tmhi p0 p1 p2 p3 p4 n0 n1 n2 n3 n4 sx4id)
    (setf (gethash sx4id *ID->4SIMPLEX*)
	  (list type tmlo tmhi (list p0 p1 p2 p3 p4) (list n0 n1 n2 n3 n4)))
    (make-3simplices type tmlo tmhi p0 p1 p2 p3 p4)
    (make-2simplices type tmlo tmhi p0 p1 p2 p3 p4)
    (make-1simplices type tmlo tmhi p0 p1 p2 p3 p4))


;; all the simplex data, not just the points, is packed into a list
(defun make-4simplex-v5 (simplex-data)
  (let ((type (first simplex-data))
	(tmlo (second simplex-data))
	(tmhi (third simplex-data))
	(pts (fourth simplex-data)))
    (make-4simplex-v2 type tmlo tmhi pts)))


;; simplex-data-list = ((typ tmlo tmhi (p0 p1 p2 p3 p4))...)
;; the ids of the simplices are returned
(defun make-4simplices-in-bulk (simplex-data-list)
  (let ((4sxids nil))
    (dolist (simplex-data simplex-data-list)
      (push (make-4simplex-v5 simplex-data) 4sxids))
    4sxids))


;;(type tmlo tmhi (p0 p1 p2 p3 p4) (n0 n1 n2 n3 n4))
(defmacro 4sx-type (sx) `(first ,sx))
(defmacro 4sx-tmlo (sx) `(second ,sx))
(defmacro 4sx-tmhi (sx) `(third ,sx))
(defmacro 4sx-points (sx) `(fourth ,sx))
(defmacro 4sx-sx4ids (sx) `(fifth ,sx))
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


(defun show-tl3simplex-store ()
  (let ((count 1))
    (maphash #'(lambda (tl3sx id)
		 (format t "[~A] ~A ~A~%" count tl3sx id)
		 (incf count))
	     *TL3SIMPLEX->ID*)))
(defun show-sl3simplex-store ()
  (let ((count 1))
    (maphash #'(lambda (sl3sx id)
		 (format t "[~A] ~A ~A~%" count sl3sx id)
		 (incf count))
	     *SL3SIMPLEX->ID*)))


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


(defun count-n3TL22 ()
  "counts the number of (2,2) timelike 3-simplices"
  (let ((count 0))
    (maphash #'(lambda (k v)
		 (declare (ignore k))
		 (if (= 2 (first v))
		     (incf count)))
	     *TL3SIMPLEX->ID*)
    count))


(defun count-n3TL31 ()
  "counts the number of (1,3) + (3,1) timelike 3-simplices"
  (let ((count 0))
    (maphash #'(lambda (k v)
		 (declare (ignore k))
		 (if (or (= 1 (first v))
			 (= 3 (first v)))
		     (incf count)))
	     *TL3SIMPLEX->ID*)
    count))


(defmacro nth-point (sx n)
  `(nth ,n (4sx-points ,sx)))
(defmacro nth-neighbor (sx n)
  `(nth ,n (4sx-sx4ids ,sx)))


(defun connect-4simplices (sx1id sx2id)
  (let ((sx1 nil) (sx2 nil))
    (when (and (setf sx1 (get-4simplex sx1id))
	       (setf sx2 (get-4simplex sx2id)))
      (let ((tetra (intersection (4sx-points sx1) (4sx-points sx2))))
	(when (= 4 (length tetra))
	  (let* ((pts1 (4sx-points sx1))
		 (pts2 (4sx-points sx2))
		 (pos1 (position (first (set-difference pts1 tetra)) pts1))
		 (pos2 (position (first (set-difference pts2 tetra)) pts2)))
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
		       (4sx-type v) (4sx-tmlo v) (4sx-tmhi v)
		       (nth-point v 0) (nth-point v 1) (nth-point v 2) 
		       (nth-point v 3) (nth-point v 4)
		       (nth-neighbor v 0) (nth-neighbor v 1)
		       (nth-neighbor v 2) (nth-neighbor v 3)
		       (nth-neighbor v 4) k))
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


#|
(defun 4simplex->vertices (4sx)
  "populates the *VERTEX-GRAPH* hashtable by taking apart the 4sx"
  (let ((p0 (nth-point 4sx 0))
	(p1 (nth-point 4sx 1))
	(p2 (nth-point 4sx 2))
	(p3 (nth-point 4sx 3))
	(p4 (nth-point 4sx 4))
	(typ (4sx-type 4sx))
	(tml (4sx-tmlo 4sx))
	(tmh (4sx-tmhi 4sx)))
    (cond ((= 1 typ) ; (p0 | p1 p2 p3 p4)
	   (pushnew `(,p1 -1) (gethash `(,p0 ,tml) *VERTEX-GRAPH*) :test 'equal )
	   (pushnew `(,p2 -1) (gethash `(,p0 ,tml) *VERTEX-GRAPH*) :test 'equal )
	   (pushnew `(,p3 -1) (gethash `(,p0 ,tml) *VERTEX-GRAPH*) :test 'equal )
	   (pushnew `(,p4 -1) (gethash `(,p0 ,tml) *VERTEX-GRAPH*) :test 'equal )
	   (pushnew `(,p0 -1) (gethash `(,p1 ,tmh) *VERTEX-GRAPH*) :test 'equal )
	   (pushnew `(,p2 1) (gethash `(,p1 ,tmh) *VERTEX-GRAPH*) :test 'equal )
	   (pushnew `(,p3 1) (gethash `(,p1 ,tmh) *VERTEX-GRAPH*) :test 'equal )
	   (pushnew `(,p4 1) (gethash `(,p1 ,tmh) *VERTEX-GRAPH*) :test 'equal )
	   (pushnew `(,p0 -1) (gethash `(,p2 ,tmh) *VERTEX-GRAPH*) :test 'equal )
	   (pushnew `(,p1 1) (gethash `(,p2 ,tmh) *VERTEX-GRAPH*) :test 'equal )
	   (pushnew `(,p3 1) (gethash `(,p2 ,tmh) *VERTEX-GRAPH*) :test 'equal )
	   (pushnew `(,p4 1) (gethash `(,p2 ,tmh) *VERTEX-GRAPH*) :test 'equal )
	   (pushnew `(,p0 -1) (gethash `(,p3 ,tmh) *VERTEX-GRAPH*) :test 'equal )
	   (pushnew `(,p1 1) (gethash `(,p3 ,tmh) *VERTEX-GRAPH*) :test 'equal )
	   (pushnew `(,p2 1) (gethash `(,p3 ,tmh) *VERTEX-GRAPH*) :test 'equal )
	   (pushnew `(,p4 1) (gethash `(,p3 ,tmh) *VERTEX-GRAPH*) :test 'equal )
	   (pushnew `(,p0 -1) (gethash `(,p4 ,tmh) *VERTEX-GRAPH*) :test 'equal )
	   (pushnew `(,p1 1) (gethash `(,p4 ,tmh) *VERTEX-GRAPH*) :test 'equal )
	   (pushnew `(,p2 1) (gethash `(,p4 ,tmh) *VERTEX-GRAPH*) :test 'equal )
	   (pushnew `(,p3 1) (gethash `(,p4 ,tmh) *VERTEX-GRAPH*) :test 'equal ))
	  ((= 2 typ) ; (p0 p1 | p2 p3 p4)
	   (pushnew `(,p1 1) (gethash `(,p0 ,tml) *VERTEX-GRAPH*) :test 'equal )
	   (pushnew `(,p2 -1) (gethash `(,p0 ,tml) *VERTEX-GRAPH*) :test 'equal )
	   (pushnew `(,p3 -1) (gethash `(,p0 ,tml) *VERTEX-GRAPH*) :test 'equal )
	   (pushnew `(,p4 -1) (gethash `(,p0 ,tml) *VERTEX-GRAPH*) :test 'equal )
	   (pushnew `(,p0 1) (gethash `(,p1 ,tmh) *VERTEX-GRAPH*) :test 'equal )
	   (pushnew `(,p2 -1) (gethash `(,p1 ,tmh) *VERTEX-GRAPH*) :test 'equal )
	   (pushnew `(,p3 -1) (gethash `(,p1 ,tmh) *VERTEX-GRAPH*) :test 'equal )
	   (pushnew `(,p4 -1) (gethash `(,p1 ,tmh) *VERTEX-GRAPH*) :test 'equal )
	   (pushnew `(,p0 -1) (gethash `(,p2 ,tmh) *VERTEX-GRAPH*) :test 'equal )
	   (pushnew `(,p1 -1) (gethash `(,p2 ,tmh) *VERTEX-GRAPH*) :test 'equal )
	   (pushnew `(,p3 1) (gethash `(,p2 ,tmh) *VERTEX-GRAPH*) :test 'equal )
	   (pushnew `(,p4 1) (gethash `(,p2 ,tmh) *VERTEX-GRAPH*) :test 'equal )
	   (pushnew `(,p0 -1) (gethash `(,p3 ,tmh) *VERTEX-GRAPH*) :test 'equal )
	   (pushnew `(,p1 -1) (gethash `(,p3 ,tmh) *VERTEX-GRAPH*) :test 'equal )
	   (pushnew `(,p2 1) (gethash `(,p3 ,tmh) *VERTEX-GRAPH*) :test 'equal )
	   (pushnew `(,p4 1) (gethash `(,p3 ,tmh) *VERTEX-GRAPH*) :test 'equal )
	   (pushnew `(,p0 -1) (gethash `(,p4 ,tmh) *VERTEX-GRAPH*) :test 'equal )
	   (pushnew `(,p1 -1) (gethash `(,p4 ,tmh) *VERTEX-GRAPH*) :test 'equal )
	   (pushnew `(,p2 1) (gethash `(,p4 ,tmh) *VERTEX-GRAPH*) :test 'equal )
	   (pushnew `(,p3 1) (gethash `(,p4 ,tmh) *VERTEX-GRAPH*) :test 'equal ))
	  ((= 3 typ) ; (p0 p1 p2 | p3 p4)
	   (pushnew `(,p1 1) (gethash `(,p0 ,tml) *VERTEX-GRAPH*) :test 'equal )
	   (pushnew `(,p2 1) (gethash `(,p0 ,tml) *VERTEX-GRAPH*) :test 'equal )
	   (pushnew `(,p3 -1) (gethash `(,p0 ,tml) *VERTEX-GRAPH*) :test 'equal )
	   (pushnew `(,p4 -1) (gethash `(,p0 ,tml) *VERTEX-GRAPH*) :test 'equal )
	   (pushnew `(,p0 1) (gethash `(,p1 ,tmh) *VERTEX-GRAPH*) :test 'equal )
	   (pushnew `(,p2 1) (gethash `(,p1 ,tmh) *VERTEX-GRAPH*) :test 'equal )
	   (pushnew `(,p3 -1) (gethash `(,p1 ,tmh) *VERTEX-GRAPH*) :test 'equal )
	   (pushnew `(,p4 -1) (gethash `(,p1 ,tmh) *VERTEX-GRAPH*) :test 'equal )
	   (pushnew `(,p0 1) (gethash `(,p2 ,tmh) *VERTEX-GRAPH*) :test 'equal )
	   (pushnew `(,p1 1) (gethash `(,p2 ,tmh) *VERTEX-GRAPH*) :test 'equal )
	   (pushnew `(,p3 -1) (gethash `(,p2 ,tmh) *VERTEX-GRAPH*) :test 'equal )
	   (pushnew `(,p4 -1) (gethash `(,p2 ,tmh) *VERTEX-GRAPH*) :test 'equal )
	   (pushnew `(,p0 -1) (gethash `(,p3 ,tmh) *VERTEX-GRAPH*) :test 'equal )
	   (pushnew `(,p1 -1) (gethash `(,p3 ,tmh) *VERTEX-GRAPH*) :test 'equal )
	   (pushnew `(,p2 -1) (gethash `(,p3 ,tmh) *VERTEX-GRAPH*) :test 'equal )
	   (pushnew `(,p4 1) (gethash `(,p3 ,tmh) *VERTEX-GRAPH*) :test 'equal )
	   (pushnew `(,p0 -1) (gethash `(,p4 ,tmh) *VERTEX-GRAPH*) :test 'equal )
	   (pushnew `(,p1 -1) (gethash `(,p4 ,tmh) *VERTEX-GRAPH*) :test 'equal )
	   (pushnew `(,p2 -1) (gethash `(,p4 ,tmh) *VERTEX-GRAPH*) :test 'equal )
	   (pushnew `(,p3 1) (gethash `(,p4 ,tmh) *VERTEX-GRAPH*) :test 'equal ))
	  ((= 4 typ) ; (p0 p1 p2 p3 | p4)
	   (pushnew `(,p1 1) (gethash `(,p0 ,tml) *VERTEX-GRAPH*) :test 'equal )
	   (pushnew `(,p2 1) (gethash `(,p0 ,tml) *VERTEX-GRAPH*) :test 'equal )
	   (pushnew `(,p3 1) (gethash `(,p0 ,tml) *VERTEX-GRAPH*) :test 'equal )
	   (pushnew `(,p4 -1) (gethash `(,p0 ,tml) *VERTEX-GRAPH*) :test 'equal )
	   (pushnew `(,p0 1) (gethash `(,p1 ,tmh) *VERTEX-GRAPH*) :test 'equal )
	   (pushnew `(,p2 1) (gethash `(,p1 ,tmh) *VERTEX-GRAPH*) :test 'equal )
	   (pushnew `(,p3 1) (gethash `(,p1 ,tmh) *VERTEX-GRAPH*) :test 'equal )
	   (pushnew `(,p4 -1) (gethash `(,p1 ,tmh) *VERTEX-GRAPH*) :test 'equal )
	   (pushnew `(,p0 1) (gethash `(,p2 ,tmh) *VERTEX-GRAPH*) :test 'equal )
	   (pushnew `(,p1 1) (gethash `(,p2 ,tmh) *VERTEX-GRAPH*) :test 'equal )
	   (pushnew `(,p3 1) (gethash `(,p2 ,tmh) *VERTEX-GRAPH*) :test 'equal )
	   (pushnew `(,p4 -1) (gethash `(,p2 ,tmh) *VERTEX-GRAPH*) :test 'equal )
	   (pushnew `(,p0 1) (gethash `(,p3 ,tmh) *VERTEX-GRAPH*) :test 'equal )
	   (pushnew `(,p1 1) (gethash `(,p3 ,tmh) *VERTEX-GRAPH*) :test 'equal )
	   (pushnew `(,p2 1) (gethash `(,p3 ,tmh) *VERTEX-GRAPH*) :test 'equal )
	   (pushnew `(,p4 -1) (gethash `(,p3 ,tmh) *VERTEX-GRAPH*) :test 'equal )
	   (pushnew `(,p0 -1) (gethash `(,p4 ,tmh) *VERTEX-GRAPH*) :test 'equal )
	   (pushnew `(,p1 -1) (gethash `(,p4 ,tmh) *VERTEX-GRAPH*) :test 'equal )
	   (pushnew `(,p2 -1) (gethash `(,p4 ,tmh) *VERTEX-GRAPH*) :test 'equal )
	   (pushnew `(,p3 -1) (gethash `(,p4 ,tmh) *VERTEX-GRAPH*) :test 'equal )))))


;--------|;--------|;--------|;--------|;--------|;--------|;--------|;--------|
(defun create-vertex-graph ()
  "after a spacetime has been loaded, calling this function will generate the
vertex graph"
  (maphash #'(lambda (k v)
	       (declare (ignore k))
	       (4simplex->vertices v))
	   *ID->4SIMPLEX*))
	       
|#
