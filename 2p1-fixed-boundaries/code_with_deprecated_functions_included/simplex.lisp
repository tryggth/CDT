;;;;---------------------------------------------------------------------;;;;
;;;;                   SIMPLEX DATA STRUCTURES                           ;;;;
;;;;---------------------------------------------------------------------;;;;

;;;; This file contains simplex information for the main data
;;;; structures of the simulation. It also contains important
;;;; parameters and functions for fixed boundary conditions.

;; Authors:
;; ------- Rajesh Kommu
;; ------- David Kamensky
;; ------- Jonah Miller (jonah.maxwell.miller@gmail.com)

;;;;---------------------------------------------------------------------

;; Sets the optimization parameters for compilation.
;; JM: I'm not sure why this is here. I'll leave it for now.
(declaim (optimize (speed 3)
		   (compilation-speed 0)
		   (debug 0)
		   (safety 0)))


;; Macro to decide whether or not to apply a modulo division,
;; depending on whether or not periodic boundary conditions are in
;; use.

;; JM: I'm not convinced I should enable periodic boundary conditions
;; at all for this program, since Rajesh's original code does periodic
;; boundary conditions just fine. It might be better to make a wrapper
;; for all of the simulation programs we have. I think this would make
;; the simulation a little faster. I'll have to think about that.
(defmacro bc-mod (num)
  `(if (or (string= BCTYPE "PERIODIC") *merge-faces*)
       (mod ,num NUM-T)
       ,num))

;;; sl 1-simplex is (tslice (p0 p1))
;;; tl 1-simplex is (type tmlo (p0 p1)) where type = 1
;;; sl 2-simplex is (tslice (p0 p1 p2))
;;; tl 2-simplex is (type tmlo (p0 p1 p2)) where type = 1,2
;;; tl 3-simplex is (type tmlo tmhi (p0 p1 p2 p3) (n0 n1 n2 n3))
;;; where type = 1,2,3,4 
;;; nj = id of the 4sx that does not have pj

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
;; Useful during the moves, when the points of the simplex are computed 
;; via append / unions / intersections in the form of a list.
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
    ;; set last used point if a new point is used
    ;; Required by David's initialization.lisp code.
    (setf *LAST-USED-POINT* (max *LAST-USED-POINT* p0 p1 p2 p3))
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


;; Some useful macros
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
	       (cond 
		 ((= 1 (3sx-type 3sx))
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

;; returns (1ids 2ids 3ids) where 1ids is a list of (1,3) ids in the sandwich
;;  etc.
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
		 (when (and (= (3sx-tmlo sx) (bc-mod tlo)) 
			    (= (3sx-tmhi sx) (bc-mod thi)))
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

;;; JM : The next few count functions were written by David Kamensky
;;; for deprecated code and I don't know if they work. They are
;;; required for initializing the simulation and setting the b-vector
;;; initially, so it is important that they work correctly. Rajesh's
;;; updated code with more data structures should make these functions
;;; much easier. Assuming the hash table format is the same, they
;;; should work fine, but use at your own risk.


;;; TODO: It seems that many of these count functions only work in the
;;; fixed boundary case. Fix this using David Kamensky's bc-mod macro.


;; JM : I think this one will still work, because it only cares about
;; 3-simplexes.
(defun count-boundary-vs-bulk ()
  "Counts the number of 3-simplexes with edges on the boundary
compared to those only in the bulk."
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


;; Function to count objects in a single time slice over the entire
;; spacetime. I know this could be written as a macro, but it works
;; passably as a function, so we'll do that. It's true that if the
;; boundary conditions are periodic, this count will look at the empty
;; NUM-T time slice, but since the slice is empty (and there's no
;; periodic boundary conditions modulo going on here), the result of
;; the count there will be zero.
(defun count-over-all-spacetime-slices (counting-function 
					&optional (count-argument nil))
  "Runs a counting function of a single time slice over all spacetime
slices and returns the sum. If the spacetime. If a count-argument is
given, this is passed to the counting-function as the second argument
after time-slice."
  (let ((count 0)
	(tmax (if (string= BCTYPE "PERIODIC")
		  (1- NUM-T)
		  NUM-T)))
    (if (not count-argument)
	(loop for time-slice from 0 to tmax do
	     (setf count (+ count (funcall counting-function time-slice))))
	(loop for time-slice from 0 to tmax do
	     (setf count (+ count (funcall 
				   counting-function 
				   time-slice 
				   count-argument)))))
    count))
       

;; Function to count objects in a sandwich over the entire
;; spacetime. I know this could be written as a macro, but it works
;; passably as a function, so we'll do that. It's true that if the
;; boundary conditions are periodic, this count will look at the empty
;; NUM-T-1, NUM-T, sandwich, but since the slice is empty (and there's no
;; periodic boundary conditions modulo going on here), the result of
;; the count there will be zero.
(defun count-over-all-spacetime-sandwiches (counting-function 
					    &optional (count-argument nil))
  "Runs a counting function of a single sandwich over all sandwiches
in the spacetime and returns the sum. If a count-argument is given,
this is passed to the counting-function as the second argument after
time-slice."
  (let ((count 0))
    (if (not count-argument)
	(loop for time-slice from 0 to (1- NUM-T) do
	     (setf count (+ count (funcall counting-function 
					   time-slice
					   (1+ time-slice)))))
	(loop for time-slice from 0 to (1- NUM-T) do
	     (setf count (+ count (funcall 
				   counting-function 
				   time-slice 
				   (1+ time-slice)
				   count-argument)))))
    count))

    
;; function to count the number of points in a particular time slice
;; JM: This function could probably benefit from Rajesh's new data
;; structures. However, I believe the benefit is minimal, since all
;; the necessary information is still in the top-level hash-table.
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


;;; JM: I've added this function to count all the points in the entire
;;; spacetime for use in the initialization phase.
(defun count-points-in-spacetime ()
  "Count points in the entire spacetime."
  (count-over-all-spacetime-slices #'count-points-at-time))


(defun get-timelike-links-in-sandwich (t-low t-high)
  "Retrieve the timelike links in a sandwich. Useful for initialization 
  and debugging."
  t-high
  (list-keys-with-trait #'(lambda (x) (= x (bc-mod t-low)))
			*TL1SIMPLEX->ID* 1))


;;; JM: I have changed count-timelike-links-in-sandwich to take
;;; advantage of Rajesh's bug-fixed code, including the new
;;; sub-simplex hash tables. Thus, the old
;;; count-timelike-links-in-sandwich is deprecated. It has been marked
;;; as such, and I leave it only for completeness sake. Use at your
;;; own risk.

;; Function to count the number of timelike links in a sandwich Since
;; we have a list of timelike links waiting for us in a hash table, we
;; simply count up the ones with the low element tmlo=t0. Order
;; matters. Since order matters for get-simplices-in-sandwich, this is
;; fine. Note: This function calls a parameter it never uses for
;; compatibility reasons. All we care about is t-low.
(defun count-timelike-links-in-sandwich (t-low t-high)
  "Count time-like links in a sandwich. 
t-high is for compatibility. We only care about t-low."
  t-high
  (count-keys-with-trait #'(lambda (x) (= x (bc-mod t-low))) 
			 *TL1SIMPLEX->ID* 1))


;; Function to count the number of timelike links in a sandwich
;; DEPRECATED: USE AT YOUR OWN RISK
(defun count-timelike-links-in-sandwich-deprecated (t0 t1)
 
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
				   :test #'(lambda (x y) 
					     (not (set-difference x y)))))))
    (length list-of-links)))

;; Count the number of time-like links in the entire spacetime. Needs
;; to be careful of indices to avoid double-counting
(defun count-timelike-links-in-spacetime ()
  "Count all the timelike links in the entire spacetime."
  (count-over-all-spacetime-sandwiches #'count-timelike-links-in-sandwich))

(defun get-spacelike-triangles-at-time (t0)
  "Get the keys for spacelike triangles at a given time (since the 
   hash table matches the simplex to the IDs, we want the simplex, 
   not the id.)"
  (list-keys-with-trait #'(lambda (x) (= x (bc-mod t0))) *SL2SIMPLEX->ID* 0))


;;; JM: I have changed count-spacelike-triangles-at-time to take
;;; advantage of Rajesh's bug-fixed code, including the new
;;; sub-simplex hash tables. Thus, the old
;;; function is deprecated. It has been marked
;;; as such, and I leave it only for completeness sake. Use at your
;;; own risk.
(defun count-spacelike-triangles-at-time (t0)
  "Count spacelike triangles (2-simplices) on a given time-slice."
  (count-keys-with-trait #'(lambda (x) (= x (bc-mod t0))) *SL2SIMPLEX->ID* 0))

;;count the number of spacelike triangles at a certain time
;; Although this is deprecated, it should work anyways.
(defun count-spacelike-triangles-at-time-deprecated (t0)

  ;;must treat one end as a special case; choosing to treat
  ;;time 0 specially. number of triangles is just related to the
  ;;number of type 1/3 simplices in nieghboring sandwiches

  (if (= t0 0)
      (count-simplices-in-sandwich-of-type 0 1 3)
      (count-simplices-in-sandwich-of-type (1- t0) t0 1)))

;; Count the number of spacelike triangles at all times
(defun count-spacelike-triangles-in-spacetime ()
  "Count all the spacelike 2-simplices in the spacetime."
  (count-over-all-spacetime-slices #'count-spacelike-triangles-at-time))


(defun get-spacelike-links-at-time (t0)
  "Returns a list of spacelike links at a given time.
   This is the key to the hash table, not the value."
  (list-keys-with-trait #'(lambda (x) (= x (bc-mod t0))) *SL1SIMPLEX->ID* 0))


;;; JM: I have changed count-spacelike-links-at-time to take advantage
;;; of Rajesh's bug-fixed code, including the new sub-simplex hash
;;; tables. Thus, the old function is deprecated. It has been removed
;;; becuase it relied on an outmoded version of euler-char, which has
;;; been since removed.

(defun count-spacelike-links-at-time (t0)
  "Count number of spacelike links on a given time slice."
  (count-keys-with-trait #'(lambda (x) (= x (bc-mod t0))) *SL1SIMPLEX->ID* 0))
   

;; Count the total number of spacelike links in the spacetime.
(defun count-spacelike-links-in-spacetime ()
  "Count the total number of spacelike links in all of space and time."
  (count-over-all-spacetime-slices #'count-spacelike-links-at-time))


(defun get-timelike-triangles-in-sandwich (t-low t-high)
  "Retrieve a list of all timelike triangles in a sandwich."
  t-high
  (list-keys-with-trait #'(lambda (x) (= x t-low)) *TL2SIMPLEX->ID* 1))

;;; JM: I have modified count-timelike-triangles-in-sandwich to take
;;; advantage of Rajesh's new and improved data structures. The first
;;; function is mine. The second one is deprecated, left in for
;;; completeness only. Use at your own risk.

;;; TODO: It seems that these count functions only work in the fixed
;;; boundary case. Fix this using David Kamensky's bc-mod macro.

(defun count-timelike-triangles-in-sandwich (t0 t1)
  "Count the number of timelike triangles in a sandwich."
  t1
  (count-keys-with-trait #'(lambda (x) (= t0 x)) *TL2SIMPLEX->ID* 1))

;; count the number of timelike triangles in a particular sandwich
;; DEPRECATED USE AT YOUR OWN RISK
(defun count-timelike-triangles-in-sandwich-deprecated (t0 t1)
  
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
				       :test #'(lambda (x y) 
						 (not 
						  (set-difference x y)))))))
    (length list-of-triangles)))


(defun count-timelike-triangles-in-spacetime ()
  "Count the total number of spacelike triangles in the entire spacetime."
  (count-over-all-spacetime-sandwiches 
   #'count-timelike-triangles-in-sandwich))

;; JM: The following two functions are useful for open boundary
;; conditions. They're mostly convenient for debugging.

;; Get the 13-simplices with faces on the final time slice.
(defun topfaces nil
  (let ((tf    NUM-T)
	(tf-1 (- NUM-T 1))
	(13-simplices 1))
    (get-simplices-in-sandwich-of-type tf-1 tf 13-simplices)))

;; Get the 31-simplices with faces on the initial time slice.
(defun bottomfaces nil
  (let ((31-simplices 3)
	(ti   0)
	(ti+1 1))
    (get-simplices-in-sandwich-of-type ti ti+1 31-simplices)))

;; JM: The next function produces a list of all 3-simplex IDs
;; currently in use. This function is SLOOOW, but useful for
;; tuning and debugging.
(defun list-all-3-simplices ()
  "Lists all 3-simplex ids in the spacetime."
  (let ((id-list nil)
	(val-list nil))
    (maphash #'(lambda (key val) (push key id-list) (push val val-list))
	     *ID->3SIMPLEX*)
    id-list))

(defun count-all-3-simplices ()
  "Counts all the 3-simplices in the spacetime. Very slow."
  (length (list-all-3-simplices)))

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

;; JM: Style rules broken here because it's clearer for the string.
(defun save-spacetime-to-file (outfile)
  (format outfile "~A ~A ~A ~A ~A ~A ~A ~A ~A ~A ~A ~A ~A ~A ~A ~A ~A ~A ~A ~A ~A ~A ~A~%" 
	  BCTYPE STOPOLOGY NUM-T N-INIT *LAST-USED-POINT* *LAST-USED-3SXID* 
	  N0 N1-SL N1-TL N2-SL N2-TL N3-TL-31 N3-TL-22 
	  *N1-SL-TOP* *N3-22-TOP* *N3-31-TOP*
	  *N1-SL-BOT* *N3-22-BOT* *N3-31-BOT*
	  *eps* *k0* *k3* *alpha*)
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
      (setf BCTYPE            (nth 0  data) 
	    STOPOLOGY         (nth 1  data) 
	    NUM-T             (nth 2  data) 
	    N-INIT            (nth 3  data) 
	    *LAST-USED-POINT* (nth 4  data) 
	    *LAST-USED-3SXID* (nth 5  data)
	    N0                (nth 6  data)
	    N1-SL             (nth 7  data) 
	    N1-TL             (nth 8  data)
	    N2-SL             (nth 9  data)
	    N2-TL             (nth 10 data)
	    N3-TL-31          (nth 11 data)
	    N3-TL-22          (nth 12 data)
	    *N1-SL-TOP*       (nth 13 data)
	    *N3-22-TOP*       (nth 14 data)
	    *N3-31-TOP*       (nth 15 data)
	    *N1-SL-BOT*       (nth 16 data)
	    *N3-22-BOT*       (nth 17 data)
	    *N3-31-BOT*       (nth 18 data)
	    *eps* (nth 19 data))
      (set-k0-k3-alpha (nth 20 data) (nth 21 data) (nth 22 data)))))

(defun parse-simplex-data-line (line)
  (with-input-from-string (s line)
    (let ((data (loop
		   :for num := (read s nil nil)
		   :while num
		   :collect num)))
      (make-3simplex-v4 (nth 0 data) (nth 1 data)  (nth 2 data) 
			(nth 3 data) (nth 4 data)  (nth 5 data) 
			(nth 6 data) (nth 7 data)  (nth 8 data)
			(nth 9 data) (nth 10 data) (nth 11 data)))))

(defun load-spacetime-from-file (infile)
  (parse-parameters-line (read-line infile nil))
  (loop for line = (read-line infile nil)
     while line do (parse-simplex-data-line line)))

;;; spatial 2-simplex is a spatial triangle; this information is
;;; needed for computing the spectral and hausdorff dimensions of the
;;; spatial slice, among other things.

;;; JM: there may be some redundancy with counting functions down here
;;; and those above. I will leave it alone so I don't break anything.
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
	  (let ((pos1 (position (first (set-difference points1 line))
				points1))
		(pos2 (position (first (set-difference points2 line))
				points2)))
	    (setf (nth pos1 (s2sx-sx2ids st1)) st2id 
		  (nth pos2 (s2sx-sx2ids st2)) st1id)))))))

(defun connect-s2simplices-within-list (sx1ids)
  (for (n 0 (1- (length sx1ids)))
       (for (m (1+ n) (1- (length sx1ids)))
	    (connect-s2simplices (nth n sx1ids) (nth m sx1ids)))))

(defun count-s2simplices-in-slice (ts)
  "counts the number of spaeclike 2-simplices
   (i.e. spatial triangles) in spatial slice ts"
  (let ((count 0))
    (maphash #'(lambda (id sx)
		 (declare (ignore id))
		 (when (= (s2sx-time sx) (bc-mod ts)) 
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
      (make-s2simplex-v2 (nth 0 data) (nth 1 data) (nth 2 data)
			 (nth 3 data) (nth 4 data) (nth 5 data)
			 (nth 6 data) (nth 7 data)))))

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
  "3sx2p1->s2sx2p1 generates the spatial 2-simplex information 
  for each spatial slice from the 3-simplex data for the entire spacetime."
  (clrhash *ID->SPATIAL-2SIMPLEX*)
  (setf *LAST-USED-S2SXID* 0)
  (for (ts 0 (1- NUM-T))
       ;; list of ids
       (let ((31simplices (get-simplices-in-sandwich-of-type ts (1+ ts) 3)) 
	     (spatial-triangles '()))
	 (dolist (31simplex 31simplices)
	   (push (make-s2simplex ts (3sx-lopts (get-3simplex 31simplex))) 
		 spatial-triangles))
	 (connect-s2simplices-within-list spatial-triangles))))

(defun generate-s2sx2p1-files (3sx2p1files)
  "3sx2p1files is a list of .3sx2p1 files. For each file in this list, this
   function generates a .s2sx2p1 file. The prefix for the .3sx2p1 and the 
   .s2sx2p1 file are identical"
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

	     
