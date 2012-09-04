;;;; spacelike_2simplices.lisp
;;;; Authors: Rajesh Kommu,
;;;;          Jonah Miller (jonah.maxwell.miller@gmail.com)

;;; spatial 2-simplex is a spatial triangle; this information is
;;; needed for computing the spectral and hausdorff dimensions of the
;;; spatial slice, among other things.


;;; ID tracking and hash tables
(defparameter *ID->SPATIAL-2SIMPLEX* (make-hash-table))
(defparameter *LAST-USED-S2SXID* 0)
(defmacro next-s2simplex-id ()
  `(incf *LAST-USED-S2SXID*))


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


;;;; OUTPUT
(defun make-spatial-2-simplex-file (filename)
  "Like make-spacetime-file aboce, but generates spatial 2-simplex
   data instead of 3-simplex data. The filename input is made
   by (generate-filename start-sweep) or somthing similar. It does not
   have the file extension."
  (3sx2p1->s2sx2p1)
  (with-open-file (datafile (concatenate 'string filename S2SXEXT)
			    :direction :output
			    :if-exists :supersede)
    (save-s2simplex-data-to-file datafile)))

;; generate-data-v3 is similar to generate-data-v2 except it also
;; creates an additional data file every SAVE-EVERY-N-SWEEPS that
;; contains the spatial 2-simplex information for each spatial slice.
(defun generate-data-v3 (&optional (start-sweep 1))
  "Like generate-data-v2, but generates spatial 2 simplex and 
   3-simplex information every SAVE-EVERY-N-SWEEPS."
  (setf SIM-START-TIME (cdt-now-str))
  (let ((end-sweep (+ start-sweep NUM-SWEEPS -1)))
    (when (= 1 start-sweep) ;; save the initial spacetime contents if
			    ;; this is a brand new run
      (make-spacetime-file (generate-filename-v2 start-sweep 0))
      (make-spatial-2-simplex-file (generate-filename-v2 start-sweep 0)))
    (for (ns start-sweep end-sweep)
      (sweep)
      (make-spacetime-file (generate-filename-v2 start-sweep ns))
      (make-spatial-2-simplex-file (generate-filename-v2 start-sweep ns)))))
