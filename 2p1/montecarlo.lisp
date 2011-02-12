; cdt-2plus1-montecarlo.lisp
(defun try-move (sxid mtype)
  (ecase mtype
    (0 (try-2->6 sxid))
    (1 (try-2->3 sxid))
    (2 (try-4->4 sxid))
    (3 (try-3->2 sxid))
    (4 (try-6->2 sxid))))

(defun random-move (nsweeps)
  (loop :for sweepnum :from 1 :to nsweeps
     do
     (let* ((id (mt19937:random *LAST-USED-3SXID*))
	    (mtype (mt19937:random 5))
	    (sx (get-3simplex id))
	    (movedata nil))

       (incf CURRENT-MOVE-NUMBER)
       (when (and sx (setf movedata (try-move id mtype)))
	 (2plus1move movedata))
       (when (= 0 (mod sweepnum 1000))
	 (format t "finished ~A of ~A sweeps with count ~A~%" sweepnum nsweeps 
		 (count-simplices-of-all-types))
	 (finish-output)))))

(defun action (num-0 num-3)
  (+ (- (* k3 num-3) (* k0 num-0)) (* eps (abs (- num-3 N-INIT)))))

(defun accept-move? (move-type)
  (let ((delta-action 0.0))
    (cond ((= move-type 26MTYPE) 
	   (setf delta-action (- (action (+ N0 1) (+ (N3) 4)) (action N0 (N3)))))
	  ((= move-type 62MTYPE) 
	   (setf delta-action (- (action (- N0 1) (- (N3) 4)) (action N0 (N3)))))
	  ((= move-type 44MTYPE)
	   (setf delta-action 0.0))
	  ((= move-type 23MTYPE) 
	   (setf delta-action (- (action N0 (+ (N3) 1)) (action N0 (N3)))))
	  ((= move-type 32MTYPE) 
	   (setf delta-action (- (action N0 (- (N3) 1)) (action N0 (N3))))))
    (< (mt19937:random 1.0) (exp (- delta-action)))))

;; a sweep is defined as N-INIT number of attempted moves
(defun sweep ()
  (let ((num-attempted 0))
    (while (< num-attempted N-INIT)
      (let* ((sxid (mt19937:random *LAST-USED-3SXID*))
	     (mtype (mt19937:random 5))
	     (sx (get-3simplex sxid))
	     (movedata nil))
	(when (and sx (setf movedata (try-move sxid mtype)))
	  (incf num-attempted)
	  (when (accept-move? mtype)
	    (2plus1move movedata)))))))

;; following is to be used for tuning the k0 and k3 parameters
(defun generate-data-console (&optional (start-sweep 1))
  (for (ns start-sweep (+ start-sweep NUM-SWEEPS -1))
       (sweep)
       (when (= 0 (mod ns 10))
	 (format t "start = ~A end = ~A current = ~A count = ~A~%"
		 start-sweep (+ start-sweep NUM-SWEEPS -1) ns (count-simplices-of-all-types)))
       (finish-output)))

;; generate-data should be called after setting the values for eps, k0, k3,
;; NUM-SWEEPS and calling one of the initialize-xx-slices.
(defun generate-data (&optional (start-sweep 1))
  (let ((datafilestr (concatenate 'string (generate-filename start-sweep) 3SXEXT))
	(progfilestr (concatenate 'string (generate-filename start-sweep) PRGEXT))
	(end-sweep (+ start-sweep NUM-SWEEPS -1)))
    (for (ns start-sweep end-sweep)
	 (sweep)
	 (when (= 0 (mod ns SAVE-EVERY-N-SWEEPS))
	   (with-open-file (datafile datafilestr 
				     :direction :output
				     :if-exists :supersede)
	     (save-spacetime-to-file datafile))
	   (with-open-file (progfile progfilestr
				     :direction :output
				     :if-exists :supersede)
	     (format progfile "start = ~A end = ~A current = ~A count = ~A~%"
		     start-sweep end-sweep ns (count-simplices-of-all-types)))))))

;; generate-movie-data saves number of simplices every SAVE-EVERY-N-SWEEPS
(defun generate-movie-data (&optional (start-sweep 1))
  (setf SAVE-EVERY-N-SWEEPS 10)
  (let ((moviefilestr (concatenate 'string (generate-filename start-sweep) MOVEXT))
	(trackfilestr (concatenate 'string (generate-filename start-sweep) ".movprg"))
	(end-sweep (+ start-sweep NUM-SWEEPS -1)))
    (with-open-file (moviefile moviefilestr 
			       :direction :output
			       :if-exists :supersede)
      ;; record the initial data only if start-sweep = 1
      (when (= start-sweep 1)
	(for (ts 0 (1- NUM-T))
	     (format moviefile "~A " (count-simplices-in-sandwich ts (1+ ts)))))
	(format moviefile "~%"))
    
      (for (ns start-sweep end-sweep)
	   (sweep)
	   (when (= 0 (mod ns SAVE-EVERY-N-SWEEPS))
	     (with-open-file (moviefile moviefilestr 
					:direction :output
					:if-exists :append)
	       (for (ts 0 (1- NUM-T))
		    (format moviefile "~A " (count-simplices-in-sandwich ts (1+ ts))))
	       (format moviefile "~%"))
	     (with-open-file (trackfile trackfilestr
					:direction :output
					:if-exists :supersede)
	       (format trackfile "start = ~A end = ~A current = ~A count = ~A~%"
		       start-sweep end-sweep ns (count-simplices-of-all-types)))))))

#|
(defun calculate-order-parameter (&optional (start-sweep 1))
  (let* ((end-sweep (+ start-sweep NUM-SWEEPS -1))
	 (order-parameter 0.0)
	 (datafilestr (format nil 
			      "~A-~A-op-T~A_V~A_eps~A_kz~A_kt~A_sweeps~Ato~A.op" 
			      *topology* *boundary-conditions*
			      NUM-T N-INIT eps k0 k3 start-sweep end-sweep))
	 (trackfilestr (format nil 
			       "~A-~A-op-T~A_V~A_eps~A_kz~A_kt~A_sweeps~Ato~A.progress" 
			       *topology* *boundary-conditions*
			       NUM-T N-INIT eps k0 k3 start-sweep end-sweep)))
    (do ((ns start-sweep (1+ ns)) 
	 (tot 0.0 (incf tot (/ N3-TL-22 (N3)))))
	((> ns end-sweep) (setf order-parameter (/ tot NUM-SWEEPS)))
      (sweep)
      (when (= 0 (mod ns SAVE-EVERY-N-SWEEPS))
	(with-open-file (trackfile trackfilestr
				   :direction :output
				   :if-exists :supersede)
	  (format trackfile "start = ~A end = ~A current = ~A~%"
		  start-sweep end-sweep ns))))
    (with-open-file (datafile datafilestr
			      :direction :output
			      :if-exists :supersede)
      (format datafile "T=~A V=~A eps=~A k0=~A k3=~A start=~A end=~A op=~A~%" 
	      NUM-T N-INIT eps k0 k3 start-sweep end-sweep order-parameter))))

(defun calculate-volume-volume-correlator (&optional (start-sweep 1))
  (let* ((end-sweep (+ start-sweep NUM-SWEEPS -1))
	 (vvparams (make-array (+ NUM-T 1) :initial-element 0.0))
	 (dfilestr (format nil 
			   "~A-~A-vv-T~A_V~A_eps~A_kz~A_kt~A_sweeps~Ato~A.vv" 
			   *topology* *boundary-conditions*
			   NUM-T N-INIT eps k0 k3 start-sweep end-sweep))
	 (tfilestr (format nil 
			   "~A-~A-vv-T~A_V~A_eps~A_kz~A_kt~A_sweeps~Ato~A.prog" 
			   *topology* *boundary-conditions*
			   NUM-T N-INIT eps k0 k3 start-sweep end-sweep)))
    (do ((ns start-sweep (1+ ns)))
	((> ns end-sweep)
	 (do ((j 0 (1+ j))) ((> j NUM-T))
	   (setf (aref vvparams j) 
		 (/ (aref vvparams j) (* NUM-T NUM-T (/ NUM-SWEEPS 100))))))
      (sweep)
      (when (= 0 (mod ns SAVE-EVERY-N-SWEEPS))
	(do ((col 0 (incf col))) ((> col NUM-T))
	  (do ((ts 1/2 (1+ ts))) ((> ts NUM-T))
	    (incf (svref vvparams col)
		  (* (count-simplices-at-time ts)
		     (count-simplices-at-time-pbc (+ ts
						     (- col (/ NUM-T 2))))))))
	(with-open-file (tfile tfilestr 
			       :direction :output 
			       :if-exists :supersede)
	  (format tfile "start = ~A end = ~A current = ~A~%"
		  start-sweep end-sweep ns))))
    (with-open-file (dfile dfilestr
			   :direction :output
			   :if-exists :supersede)
      (format dfile "T=~A V=~A eps=~A k0=~A k3=~A start=~A end=~A vvp=~A~%" 
	      NUM-T N-INIT eps k0 k3 start-sweep end-sweep vvparams))))

(defun compute-spatial-slice-hausdorff-dimension ()
  "compute the hausdorff dimension of all the spatial slices")

(defun compute-thin-sandwich-hausdorff-dimension ()
  "a thin sandwich consists of two adjacent spatial slices")

(defun compute-spacetime-hausdorff-dimension ()
  "the hausdorff dimension of the entire spacetime")

(defun compute-spatial-slice-spectral-dimension ()
  "compute the spectral dimension of all the spatial slices")

(defun compute-thin-sandwich-spectral-dimension ()
  "a thin sandwich consists of two adjacent spatial slices")

(defun compute-spacetime-spectral-dimension ()
  "the spectral dimension of the entire spacetime")

|#