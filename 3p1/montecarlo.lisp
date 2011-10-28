;; cdt-3plus1-montecarlo.lisp
(defun try-move (sxid mtype)
  (cond ((= 0 mtype) (try-2->8 sxid))
	((= 1 mtype) (try-4->6 sxid))
	((= 2 mtype) (try-2->4-v1 sxid))
	((= 3 mtype) (try-2->4-v2 sxid))
	((= 4 mtype) (try-3->3-v1 sxid))
	((= 5 mtype) (try-3->3-v2 sxid))
	((= 6 mtype) (try-4->2-v2 sxid))
	((= 7 mtype) (try-4->2-v1 sxid))
	((= 8 mtype) (try-6->4 sxid))
	((= 9 mtype) (try-8->2 sxid))))

(defun move->string (mtype)
  (cond ((= 0 mtype) "2->6")
	((= 1 mtype) "4->6")
	((= 2 mtype) "2->4-v1")
	((= 3 mtype) "2->4-v2")
	((= 4 mtype) "3->3-v1")
	((= 5 mtype) "3->3-v2")
	((= 6 mtype) "4->2-v2")
	((= 7 mtype) "4->2-v1")
	((= 8 mtype) "6->4")
	((= 9 mtype) "8->2")))

(defun random-move (nsweeps)
  (loop :for sweepnum :from 1 :to nsweeps
     do
     (let* ((id (random *LAST-USED-4SXID*))
	    (mtype (select-move))
	    (sx (get-4simplex id))
	    (movedata nil))

       (incf CURRENT-MOVE-NUMBER)
       (when (and sx (setf movedata (try-move id mtype)))
	 (3plus1move movedata))
       (when (= 0 (mod sweepnum 10000))
	 (format t "finished ~A of ~A sweeps with count ~A~%" sweepnum nsweeps 
		 (count-simplices-of-all-types))
	 (finish-output)))))

(defun accept-move? (mtype)
  (let ((delta-action 0.0)
	(delta-damping 0.0))
    (cond ((= 0 mtype)
	   (setf delta-damping (- (damping (+ N4-TL-41 6) (+ N4-TL-32 0)) (damping N4-TL-41 N4-TL-32)))
	   (setf delta-action (- (action (+ N0 1) (+ N4-TL-41 6) (+ N4-TL-32 0))
				 (action N0 N4-TL-41 N4-TL-32))))
	  ((= 1 mtype)
	   (setf delta-damping (- (damping (+ N4-TL-41 2) (+ N4-TL-32 0)) (damping N4-TL-41 N4-TL-32)))
	   (setf delta-action (- (action (+ N0 0) (+ N4-TL-41 2) (+ N4-TL-32 0))
				 (action N0 N4-TL-41 N4-TL-32))))
	  ((= 2 mtype)
	   (setf delta-damping (- (damping (+ N4-TL-41 0) (+ N4-TL-32 2)) (damping N4-TL-41 N4-TL-32)))
	   (setf delta-action (- (action (+ N0 0) (+ N4-TL-41 0) (+ N4-TL-32 2))
				 (action N0 N4-TL-41 N4-TL-32))))
	  ((= 3 mtype)
	   (setf delta-damping (- (damping (+ N4-TL-41 0) (+ N4-TL-32 2)) (damping N4-TL-41 N4-TL-32)))
	   (setf delta-action (- (action (+ N0 0) (+ N4-TL-41 0) (+ N4-TL-32 2))
				 (action N0 N4-TL-41 N4-TL-32))))
	  ((= 4 mtype)
	   (setf delta-damping (- (damping (+ N4-TL-41 0) (+ N4-TL-32 0)) (damping N4-TL-41 N4-TL-32)))
	   (setf delta-action (- (action (+ N0 0) (+ N4-TL-41 0) (+ N4-TL-32 0))
				 (action N0 N4-TL-41 N4-TL-32))))
	  ((= 5 mtype)
	   (setf delta-damping (- (damping (+ N4-TL-41 0) (+ N4-TL-32 0)) (damping N4-TL-41 N4-TL-32)))
	   (setf delta-action (- (action (+ N0 0) (+ N4-TL-41 0) (+ N4-TL-32 0))
				 (action N0 N4-TL-41 N4-TL-32))))
	  ((= 6 mtype)
	   (setf delta-damping (- (damping (+ N4-TL-41 0) (+ N4-TL-32 -2)) (damping N4-TL-41 N4-TL-32)))
	   (setf delta-action (- (action (+ N0 0) (+ N4-TL-41 0) (+ N4-TL-32 -2))
				 (action N0 N4-TL-41 N4-TL-32))))
	  ((= 7 mtype)
	   (setf delta-damping (- (damping (+ N4-TL-41 0) (+ N4-TL-32 -2)) (damping N4-TL-41 N4-TL-32)))
	   (setf delta-action (- (action (+ N0 0) (+ N4-TL-41 0) (+ N4-TL-32 -2))
				 (action N0 N4-TL-41 N4-TL-32))))
	  ((= 8 mtype)
	   (setf delta-damping (- (damping (+ N4-TL-41 -2) (+ N4-TL-32 0)) (damping N4-TL-41 N4-TL-32)))
	   (setf delta-action (- (action (+ N0 0) (+ N4-TL-41 -2) (+ N4-TL-32 0))
				 (action N0 N4-TL-41 N4-TL-32))))
	  ((= 9 mtype)
	   (setf delta-damping (- (damping (+ N4-TL-41 -6) (+ N4-TL-32 0)) (damping N4-TL-41 N4-TL-32)))
	   (setf delta-action (- (action (+ N0 -1) (+ N4-TL-41 -6) (+ N4-TL-32 0))
				 (action N0 N4-TL-41 N4-TL-32)))))
    (< (random 1.0) (* (exp (* -1.0 delta-action))
		       (exp (* -1.0 delta-damping))))))



(defun sweep ()
  (let ((num-attempted 0))
    (while (< num-attempted N-INIT)
      (let* ((sxid (random *LAST-USED-4SXID*))
	     (mtype (select-move))
	     (movedata (try-move sxid mtype)))
	(while (null movedata) 
	  (setf sxid (random *LAST-USED-4SXID*) 
		mtype (select-move) 
		movedata (try-move sxid mtype)))
	(incf num-attempted) ;; number-of-attempted-moves-counter for this sweep
	(incf (nth mtype ATTEMPTED-MOVES)) ;; number of moves of mtype that have been attempted
	(when (accept-move? mtype)
	  (incf (nth mtype SUCCESSFUL-MOVES)) ;; number of moves of mtype that have succeeded
	  (3plus1move movedata))))))

;; the following is used for tuning KAPPA-0, DELTA, KAPPA-4 parameters
(defun generate-data-console (&optional (start-sweep 1))
  (for (ns start-sweep (+ start-sweep NUM-SWEEPS -1))
       (sweep)
       (when (= 0 (mod ns 10))
	 (format t "~A/~A/~A ~A ~A~%"
		 start-sweep ns (+ start-sweep NUM-SWEEPS -1) 
		 (count-simplices-of-all-types) (accept-ratios)))
       (finish-output)))

(defun generate-data (&optional (start-sweep 1))
  (let ((datafilestr (concatenate 'string (generate-filename start-sweep) 4SXEXT))
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
	     (format progfile "~A/~A/~A ~A~%"
		     start-sweep ns end-sweep (count-simplices-of-all-types)))))))

;; generate-data-v2 is similar to generate data except it creates a fresh data file every 
;; SAVE-EVERY-N-SWEEPS. since a fresh datafile is created, there is no need to maintain a seprate progress
;; file.
(defun generate-data-v2 (&optional (start-sweep 1))
  (setf SIM-START-TIME (cdt-now-str))
  (let ((end-sweep (+ start-sweep NUM-SWEEPS -1)))
    (when (= 1 start-sweep) ;; save the initial spacetime if this is a brand new run
      (with-open-file (datafile (concatenate 'string (generate-filename-v2 start-sweep 0) 4SXEXT)
				:direction :output
				:if-exists :supersede)
	(save-spacetime-to-file datafile)))
    (for (ns start-sweep end-sweep)
	 (sweep)
	 (when (= 0 (mod ns SAVE-EVERY-N-SWEEPS))
	   (with-open-file (datafile (concatenate 'string (generate-filename-v2 start-sweep ns) 4SXEXT)
				     :direction :output
				     :if-exists :supersede)
	     (save-spacetime-to-file datafile))))))

;; generate-data-v3 is similar to generate-data-v2 except it also creates an additional data file every 
;; SAVE-EVERY-N-SWEEPS that contains the spatial 3-simplex information for each spatial slice.
(defun generate-data-v3 (&optional (start-sweep 1))
  (setf SIM-START-TIME (cdt-now-str))
  (let ((end-sweep (+ start-sweep NUM-SWEEPS -1)))
    (when (= 1 start-sweep) ;; save the initial spacetime contents if this is a brand new run
      (with-open-file (datafile (concatenate 'string (generate-filename-v2 start-sweep 0) 4SXEXT)
				:direction :output
				:if-exists :supersede)
	(save-spacetime-to-file datafile)))
    (for (ns start-sweep end-sweep)
	 (sweep)
	 (when (= 0 (mod ns SAVE-EVERY-N-SWEEPS))
	   (let ((filename (generate-filename-v2 start-sweep ns)))
	     (with-open-file (datafile (concatenate 'string filename 4SXEXT)
				       :direction :output
				       :if-exists :supersede)
	       (save-spacetime-to-file datafile))
	     (4sx3p1->s3sx3p1)
	     (with-open-file (datafile (concatenate 'string filename S3SXEXT)
				       :direction :output
				       :if-exists :supersede)
	       (save-s3simplex-data-to-file datafile)))))))

(defun generate-movie-data (&optional (start-sweep 1))
  (setf SAVE-EVERY-N-SWEEPS 10)
  (let ((moviefilestr (concatenate 'string (generate-filename start-sweep) MOVEXT))
	(trackfilestr (concatenate 'string (generate-filename start-sweep) PRGEXT))
	(end-sweep (+ start-sweep NUM-SWEEPS -1)))

    ;; open and close the file for :append to work
    (with-open-file (moviefile moviefilestr 
			       :direction :output
			       :if-exists :supersede)
      ;; record the initial data only if start-sweep = 1
      (when (= start-sweep 1)

	(for (ts 0 (1- NUM-T))
	     (format moviefile "~A " (count-simplices-in-sandwich ts (1+ ts))))
	(format moviefile "~%")))
    
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
	     (format trackfile "~A/~A/~A ~A~%"
		     start-sweep ns end-sweep (count-simplices-of-all-types)))))))

(defun generate-movie-data-console (&optional (start-sweep 1))
  (when (= 1 start-sweep)
    (reset-move-counts))
  (let ((end-sweep (+ start-sweep NUM-SWEEPS -1)))
       (for (ns start-sweep end-sweep)
	   (sweep)
	   (when (= 0 (mod ns 10))
	     (format t "~A/~A " ns end-sweep)
	     (for (ts 0 (1- NUM-T))
		  (format t "~A " (count-simplices-in-sandwich ts (1+ ts))))
	     (format t "~A ~A~%" (count-simplices-of-all-types) (accept-ratios))))))

(defun generate-spacetime-and-movie-data (&optional (start-sweep 1))
  (let* ((end-sweep (+ start-sweep NUM-SWEEPS -1))
	 (datafilestr (concatenate 'string (generate-filename start-sweep) 4SXEXT))
	 (trackfilestr (concatenate 'string (generate-filename start-sweep) PRGEXT))
	 (moviefilestr (concatenate 'string (generate-filename start-sweep) MOVEXT)))
    
    ;; open and close the file, for :append to work properly
    (with-open-file (moviefile moviefilestr 
			       :direction :output
			       :if-exists :supersede)
      ;; record the initial data only if start-sweep = 1
      (when (= 1 start-sweep)
	(for (ts 0 (1- NUM-T))
	     (format moviefile "~A " (count-simplices-in-sandwich ts (1+ ts))))
	(format moviefile "~%")))
    
    (for (ns start-sweep end-sweep)
	 (sweep)
	 (when (= 0 (mod ns SAVE-EVERY-N-SWEEPS))
	   (with-open-file (datafile datafilestr 
				     :direction :output
				     :if-exists :supersede)
	     (save-spacetime-to-file datafile))
	   
	   (with-open-file (trackfile trackfilestr
				      :direction :output
				      :if-exists :supersede)
	     (format trackfile "~A/~A/~A ~A~%" start-sweep ns end-sweep (count-simplices-of-all-types)))
	   
	   (with-open-file (moviefile moviefilestr 
				      :direction :output
				      :if-exists :append)
	     (for (ts 0 (1- NUM-T))
		  (format moviefile "~A " (count-simplices-in-sandwich ts (1+ ts))))
	     (format moviefile "~%"))))))

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