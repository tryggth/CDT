; cdt-2plus1-montecarlo.lisp

;; translate a move type into an actual move
(defun try-move (sxid mtype)
  (let ((sx (get-3simplex sxid)))
    (if (and sx (is-real-simplex sx))
	(ecase mtype
	  (0 (try-2->6 sxid))
	  (1 (try-2->3 sxid))
	  (2 (try-4->4 sxid))
	  (3 (try-3->2 sxid))
	  (4 (try-6->2 sxid)))
	nil)))

(defun random-move (nsweeps)
  (loop :for sweepnum :from 1 :to nsweeps
     do
     (let* ((id (random *LAST-USED-3SXID*))
	    (mtype (select-move))
	    (sx (get-3simplex id))
	    (movedata nil))

       (incf CURRENT-MOVE-NUMBER)
       (when (and sx (setf movedata (try-move id mtype)))
	 (2plus1move movedata))
       (when (= 0 (mod sweepnum 1000))
	 (format t "finished ~A of ~A sweeps with count ~A~%" sweepnum nsweeps 
		 (count-simplices-of-all-types))
	 (finish-output)))))

(defun accept-move? (mtype sxid)

  (let* (;;determine what change vectors to use, based on the type of move
	
	 ;;this is the change in the f-vector due to the move (see the DF*
	 ;;parameters and/or macros in "globals.lisp")
	 (DF
	  (ecase mtype
	    (0;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; 2->6
	     (DF26 sxid))
	    (1;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; 2->3
	     DF23)
	    (2;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; 4->4
	     DF44)
	    (3;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; 3->2
	     DF32)
	    (4;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; 6->2
	     (DF62 sxid))))
	 
	 ;;this is the change in some relevant quantities at the boundary
	 ;;(see the macros/paramters DB* in "globals.lisp")
	 (DB
	  (ecase mtype
	    (0;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; 2->6
	     (DB26 sxid))
	    (1;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; 2->3
	     (DB23 sxid))
	    (2;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; 4->4
	     DB44)
	    (3;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; 3->2
	     (DB32 sxid))
	    (4;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; 6->2
	     (DB62 sxid))))
		 
	 ;;determine changes in different numbers of geometrical objects
	 (d-n3       (+ (seventh DF) (sixth DF)))
	 (d-n1-sl    (second  DF))
	 (d-n3-tl-22 (seventh DF))
	 (d-n3-tl-31 (sixth   DF))
	 (d-n1-tl    (third   DF))
	 (d-n1-sl-b  (first   DB))
	 (d-n3-22-b  (second  DB))
	 (d-n3-31-b  (third   DB))
	 
	 ;;determine deltas in damping and action
	 (delta-damping (- (damping (+ (N3) d-n3)) (damping (N3))))
	 (delta-action  (action  d-n1-sl   d-n1-tl   d-n3-tl-22  d-n3-tl-31 
				 d-n1-sl-b 
				 d-n3-22-b 
				 d-n3-31-b)))
    
    ;;accept with probability as function of changes in action and damping

    ;;this expression assumes that the "delta-action" is the change in
    ;;real-valued euclidean action
    (< (random 1.0) (* (exp (- delta-action))
                       (exp (- delta-damping))))))


;; a sweep is defined as N-INIT number of attempted moves
(defun sweep ()
  (let ((num-attempted 0))
    (while (< num-attempted N-INIT)
      (let* ((sxid (random *LAST-USED-3SXID*))
	     (mtype (select-move))
	     (movedata (try-move sxid mtype)))
	(while (null movedata)
	  (setf sxid (random *LAST-USED-3SXID*)
		mtype (select-move) 
		movedata (try-move sxid mtype)))
	(incf num-attempted) ;; number-of-attempted-moves-counter for this sweep
	(incf (nth mtype ATTEMPTED-MOVES)) ;; number of moves of mtype that have been attempted
	(when (accept-move? mtype sxid)
	  (incf (nth mtype SUCCESSFUL-MOVES)) ;; number of moves of mtype that have succeeded
	  (2plus1move movedata))))))

;; following is to be used for tuning the k0 and k3 parameters
(defun generate-data-console (&optional (start-sweep 1))
  (for (ns start-sweep (+ start-sweep NUM-SWEEPS -1))
       (sweep)
       (when (= 0 (mod ns 10))

;	 (format t "start = ~A end = ~A current = ~A count = ~A ~A\%~%"
;		 start-sweep (+ start-sweep NUM-SWEEPS -1) ns 
;		 (count-simplices-of-all-types) (percent-tv)));SUCCESSFUL-MOVES));(accept-ratios)))
;		 ;(count-simplices-of-all-types) (accept-ratios)))

	 (format t "start = ~A end = ~A current = ~A count = ~A ~$\%~%"
		 start-sweep (+ start-sweep NUM-SWEEPS -1) ns 
		 (count-boundary-vs-bulk) (percent-tv)))

       (finish-output)))

;; generate-data should be called after setting the values for eps, k0, k3,
;; NUM-SWEEPS and calling one of the initialize-xx-slices.
(defun generate-data (&optional (start-sweep 1))
  (setf SIM-START-TIME (cdt-now-str))
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

;; generate-data-v2 is similar to generate data except it creates a fresh data file every 
;; SAVE-EVERY-N-SWEEPS. since a fresh datafile is created, there is no need to maintain a seprate progress
;; file.
(defun generate-data-v2 (&optional (start-sweep 1))
  (setf SIM-START-TIME (cdt-now-str))
  (let ((end-sweep (+ start-sweep NUM-SWEEPS -1)))
    (when (= 1 start-sweep)
      (with-open-file (datafile (concatenate 'string (generate-filename-v2 start-sweep 0) 3SXEXT)
				:direction :output
				:if-exists :supersede)
	(save-spacetime-to-file datafile))
      (for (ns start-sweep end-sweep)
	   (sweep)
	   (when (= 0 (mod ns SAVE-EVERY-N-SWEEPS))
	     (with-open-file (datafile (concatenate 'string (generate-filename-v2 start-sweep ns) 3SXEXT)
				       :direction :output
				       :if-exists :supersede)
	       (save-spacetime-to-file datafile)))))))

;; generate-movie-data saves number of simplices every SAVE-EVERY-N-SWEEPS
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
	     (format trackfile "start = ~A end = ~A current = ~A count = ~A~%"
		     start-sweep end-sweep ns (count-simplices-of-all-types)))))))

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
	 (datafilestr (concatenate 'string (generate-filename start-sweep) 3SXEXT))
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

(defun 3sx2p1->2sx2p1 (infile outfile)
  "3sx2p1->2sx2p1 generates the 2-simplex information for each spatial slice from the 3-simplex data for
the entire spacetime. The generated information is written to outfile"
  (load-spacetime-from-file infile)
  (clrhash *ID->SPATIAL-2SIMPLEX*)
  (for (ts 0 (1- NUM-T))
       (let ((31simplices (get-simplices-in-sandwich-of-type ts (1+ ts) 3))
	     (spatial-triangles '()))
	 (dolist (31simplex 31simplices)
	   (push (make-s2simplex ts (3sx-lopts (get-3simplex 31simplex))) spatial-triangles))
	 (connect-spatial-2simplices-within-list spatial-triangles)))
  (save-s2simplex-data-to-file outfile))

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