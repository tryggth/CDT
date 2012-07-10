; cdt-2plus1-montecarlo.lisp

;; translate a move type into an actual move
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


;; JM: The new accept-move? function which calculates whether or not
;; to accept a move based on how the move affects the action. This is
;; seperate from the rejection of a move if try-move returns nil. This
;; function relies on the fact that the action action is linear in
;; N_simplexes such that a change in the number of simplexes we keep
;; track of translates directly into a change in the
;; action. accept-move? is now a function of the sxid because DB and
;; DF behave differently depending on simplex location.
(defun accept-move? (mtype sxid)
  (let* (;; Determine what vectors to use, based on the type of move

	 ;; This is the change in the f-vector due to the move (see
	 ;; the DF* parameters in "globals.lisp")
	 (DF 
	  (ecase mtype
	    (0 DF26)   ; 2->6 move
	    (1 DF23)   ; 2->3 move
	    (2 DF44)   ; 2->4 move
	    (3 DF32)   ; 3->2 move
	    (4 DF62))) ; 6->2 move
	 
	 ;; This is the change in some relevant quantities at the
	 ;; boundary (see DB* parameters in "globals.lisp")
	 (DB
	  (ecase mtype
	    (0 DB26)        ; 2->6 move
	    (1 (DB23 sxid)) ; 2->3 move
	    (2 DB44)        ; 4->4 move
	    (3 (DB32 sxid)) ; 3->2 move
	    (4 DB62)))      ; 6->2 move 
	 
	 ;; Determine the changes in different numbers of geometrical
	 ;; objects If your action doesn't require all of these
	 ;; quantities, comment some out. For clarity, please do not
	 ;; erase them.
;	 (d-n0         (nth 0 DF)) ; Not currently used
	 (d-n1-sl      (nth 1 DF)) ; Change in  N1-SL
	 (d-n1-tl      (nth 2 DF)) ; Chnange in N1-TL. etc.
;	 (d-n2-sl      (nth 3 DF)) ; Not currently used
;	 (d-n2-tl      (nth 4 DF)) ; Not currently used
	 (d-n3-tl-31   (nth 5 DF))
	 (d-n3-tl-22   (nth 6 DF))
	 (d-n3         (+ d-n3-tl-31 d-n3-tl-22))
	 (d-n1-sl-top  (nth 0 DB))
	 (d-n3-22-top  (nth 1 DB))
	 (d-n3-31-top  (nth 2 DB))
	 (d-n1-sl-bot  (nth 3 DB))
	 (d-n3-22-bot  (nth 4 DB))
	 (d-n3-31-bot  (nth 5 DB))

	 ;; N3 is not defined as a global constant. So, for
	 ;; convenience, it is defined here.  
	 ;; JM: There's a macro for this. 
	 ;; However, I use the let statement here for clarity.
	 (N3 (+ N3-TL-31 N3-TL-22))

	 ;; Determine the change in the damping
	 (delta-damping (- (damping (+ N3 d-n3)) (damping N3)))
	 
	 ;; Determine the change in the action
	 (old-action (action N1-SL N1-TL N3-TL-31 N3-TL-22 
			     *N1-SL-TOP* *N3-22-TOP* *N3-31-TOP*
			     *N1-SL-BOT* *N3-22-BOT* *N3-31-BOT*))
	 (new-action (action (+ N1-SL d-n1-sl)
			     (+ N1-TL d-n1-tl)
			     (+ N3-TL-31 d-n3-tl-31)
			     (+ N3-TL-22 d-n3-tl-22)
			     (+ *N1-SL-TOP* d-n1-sl-top)
			     (+ *N3-22-TOP* d-n3-22-top)
			     (+ *N3-31-TOP* d-n3-31-top)
			     (+ *N1-SL-BOT* d-n1-sl-bot)
			     (+ *N3-22-BOT* d-n3-22-bot)
			     (+ *N3-31-BOT* d-n3-31-bot)))

	 ;; Complete the Wick rotation by multiplying by i to form the
	 ;; Euclidean action.
	 (delta-action (* *i* (- new-action old-action)))
	 
	 ;; We need the action to be real and <1 to get a real
	 ;; probability distribution out:
	 (action-okay (zerop (imagpart delta-action))))

    ;; Raise an error message if the action is not okay.
    ;; (print delta-action) ; For debugging the action. 
    ;; (print (N3)) ; For debugging the action
    (if (not action-okay)
	(prog nil
	   (print "Data:")
	   (print (list d-n1-sl d-n1-tl d-n3-tl-31 d-n3-tl-22
			d-n1-sl-top d-n3-22-top d-n3-31-top
			d-n1-sl-bot d-n3-22-bot d-n3-31-bot
			*alpha* *k* *litL*))
	   (print "i*action:")
	   (print delta-action)
	   (error "Action must be completely imaginary."))
	(< (random 1.0) (* (exp (realpart delta-action))
			   (exp (- delta-damping)))))))

;; JM: The old accept-move? function. Left for completeness. Doesn't
;; take into account the additional terms to the action. DEPRECATED
;; USE AT YOUR OWN RISK.
(defun accept-move?-deprecated (mtype)
  (let ((delta-action 0.0)
	(delta-damping 0.0))
    (cond ((= 0 mtype) ;; 2->6 move
	   (setf delta-damping (- (damping (+ (N3) 4)) (damping (N3))))
	   (setf delta-action 
		 (- (action (+ N1-SL 3) (+ N1-TL 2) 
			    (+ N3-TL-31 4) (+ N3-TL-22 0)) 
		    (action N1-SL N1-TL N3-TL-31 N3-TL-22)))) 
	  ((= 1 mtype) ;; 2->3 move
	   (setf delta-damping (- (damping (+ (N3) 1)) (damping (N3))))
	   (setf delta-action (- (action (+ N1-SL 0) (+ N1-TL 1) (+ N3-TL-31 0) (+ N3-TL-22 1)) 
				 (action N1-SL N1-TL N3-TL-31 N3-TL-22))))
	  ((= 2 mtype) ;; 4->4 move
	   (setf delta-damping (- (damping (+ (N3) 0)) (damping (N3))))
	   (setf delta-action (- (action (+ N1-SL 0) (+ N1-TL 0) (+ N3-TL-31 0) (+ N3-TL-22 0)) 
				 (action N1-SL N1-TL N3-TL-31 N3-TL-22))))
	  ((= 3 mtype) ;; 3->2 move
	   (setf delta-damping (- (damping (+ (N3) -1)) (damping (N3))))
	   (setf delta-action (- (action (+ N1-SL 0) (+ N1-TL -1) (+ N3-TL-31 0) (+ N3-TL-22 -1)) 
				 (action N1-SL N1-TL N3-TL-31 N3-TL-22))))
	  ((= 4 mtype) ;; 6->2 move
	   (setf delta-damping (- (damping (+ (N3) -4)) (damping (N3))))
	   (setf delta-action (- (action (+ N1-SL -3) (+ N1-TL -2) (+ N3-TL-31 -4) (+ N3-TL-22 0)) 
				 (action N1-SL N1-TL N3-TL-31 N3-TL-22)))))
    (< (random 1.0) (* (exp (realpart (* *i* delta-action)))
		       (exp (* -1.0 delta-damping))))))


;; a sweep is defined as N-INIT number of attempted moves 

;; JM: Because accept-move? now requires sxid as an input, sweep is a
;; slightly different function than it was before.
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
	  ;; (printmove mtype) ; for debugging purposes.
	  (incf (nth mtype SUCCESSFUL-MOVES)) ;; number of moves of mtype that have succeeded
	  (2plus1move movedata))))))

;; following is to be used for tuning the k0 and k3 parameters
(defun generate-data-console (&optional (start-sweep 1))
  (for (ns start-sweep (+ start-sweep NUM-SWEEPS -1))
       (sweep)
;;       (format t "~%Top Boundary: ~a~%Bottom Boundary: ~a~%" (length (topfaces)) (length (bottomfaces))) ; For  debugging
       (when (= 0 (mod ns 10))
	 (format t "start = ~A end = ~A current = ~A count = ~A ~A\%~%"
		 start-sweep (+ start-sweep NUM-SWEEPS -1) ns 
		 ;(count-simplices-of-all-types) (percent-tv)));SUCCESSFUL-MOVES));(accept-ratios)))
		 (count-simplices-of-all-types) (accept-ratios)))
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
	     (format progfile "~A/~A/~A ~A~%"
		     start-sweep ns end-sweep (count-simplices-of-all-types)))))))

;; generate-data-v2 is similar to generate-data except it creates a fresh data file every 
;; SAVE-EVERY-N-SWEEPS. since a fresh datafile is created, there is no need to maintain a seprate progress
;; file.
(defun generate-data-v2 (&optional (start-sweep 1))
  (setf SIM-START-TIME (cdt-now-str))
  (let ((end-sweep (+ start-sweep NUM-SWEEPS -1)))
    (when (= 1 start-sweep) ;; save the initial spacetime contents if this is a brand new run
      (with-open-file (datafile (concatenate 'string (generate-filename-v2 start-sweep 0) 3SXEXT)
				:direction :output
				:if-exists :supersede)
	(save-spacetime-to-file datafile)))
    (for (ns start-sweep end-sweep)
	 (sweep)
	 (when (= 0 (mod ns SAVE-EVERY-N-SWEEPS))
	   (with-open-file (datafile (concatenate 'string (generate-filename-v2 start-sweep ns) 3SXEXT)
				     :direction :output
				     :if-exists :supersede)
	     (save-spacetime-to-file datafile))))))

;; generate-data-v3 is similar to generate-data-v2 except it also creates an 
;; additional data file every 
;; SAVE-EVERY-N-SWEEPS that contains the spatial 2-simplex information for 
;; each spatial slice.
(defun generate-data-v3 (&optional (start-sweep 1))
  (setf SIM-START-TIME (cdt-now-str))
  (let ((end-sweep (+ start-sweep NUM-SWEEPS -1)))
    (when (= 1 start-sweep) ;; save the initial spacetime contents if this is a brand new run
      (with-open-file 
	  (datafile 
	   (concatenate 'string (generate-filename-v2 start-sweep 0) 3SXEXT)
	   :direction :output
	   :if-exists :supersede)
	(save-spacetime-to-file datafile)))
    (for (ns start-sweep end-sweep)
	 (sweep)
	 (when (= 0 (mod ns SAVE-EVERY-N-SWEEPS))
	   (let ((filename (generate-filename-v2 start-sweep ns)))
	     (with-open-file (datafile (concatenate 'string filename 3SXEXT)
				       :direction :output
				       :if-exists :supersede)
	       (save-spacetime-to-file datafile))
	     (3sx2p1->s2sx2p1)
	     (with-open-file (datafile (concatenate 'string filename S2SXEXT)
				       :direction :output
				       :if-exists :supersede)
	       (save-s2simplex-data-to-file datafile)))))))

;; generate-movie-data saves number of simplices every SAVE-EVERY-N-SWEEPS
(defun generate-movie-data (&optional (start-sweep 1))
  (setf SAVE-EVERY-N-SWEEPS 10)
  (let ((moviefilestr 
	 (concatenate 'string (generate-filename start-sweep) MOVEXT))
	(trackfilestr 
	 (concatenate 'string (generate-filename start-sweep) PRGEXT))
	(end-sweep (+ start-sweep NUM-SWEEPS -1))
	(max-ts (if (string= BCTYPE "OPEN") (1- NUM-T) (- NUM-T 2))))

    ;; open and close the file for :append to work
    (with-open-file (moviefile moviefilestr 
			       :direction :output
			       :if-exists :supersede)
      ;; record the initial data only if start-sweep = 1
      (when (= start-sweep 1)

	(for (ts 0 max-ts)
	     (format moviefile "~A " (count-simplices-in-sandwich ts (1+ ts))))
	(format moviefile "~%")))
    
    (for (ns start-sweep end-sweep)
	 (sweep)
	 (when (= 0 (mod ns SAVE-EVERY-N-SWEEPS))
	   (with-open-file (moviefile moviefilestr 
				      :direction :output
				      :if-exists :append)
	     (for (ts 0 max-ts)
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
	 (datafilestr (concatenate 'string (generate-filename start-sweep) 3SXEXT))
	 (trackfilestr (concatenate 'string (generate-filename start-sweep) PRGEXT))
	 (moviefilestr (concatenate 'string (generate-filename start-sweep) MOVEXT))
	 (ts-max (if (string= BCTYPE "OPEN") (1- NUM-T) (- NUM-T 2))))
    
    ;; open and close the file, for :append to work properly
    (with-open-file (moviefile moviefilestr 
			       :direction :output
			       :if-exists :supersede)
      ;; record the initial data only if start-sweep = 1
      (when (= 1 start-sweep)
	(for (ts 0 ts-max)
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
	     (for (ts 0 ts-max)
		  (format moviefile "~A " (count-simplices-in-sandwich ts (1+ ts))))
	     (format moviefile "~%"))))))

#|
(defun calculate-order-parameter (&optional (start-sweep 1))
  (let* ((end-sweep (+ start-sweep NUM-SWEEPS -1))
	 (order-parameter 0.0)
	 (datafilestr (format nil 
			      "~A-~A-op-T~A_V~A_eps~A_kz~A_kt~A_sweeps~Ato~A.op" 
			      *topology* *boundary-conditions*
			      NUM-T N-INIT *eps* *k0* *k3* start-sweep end-sweep))
	 (trackfilestr (format nil 
			       "~A-~A-op-T~A_V~A_eps~A_kz~A_kt~A_sweeps~Ato~A.progress" 
			       *topology* *boundary-conditions*
			       NUM-T N-INIT *eps* *k0* *k3* start-sweep end-sweep)))
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
	      NUM-T N-INIT *eps* *k0* *k3* start-sweep end-sweep order-parameter))))

(defun calculate-volume-volume-correlator (&optional (start-sweep 1))
  (let* ((end-sweep (+ start-sweep NUM-SWEEPS -1))
	 (vvparams (make-array (+ NUM-T 1) :initial-element 0.0))
	 (dfilestr (format nil 
			   "~A-~A-vv-T~A_V~A_eps~A_kz~A_kt~A_sweeps~Ato~A.vv" 
			   *topology* *boundary-conditions*
			   NUM-T N-INIT *eps* *k0* *k3* start-sweep end-sweep))
	 (tfilestr (format nil 
			   "~A-~A-vv-T~A_V~A_eps~A_kz~A_kt~A_sweeps~Ato~A.prog" 
			   *topology* *boundary-conditions*
			   NUM-T N-INIT *eps* *k0* *k3* start-sweep end-sweep)))
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
	      NUM-T N-INIT *eps* *k0* *k3* start-sweep end-sweep vvparams))))

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
