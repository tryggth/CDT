; cdt-2plus1-montecarlo.lisp

;; Module for cdt in 2+1 dimensions. Runs the actual monte carlo sweeps.

;; Authors:
;; -------- Rajesh Kommu
;; -------- David Kemansky
;; -------- Jonah Miller (jonah.maxwell.miller@gmail.com)

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
	 (format t "finished ~A of ~A sweeps with count ~A~%" 
		 sweepnum nsweeps 
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
	    (0 (DB26 sxid))   ; 2->6 move
	    (1 (DB23 sxid))   ; 2->3 move
	    (2 DB44)          ; 4->4 move
	    (3 (DB32 sxid))   ; 3->2 move
	    (4 (DB62 sxid)))) ; 6->2 move 
	 
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
	 (new-action (action (+ N1-SL       d-n1-sl)
			     (+ N1-TL       d-n1-tl)
			     (+ N3-TL-31    d-n3-tl-31)
			     (+ N3-TL-22    d-n3-tl-22)
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

    ;; DEBUGGING!
;;    (format t "Change in V: ~A~%Change in action: ~A~%Change in damping: ~A~%Probability of acceptance: ~A~%" 
;;	    d-n3 delta-action delta-damping (* (exp (realpart delta-action))
;;					       (exp (- delta-damping))))
;;    (format t "Change in V: ~A~%P of acceptance: ~A~%" d-n3 (* (exp (realpart delta-action)) (exp (- delta-damping))))

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
;(defun accept-move?-deprecated (mtype)
;  (let ((delta-action 0.0)
;	(delta-damping 0.0))
;    (cond ((= 0 mtype) ;; 2->6 move
;	   (setf delta-damping (- (damping (+ (N3) 4)) (damping (N3))))
;	   (setf delta-action 
;		 (- (action (+ N1-SL 3) (+ N1-TL 2) 
;			    (+ N3-TL-31 4) (+ N3-TL-22 0)) 
;		    (action N1-SL N1-TL N3-TL-31 N3-TL-22)))) 
;	  ((= 1 mtype) ;; 2->3 move
;	   (setf delta-damping (- (damping (+ (N3) 1)) (damping (N3))))
;	   (setf delta-action (- (action (+ N1-SL 0) (+ N1-TL 1) 
;					 (+ N3-TL-31 0) (+ N3-TL-22 1)) 
;				 (action N1-SL N1-TL N3-TL-31 N3-TL-22))))
;	  ((= 2 mtype) ;; 4->4 move
;	   (setf delta-damping (- (damping (+ (N3) 0)) (damping (N3))))
;	   (setf delta-action (- (action (+ N1-SL 0) (+ N1-TL 0) 
;					 (+ N3-TL-31 0) (+ N3-TL-22 0)) 
;				 (action N1-SL N1-TL N3-TL-31 N3-TL-22))))
;	  ((= 3 mtype) ;; 3->2 move
;	   (setf delta-damping (- (damping (+ (N3) -1)) (damping (N3))))
;	   (setf delta-action (- (action (+ N1-SL 0) (+ N1-TL -1) 
;					 (+ N3-TL-31 0) (+ N3-TL-22 -1)) 
;				 (action N1-SL N1-TL N3-TL-31 N3-TL-22))))
;	  ((= 4 mtype) ;; 6->2 move
;	   (setf delta-damping (- (damping (+ (N3) -4)) (damping (N3))))
;	   (setf delta-action (- (action (+ N1-SL -3) (+ N1-TL -2) 
;					 (+ N3-TL-31 -4) (+ N3-TL-22 0)) 
;				 (action N1-SL N1-TL N3-TL-31 N3-TL-22)))))
 ;   (< (random 1.0) (* (exp (realpart (* *i* delta-action)))
;		       (exp (* -1.0 delta-damping))))))


;; a sweep is defined as N-INIT number of attempted moves 

;; JM: Because accept-move? now requires sxid as an input, sweep is a
;; slightly different function than it was before.
(defun sweep ()
  "N3 iterations of the Metropolis-Hastings algorithm."
  (let ((num-attempted 0))
    (while (< num-attempted N-INIT)
      (let* ((sxid (random *LAST-USED-3SXID*))
	     (mtype (select-move))
	     (movedata (try-move sxid mtype)))
	(while (null movedata) 
	  (setf sxid (random *LAST-USED-3SXID*)
		mtype (select-move) 
		movedata (try-move sxid mtype)))
	(incf num-attempted) ;; number-of-attempted-moves-counter for
			     ;; this sweep
	(incf (nth mtype ATTEMPTED-MOVES)) ;; number of moves of mtype
					   ;; that have been attempted
	(when (accept-move? mtype sxid)
	  ;; (printmove mtype) ; for debugging purposes.
	  (incf (nth mtype SUCCESSFUL-MOVES)) ;; number of moves of
					      ;; mtype that have
					      ;; succeeded
	  (2plus1move movedata))))))

