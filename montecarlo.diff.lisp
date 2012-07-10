2a3
> ;; translate a move type into an actual move
4,9c5,13
<   (ecase mtype
<     (0 (try-2->6 sxid))
<     (1 (try-2->3 sxid))
<     (2 (try-4->4 sxid))
<     (3 (try-3->2 sxid))
<     (4 (try-6->2 sxid))))
---
>   (let ((sx (get-3simplex sxid)))
>     (if (and sx (is-real-simplex sx))
> 	(ecase mtype
> 	  (0 (try-2->6 sxid))
> 	  (1 (try-2->3 sxid))
> 	  (2 (try-4->4 sxid))
> 	  (3 (try-3->2 sxid))
> 	  (4 (try-6->2 sxid)))
> 	nil)))
27,53c31,103
< (defun accept-move? (mtype)
<   (let ((delta-action 0.0)
< 	(delta-damping 0.0))
<     (cond ((= 0 mtype) ;; 2->6 move
< 	   (setf delta-damping (- (damping (+ (N3) 4)) (damping (N3))))
< 	   (setf delta-action 
< 		 (- (action (+ N1-SL 3) (+ N1-TL 2) 
< 			    (+ N3-TL-31 4) (+ N3-TL-22 0)) 
< 		    (action N1-SL N1-TL N3-TL-31 N3-TL-22)))) 
< 	  ((= 1 mtype) ;; 2->3 move
< 	   (setf delta-damping (- (damping (+ (N3) 1)) (damping (N3))))
< 	   (setf delta-action (- (action (+ N1-SL 0) (+ N1-TL 1) (+ N3-TL-31 0) (+ N3-TL-22 1)) 
< 				 (action N1-SL N1-TL N3-TL-31 N3-TL-22))))
< 	  ((= 2 mtype) ;; 4->4 move
< 	   (setf delta-damping (- (damping (+ (N3) 0)) (damping (N3))))
< 	   (setf delta-action (- (action (+ N1-SL 0) (+ N1-TL 0) (+ N3-TL-31 0) (+ N3-TL-22 0)) 
< 				 (action N1-SL N1-TL N3-TL-31 N3-TL-22))))
< 	  ((= 3 mtype) ;; 3->2 move
< 	   (setf delta-damping (- (damping (+ (N3) -1)) (damping (N3))))
< 	   (setf delta-action (- (action (+ N1-SL 0) (+ N1-TL -1) (+ N3-TL-31 0) (+ N3-TL-22 -1)) 
< 				 (action N1-SL N1-TL N3-TL-31 N3-TL-22))))
< 	  ((= 4 mtype) ;; 6->2 move
< 	   (setf delta-damping (- (damping (+ (N3) -4)) (damping (N3))))
< 	   (setf delta-action (- (action (+ N1-SL -3) (+ N1-TL -2) (+ N3-TL-31 -4) (+ N3-TL-22 0)) 
< 				 (action N1-SL N1-TL N3-TL-31 N3-TL-22)))))
<     (< (random 1.0) (* (exp (realpart (* *i* delta-action)))
< 		       (exp (* -1.0 delta-damping))))))
---
> (defun accept-move? (mtype sxid)
> 
>   (let* (;;determine what change vectors to use, based on the type of move
> 	
> 	 ;;this is the change in the f-vector due to the move (see the DF*
> 	 ;;parameters and/or macros in "globals.lisp")
> 	 (DF
> 	  (ecase mtype
> 	    (0;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; 2->6
> 	     (DF26 sxid))
> 	    (1;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; 2->3
> 	     DF23)
> 	    (2;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; 4->4
> 	     DF44)
> 	    (3;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; 3->2
> 	     DF32)
> 	    (4;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; 6->2
> 	     (DF62 sxid))))
> 	 
> 	 ;;this is the change in some relevant quantities at the boundary
> 	 ;;(see the macros/paramters DB* in "globals.lisp")
> 	 (DB
> 	  (ecase mtype
> 	    (0;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; 2->6
> 	     (DB26 sxid))
> 	    (1;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; 2->3
> 	     (DB23 sxid))
> 	    (2;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; 4->4
> 	     DB44)
> 	    (3;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; 3->2
> 	     (DB32 sxid))
> 	    (4;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; 6->2
> 	     (DB62 sxid))))
> 		 
> 	 ;;determine changes in different numbers of geometrical objects
> 	 (d-n3       (+ (seventh DF) (sixth DF)))
> 	 (d-n1-sl    (second  DF))
> 	 (d-n3-tl-22 (seventh DF))
> 	 (d-n3-tl-31 (sixth   DF))
> 	 (d-n1-tl    (third   DF))
> 	 (d-n1-sl-b  (first   DB))
> 	 (d-n3-22-b  (second  DB))
> 	 (d-n3-31-b  (third   DB))
> 	 
> 	 ;;determine deltas in damping and action
> 	 (delta-damping (- (damping (+ (N3) d-n3)) (damping (N3))))
> 	 (delta-action  (action  d-n1-sl   d-n1-tl d-n3-tl-31 d-n3-tl-22 
> 				 d-n1-sl-b 
> 				 d-n3-22-b 
> 				 d-n3-31-b
> 				 *alpha* *k* *litl*))
> 	 ;; Determine whether the action is real or imaginary
> 	 (action-is-imaginary (zerop (realpart delta-action))))
>     
>     ;;accept with probability as function of changes in action and damping
> 
>     ;;this expression assumes that the "delta-action" is the change in
>     ;;real-valued euclidean action
>     (when (not action-is-imaginary)
> 	(prog nil
> 	   (print "Data:")
> 	   (print (list d-n1-sl d-n1-tl d-n3-tl-31 d-n3-tl-22 
> 			d-n1-sl-b d-n3-22-b d-n3-31-b *alpha* *k* *litl*))
> 	   (print "Action:")
> 	   (print (action  d-n1-sl   d-n1-tl d-n3-tl-31 d-n3-tl-22 
> 				 d-n1-sl-b 
> 				 d-n3-22-b 
> 				 d-n3-31-b
> 				 *alpha* *k* *litl*))
> 	   (error "Action must be completely imaginary. Something is wrong."))
> 	(< (random 1.0) (* (exp (realpart (* *i* delta-action)))
> 			   (exp (- delta-damping)))))))
> 
62c112
< 	(while (null movedata) 
---
> 	(while (null movedata)
68c118
< 	(when (accept-move? mtype)
---
> 	(when (accept-move? mtype sxid)
77c127,133
< 	 (format t "start = ~A end = ~A current = ~A count = ~A ~A\%~%"
---
> 
> ;	 (format t "start = ~A end = ~A current = ~A count = ~A ~A\%~%"
> ;		 start-sweep (+ start-sweep NUM-SWEEPS -1) ns 
> ;		 (count-simplices-of-all-types) (percent-tv)));SUCCESSFUL-MOVES));(accept-ratios)))
> ;		 ;(count-simplices-of-all-types) (accept-ratios)))
> 
> 	 (format t "start = ~A end = ~A current = ~A count = ~A ~$\%~%"
79,80c135,136
< 		 ;(count-simplices-of-all-types) (percent-tv)));SUCCESSFUL-MOVES));(accept-ratios)))
< 		 (count-simplices-of-all-types) (accept-ratios)))
---
> 		 (count-boundary-vs-bulk) (percent-tv)))
> 
85a142
>   (setf SIM-START-TIME (cdt-now-str))
99,100c156,157
< 	     (format progfile "~A/~A/~A ~A~%"
< 		     start-sweep ns end-sweep (count-simplices-of-all-types)))))))
---
> 	     (format progfile "start = ~A end = ~A current = ~A count = ~A~%"
> 		     start-sweep end-sweep ns (count-simplices-of-all-types)))))))
102c159
< ;; generate-data-v2 is similar to generate-data except it creates a fresh data file every 
---
> ;; generate-data-v2 is similar to generate data except it creates a fresh data file every 
108c165
<     (when (= 1 start-sweep) ;; save the initial spacetime contents if this is a brand new run
---
>     (when (= 1 start-sweep)
112,144c169,173
< 	(save-spacetime-to-file datafile)))
<     (for (ns start-sweep end-sweep)
< 	 (sweep)
< 	 (when (= 0 (mod ns SAVE-EVERY-N-SWEEPS))
< 	   (with-open-file (datafile (concatenate 'string (generate-filename-v2 start-sweep ns) 3SXEXT)
< 				     :direction :output
< 				     :if-exists :supersede)
< 	     (save-spacetime-to-file datafile))))))
< 
< ;; generate-data-v3 is similar to generate-data-v2 except it also creates an 
< ;; additional data file every 
< ;; SAVE-EVERY-N-SWEEPS that contains the spatial 2-simplex information for 
< ;; each spatial slice.
< (defun generate-data-v3 (&optional (start-sweep 1))
<   (setf SIM-START-TIME (cdt-now-str))
<   (let ((end-sweep (+ start-sweep NUM-SWEEPS -1)))
<     (when (= 1 start-sweep) ;; save the initial spacetime contents if this is a brand new run
<       (with-open-file 
< 	  (datafile 
< 	   (concatenate 'string (generate-filename-v2 start-sweep 0) 3SXEXT)
< 	   :direction :output
< 	   :if-exists :supersede)
< 	(save-spacetime-to-file datafile)))
<     (for (ns start-sweep end-sweep)
< 	 (sweep)
< 	 (when (= 0 (mod ns SAVE-EVERY-N-SWEEPS))
< 	   (let ((filename (generate-filename-v2 start-sweep ns)))
< 	     (with-open-file (datafile (concatenate 'string filename 3SXEXT)
< 				       :direction :output
< 				       :if-exists :supersede)
< 	       (save-spacetime-to-file datafile))
< 	     (3sx2p1->s2sx2p1)
< 	     (with-open-file (datafile (concatenate 'string filename S2SXEXT)
---
> 	(save-spacetime-to-file datafile))
>       (for (ns start-sweep end-sweep)
> 	   (sweep)
> 	   (when (= 0 (mod ns SAVE-EVERY-N-SWEEPS))
> 	     (with-open-file (datafile (concatenate 'string (generate-filename-v2 start-sweep ns) 3SXEXT)
147c176
< 	       (save-s2simplex-data-to-file datafile)))))))
---
> 	       (save-spacetime-to-file datafile)))))))
152,155c181,182
<   (let ((moviefilestr 
< 	 (concatenate 'string (generate-filename start-sweep) MOVEXT))
< 	(trackfilestr 
< 	 (concatenate 'string (generate-filename start-sweep) PRGEXT))
---
>   (let ((moviefilestr (concatenate 'string (generate-filename start-sweep) MOVEXT))
> 	(trackfilestr (concatenate 'string (generate-filename start-sweep) PRGEXT))
181,182c208,209
< 	     (format trackfile "~A/~A/~A ~A~%"
< 		     start-sweep ns end-sweep (count-simplices-of-all-types)))))))
---
> 	     (format trackfile "start = ~A end = ~A current = ~A count = ~A~%"
> 		     start-sweep end-sweep ns (count-simplices-of-all-types)))))))
231a259,271
> (defun 3sx2p1->2sx2p1 (infile outfile)
>   "3sx2p1->2sx2p1 generates the 2-simplex information for each spatial slice from the 3-simplex data for
> the entire spacetime. The generated information is written to outfile"
>   (load-spacetime-from-file infile)
>   (clrhash *ID->SPATIAL-2SIMPLEX*)
>   (for (ts 0 (1- NUM-T))
>        (let ((31simplices (get-simplices-in-sandwich-of-type ts (1+ ts) 3))
> 	     (spatial-triangles '()))
> 	 (dolist (31simplex 31simplices)
> 	   (push (make-s2simplex ts (3sx-lopts (get-3simplex 31simplex))) spatial-triangles))
> 	 (connect-spatial-2simplices-within-list spatial-triangles)))
>   (save-s2simplex-data-to-file outfile))
> 
239c279
< 			      NUM-T N-INIT *eps* *k0* *k3* start-sweep end-sweep))
---
> 			      NUM-T N-INIT eps k0 k3 start-sweep end-sweep))
243c283
< 			       NUM-T N-INIT *eps* *k0* *k3* start-sweep end-sweep)))
---
> 			       NUM-T N-INIT eps k0 k3 start-sweep end-sweep)))
258c298
< 	      NUM-T N-INIT *eps* *k0* *k3* start-sweep end-sweep order-parameter))))
---
> 	      NUM-T N-INIT eps k0 k3 start-sweep end-sweep order-parameter))))
266c306
< 			   NUM-T N-INIT *eps* *k0* *k3* start-sweep end-sweep))
---
> 			   NUM-T N-INIT eps k0 k3 start-sweep end-sweep))
270c310
< 			   NUM-T N-INIT *eps* *k0* *k3* start-sweep end-sweep)))
---
> 			   NUM-T N-INIT eps k0 k3 start-sweep end-sweep)))
293c333
< 	      NUM-T N-INIT *eps* *k0* *k3* start-sweep end-sweep vvparams))))
---
> 	      NUM-T N-INIT eps k0 k3 start-sweep end-sweep vvparams))))
313c353
< |#
\ No newline at end of file
---
> |#
