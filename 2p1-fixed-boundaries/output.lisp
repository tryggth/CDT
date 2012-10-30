;; output.lisp
;; Author: Jonah Miller (jonah.maxwell.miller@gmail.com)
;; Date: August 5, 2012

;; This program is a primary module for the 2+1-dimensional
;; fixed-boundary CDT simulations. It replaces the later part of
;; 2p1/montecarlo.lisp and streamlines the generate-*-data functions that
;; actually produce output.

;; Requires the following modules to be loaded first:

;; ../utilities.lisp
;; globals.lisp
;; generalized-hash-table-counting-functions.lisp
;; simplex.lisp
;; moves.lisp
;; montecarlo.lisp



;;;; Functions that make file names:
;;;;-------------------------------------------------------------------------
;; When running scripts, set this to change the default output directory.
(defparameter *output-directory* ""
  "The directory files are saved to. By default, uses ./")

;; STOPOLOGY-BCTYPE-NUMT-NINIT-k0-k3-eps-alpha-startsweep-endsweep-initialBoundary-finalBoundary-hostname-starttime
(defun generate-filename (&optional (start-sweep 1) 
			    (end-sweep (+ start-sweep NUM-SWEEPS -1))
			    (signifier nil))
  (if signifier
      (concatenate 'string *output-directory*
		   (format nil "~A-~A-T~3,'0d-V~6,'0d-~A-~A-~A-~A-~9,'0d-~9,'0d-IB0~A-FB0~A-S0~A-on-~A-started~A" 
			   STOPOLOGY BCTYPE NUM-T N-INIT
			   *k0* *k3* *eps* *alpha*
			   start-sweep end-sweep
			   (count-spacelike-triangles-at-time 0)
			   (count-spacelike-triangles-at-time NUM-T)
			   signifier
			   (hostname) (cdt-now-str)))
      (concatenate 'string *output-directory*
		   (format nil "~A-~A-T~3,'0d-V~6,'0d-~A-~A-~A-~A-~9,'0d-~9,'0d-IB0~A-FB0~A-on-~A-started~A" 
			   STOPOLOGY BCTYPE NUM-T N-INIT
			   *k0* *k3* *eps* *alpha*
			   start-sweep end-sweep
			   (count-spacelike-triangles-at-time 0)
			   (count-spacelike-triangles-at-time NUM-T)
			   (hostname) (cdt-now-str)))))

;; STOPOLOGY-BCTYPE-NUMT-NINIT-k0-k3-eps-alpha-INITIALBOUNDARY-FINALBOUNDARY-startsweep-currsweep-endsweep-initialBoundary-finalBoundary-hostname-starttime-currenttime
(defun generate-filename-v2 (initial-boundary final-boundary
			     &optional (ssweep 1) (csweep 0) 
			       (esweep (+ ssweep NUM-SWEEPS -1))
			       (signifier nil))
  (if signifier
      (concatenate 'string *output-directory*
		   (format nil "~A-~A-T~3,'0d-V~6,'0d-~A-~A-~A-~A-~9,'0d-~9,'0d-~9,'0d-IB0~A-FB0~A-S0~A-on-~A-start~A-curr~A" 
			   STOPOLOGY BCTYPE NUM-T N-INIT
			   *k0* *k3* *eps* *alpha* 
			   ssweep csweep esweep
			   initial-boundary final-boundary
			   signifier
			   (hostname) SIM-START-TIME (cdt-now-str)))
      (concatenate 'string *output-directory*
		   (format nil "~A-~A-T~3,'0d-V~6,'0d-~A-~A-~A-~A-~9,'0d-~9,'0d-~9,'0d-IB0~A-FB0~A-on-~A-start~A-curr~A" 
			   STOPOLOGY BCTYPE NUM-T N-INIT
			   *k0* *k3* *eps* *alpha* 
			   ssweep csweep esweep
			   initial-boundary final-boundary
			   (hostname) SIM-START-TIME (cdt-now-str)))))

;; STOPOLOGY-BCTYPE-NUMT-NINIT-k0-k3-eps-alpha-INITIALBOUNDARY-FINALBOUNDARY-hostname-starttime
(defun generate-filename-v3 (&optional (signifier nil))
  (if signifier
      (concatenate 'string *output-directory*
		   (format
		    nil
		    "~A-~A-T~3,'0d-V~6,'0d-~A-~A-~A-~A-IB0~A-FB0~A-S0~A-on-~A-start~A" 
		    STOPOLOGY BCTYPE NUM-T N-INIT *k0* *k3* *eps* *alpha*
		    (count-spacelike-triangles-at-time 0)
		    (count-spacelike-triangles-at-time NUM-T)
		    signifier
		    (hostname) SIM-START-TIME))
            (concatenate 'string *output-directory*
		   (format
		    nil
		    "~A-~A-T~3,'0d-V~6,'0d-~A-~A-~A-~A-IB0~A-FB0~A-on-~A-start~A" 
		    STOPOLOGY BCTYPE NUM-T N-INIT *k0* *k3* *eps* *alpha*
		    (count-spacelike-triangles-at-time 0)
		    (count-spacelike-triangles-at-time NUM-T)
		    (hostname) SIM-START-TIME))))
      

;;;;-------------------------------------------------------------------------


;;;; Functions which print data to a specified iostream.
;;;;-------------------------------------------------------------------------
(defun print-movie-data (iostream)
  "Prints one frame of movie data to iostream."
  (let ((t-max (if (string= BCTYPE "OPEN") (1- NUM-T) (- NUM-T 2))))
    ; The simplex count for all but the last time slice has white
    ; space after it.
    (for (ts 0 (1- t-max))
      (format iostream "~A " (count-simplices-in-sandwich ts (1+ ts))))
    ; Instead of a space after it, the last simplex count has a
    ; newline character.
    (format iostream "~A~%" (count-simplices-in-sandwich t-max (1+ t-max)))))
;;;;-------------------------------------------------------------------------



;;;; Functions which make new files
;;;;-------------------------------------------------------------------------
(defun make-spacetime-file (filename)
  "A slightly more convenient wrapper for save-spacetime-to-file. The
   filename input is made by (generate-filename start-sweep) or somthing
   similar. It does not have the file extension."
  (with-open-file (datafile (concatenate 'string filename 3SXEXT)
			    :direction :output
			    :if-exists :supersede)
    (save-spacetime-to-file datafile)))

(defun make-progress-file (filename start-sweep ns end-sweep)
  "Makes a progress file that tells us the start sweep, the current
   sweep, the end sweep, and the numbers of simplices of all
   types. Mostly useful for people who want to know the progress of a
   simulation. However, it also restores the program in the case of a
   catastrophic crash. The filename input is made
   by (generate-filename start-sweep) or somthing similar. It does not
   have the file extension."
  (with-open-file (progfile (concatenate 'string filename PRGEXT)
			    :direction :output
			    :if-exists :supersede)
    (format progfile "~A/~A/~A ~A~%"
	    start-sweep ns end-sweep (count-simplices-of-all-types))))

(defun make-movie-file (filename start-sweep)
  "Make a movie file and add data to it if start-sweep = 1. The
   filename input is made by (generate-filename start-sweep) or
   somthing similar. It does not have the file extension."
  (with-open-file (moviefile (concatenate 'string filename MOVEXT)
					  :direction :output
					  :if-exists :supersede)
    (when (= 1 start-sweep)
      (print-movie-data moviefile))))
;;;;-------------------------------------------------------------------------



;;;; Functions which append data to a file
;;;;-------------------------------------------------------------------------
(defun append-to-movie-file (filename)
  "Appends movie data to a file made by make-movie-file above. The
   filename input is made by (generate-filename start-sweep) or
   somthing similar. It does not have the file extension."
  (with-open-file (moviefile (concatenate 'string filename MOVEXT)
			     :direction :output
			     :if-exists :append)
    (print-movie-data moviefile)))
;;;;-------------------------------------------------------------------------



;;;; Functions which generate data
;;;;-------------------------------------------------------------------------
;; The following function is to be used for tuning the k0 and k3 parameters
(defun generate-data-console (&optional (start-sweep 1))
  "Prints a number of useful quantities to the console as sweeps are
   performed. Not useful for long simulations."
  (for (ns start-sweep (+ start-sweep NUM-SWEEPS -1))
       (sweep)
       (format t "~%Top Boundary: ~a~%Bottom Boundary: ~a~%" 
	       (length (topfaces)) (length (bottomfaces)))
       (when (= 0 (mod ns 10))
	 (format t "start = ~A end = ~A current = ~A count = ~A ~A\%~%"
		 start-sweep (+ start-sweep NUM-SWEEPS -1) ns 
		 (count-simplices-of-all-types) (accept-ratios)))
       (finish-output)))

;; generate-data should be called after setting the values for eps,
;; k0, k3, NUM-SWEEPS and calling one of the initialize-xx-slices. It
;; generates a single datafile and a single progress file after every
;; save-every-n-sweeps.
(defun generate-data (&optional (start-sweep 1) (signifier nil))
  "Generates a single datafile and a single progress file after every
   SAVE-EVERY-N-SWEEPS."
  (let* ((end-sweep (+ start-sweep NUM-SWEEPS -1))
	 (filename (generate-filename start-sweep end-sweep signifier)))
	
    (for (ns start-sweep end-sweep)
	 (sweep)
	 (when (= 0 (mod ns SAVE-EVERY-N-SWEEPS))
	   (make-spacetime-file filename)
	   (make-progress-file filename start-sweep ns end-sweep)))))

;; generate-data-in-time should be called after setting the values for
;; eps, k0, k3, NUM-SWEEPS and calling one of the
;; initialize-xx-slices. It generates a single data file and a
;; progress file after every save-every-n-sweeps and runs for runtime
;; seconds. You don't need to set NUM-SWEEPS, but it might be useful
;; to do so to keep track of how many sweeps you hope will be
;; finished.
(defun generate-data-in-time (runtime &optional (start-sweep 1)
					(signifier nil))
  "Generates a single datafile and a single progress file after every
  SAVE-EVERY-N-SWEEPS. Stops after runtime seconds. An hour is 3600
   seconds."
  (let ((filename (generate-filename-v3 signifier))
	(end-sweep (+ start-sweep NUM-SWEEPS -1))
	(endtime (+ (get-universal-time) runtime)) ; time measured in seconds
	(current-sweep 0))
    (while (< (get-universal-time) endtime)
      (sweep)
      (incf current-sweep)
      (when (= 0 (mod current-sweep SAVE-EVERY-N-SWEEPS))
	(make-spacetime-file filename)
	(make-progress-file filename start-sweep current-sweep end-sweep)))
    (make-spacetime-file filename)
    (make-progress-file filename start-sweep current-sweep end-sweep)))

;; generate-data-v2 is similar to generate-data except it creates a
;; fresh data file every SAVE-EVERY-N-SWEEPS. since a fresh datafile
;; is created, there is no need to maintain a seprate progress file.
(defun generate-data-v2 (&optional (start-sweep 1) (signifier nil))
  "Generate a spacetime file every SAVE-EVERY-N-SWEEPS. No progress file is
   generated."
  (setf SIM-START-TIME (cdt-now-str))
  (let ((end-sweep (+ start-sweep NUM-SWEEPS -1))
	(initial-boundary (count-spacelike-triangles-at-time 0))
	(final-boundary (count-spacelike-triangles-at-time NUM-T)))
    (when (= 1 start-sweep) ;; save the initial spacetime contents if
			    ;; this is a brand new run
      (make-spacetime-file (generate-filename-v2 initial-boundary
						 final-boundary
						 start-sweep 0
						 end-sweep
						 signifier))
    ; Start the sweeps
    (for (ns start-sweep end-sweep)
	 (sweep)
	 (when (= 0 (mod ns SAVE-EVERY-N-SWEEPS))
	   (make-spacetime-file (generate-filename-v2 initial-boundary
						      final-boundary
						      start-sweep ns
						      end-sweep
						      signifier)))))))

;; generate-data-in-time-v2 is similar to generate-data-in-time except
;; it creates a fresh data file every SAVE-EVERY-N-SWEEPS. since a
;; fresh datafile is created, there is no need to maintain a seprate
;; progress file.
(defun generate-data-in-time-v2 (runtime &optional (start-sweep 1)
					   (signifier nil))
  "Generate a spacetime file every SAVE-EVERY-N-SWEEPS. No progress file is
   generated."
  (setf SIM-START-TIME (cdt-now-str))
  (let ((end-sweep (+ start-sweep NUM-SWEEPS -1))
	(initial-boundary (count-spacelike-triangles-at-time 0))
	(final-boundary (count-spacelike-triangles-at-time NUM-T))
	(endtime (+ (get-universal-time) runtime))
	(current-sweep 0))
    
    ;; save the initial spacetime contents if this is a brand new run
    (make-spacetime-file (generate-filename-v2 initial-boundary
					       final-boundary
					       start-sweep
					       current-sweep
					       end-sweep
					       signifier))

    ; Start the sweeps
    (while (< (get-universal-time) endtime)
      (sweep)
      (incf current-sweep)
      (when (= 0 (mod current-sweep SAVE-EVERY-N-SWEEPS))
	(make-spacetime-file (generate-filename-v2 initial-boundary
						   final-boundary
						   start-sweep
						   current-sweep
						   end-sweep
						   signifier))))

    ; When time runs out, save immediately once more.
    (make-spacetime-file (generate-filename-v2 initial-boundary
					       final-boundary
					       start-sweep
					       current-sweep
					       end-sweep
					       signifier))))
    

;; generate-movie-data saves number of simplices every SAVE-EVERY-N-SWEEPS
(defun generate-movie-data (&optional (start-sweep 1) (signifier nil))
  "Makes a movie file and a progress file. No spacetime data is generated."
  (setf SAVE-EVERY-N-SWEEPS 10)
  (let* ((end-sweep (+ start-sweep NUM-SWEEPS -1))
	 (filename (generate-filename start-sweep end-sweep signifier)))
	

    ;; Initialize the movie file
    (make-movie-file filename start-sweep)
    
    ;; Run the sweeps
    (for (ns start-sweep end-sweep)
	 (sweep)
	 (when (= 0 (mod ns SAVE-EVERY-N-SWEEPS))
	   (append-to-movie-file filename)
	   (make-progress-file filename start-sweep ns end-sweep)))))

(defun generate-movie-data-console (&optional (start-sweep 1))
  "Generates movie data and console data."
  (when (= 1 start-sweep)
    (reset-move-counts))
  (let ((end-sweep (+ start-sweep NUM-SWEEPS -1))
	(filename (generate-filename start-sweep)))
    (make-movie-file filename start-sweep) ; initialize spacetime file
    (for (ns start-sweep end-sweep) ;sweep
      (sweep)
      (when (= 0 (mod ns 10))
	(format t "~A/~A " ns end-sweep)
	(for (ts 0 (1- NUM-T))
	  (format t "~A " (count-simplices-in-sandwich ts (1+ ts))))
	(format t "~A ~A~%" (count-simplices-of-all-types) 
		(accept-ratios))
	(make-progress-file filename start-sweep ns end-sweep)
	(append-to-movie-file filename)))))

(defun generate-spacetime-and-movie-data (&optional (start-sweep 1)
					    (signifier nil))
  "Generate 3-simplex data and movie data. Keeps a progress file."
  (let* ((end-sweep (+ start-sweep NUM-SWEEPS -1))
	 (filename (generate-filename start-sweep end-sweep signifier)))
    
    ;; open and close the file, for :append to work properly
    ;; record the initial data only if start-sweep = 1
    (make-movie-file filename start-sweep)
    
    ;; Sweep
    (for (ns start-sweep end-sweep)
	 (sweep)
	 (when (= 0 (mod ns SAVE-EVERY-N-SWEEPS))
	   (make-spacetime-file filename)
	   (make-progress-file filename start-sweep ns end-sweep)
	   (append-to-movie-file filename)))))
;;;;-------------------------------------------------------------------------

