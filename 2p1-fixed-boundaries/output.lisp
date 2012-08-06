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
(defun generate-data (&optional (start-sweep 1))
  "Generates a single datafile and a single progress file after every
   SAVE-EVERY-N-SWEEPS."
  (let ((filename (generate-filename start-sweep))
	(end-sweep (+ start-sweep NUM-SWEEPS -1)))
    (for (ns start-sweep end-sweep)
	 (sweep)
	 (when (= 0 (mod ns SAVE-EVERY-N-SWEEPS))
	   (make-spacetime-file filename)
	   (make-progress-file filename start-sweep ns end-sweep)))))

;; generate-data-v2 is similar to generate-data except it creates a
;; fresh data file every SAVE-EVERY-N-SWEEPS. since a fresh datafile
;; is created, there is no need to maintain a seprate progress file.
(defun generate-data-v2 (&optional (start-sweep 1))
  "Generate a spacetime file every SAVE-EVERY-N-SWEEPS. No progress file is
   generated."
  (setf SIM-START-TIME (cdt-now-str))
  (let ((end-sweep (+ start-sweep NUM-SWEEPS -1)))
    (when (= 1 start-sweep) ;; save the initial spacetime contents if
			    ;; this is a brand new run
      (make-spacetime-file (generate-filename-v2 start-sweep 0)))
    ; Start the sweeps
    (for (ns start-sweep end-sweep)
	 (sweep)
	 (when (= 0 (mod ns SAVE-EVERY-N-SWEEPS))
	   (make-spacetime-file (generate-filename-v2 start-sweep ns))))))


;; generate-data-v3 is similar to generate-data-v2 except it also
;; creates an additional data file every SAVE-EVERY-N-SWEEPS that
;; contains the spatial 2-simplex information for each spatial slice.
(defun generate-data-v3 (&optional (start-sweep 1))
  "Like generate-data-v2, but generates spatial 2 simplex and spatial 
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

;; generate-movie-data saves number of simplices every SAVE-EVERY-N-SWEEPS
(defun generate-movie-data (&optional (start-sweep 1))
  "Makes a movie file and a progress file. No spacetime data is generated."
  (setf SAVE-EVERY-N-SWEEPS 10)
  (let ((filename (generate-filename start-sweep))
	(end-sweep (+ start-sweep NUM-SWEEPS -1)))

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

(defun generate-spacetime-and-movie-data (&optional (start-sweep 1))
  "Generate 3-simplex data and movie data. Keeps a progress file."
  (let* ((end-sweep (+ start-sweep NUM-SWEEPS -1))
	 (filename (generate-filename start-sweep)))
    
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

