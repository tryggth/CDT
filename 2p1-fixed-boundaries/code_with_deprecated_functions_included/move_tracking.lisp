;; move_tracking.lisp
;; Author: Jonah Miller (jonah.maxwell.miller@gmail.com)
;; Date: August 3, 2012

;; This program is a module meant to be used with the 2+1-dimensional
;; CDT simulations to track how many moves of a given type may be
;; attempted at each time slice.

;; Requires the following modules to be loaded first:

;; ../utilities.lisp
;; globals.lisp
;; generalized-hash-table-counting-functions.lisp
;; simplex.lisp
;; moves.lisp

;; If you want to plot what you get, you also need
;; ascii_plotting_tools.lisp

;; If you want to make movies, you also need montecarlo.lisp (the
;; program will raise errors at compile time if you do not include
;; montecarlo.lisp, but this won't cause any actual errors).

;;;; Constants
;;;;-------------------------------------------------------------------------
;; Move count file name endings
(defvar *move-ext* ".count2p1" 
  "Base file extension used for looking at the number of simplices that are 
   topologically acceptable for a move to act on them.")
(defvar *26-move-ext* (concatenate 'string ".26move" *move-ext*)
  "Used for storing percentage of simplices per time slice that are
   topologically acceptable for a 2->6 move to operate on them.")
(defvar *23-move-ext* (concatenate 'string ".23move" *move-ext*)
  "Used for storing percentage of simplices per time slice that are
   topologically acceptable for a 2->3 move to operate on them.")
(defvar *44-move-ext* (concatenate 'string ".44move" *move-ext*)
  "Used for storing percentage of simplices per time slice that are
   topologically acceptable for a 4->4 move to operate on them.")
(defvar *32-move-ext* (concatenate 'string ".32move" *move-ext*)
  "Used for storing percentage of simplices per time slice that are
   topologically acceptable for a 3->2 move to operate on them.")
(defvar *62-move-ext* (concatenate 'string ".62move" *move-ext*)
  "Used for storing percentage of simplices per time slice that are
   topologically acceptable for a 6->2 move to operate on them.")
(defvar *default-move-ext* (concatenate 'string ".move" *move-ext*)
  "Just a default file extension.")

(defun move-count-extension (movefunction)
  "While this is not strictly a constant, it maps a move function
   (such as try-2->3) to a file name extension. If you add more 
   extensions or moves, change this too. It's at the top of the file 
   for that reason."
  (cond
    ((equalp movefunction #'try-2->6) *26-move-ext*)
    ((equalp movefunction #'try-2->3) *23-move-ext*)
    ((equalp movefunction #'try-4->4) *44-move-ext*)
    ((equalp movefunction #'try-3->2) *32-move-ext*)
    ((equalp movefunction #'try-6->2) *62-move-ext*)
    (t *default-move-ext*)))

(defun move-type-names (movefunction)
  "Returns the type of move as a string. A primitive type cast. Useful 
   for debugging tools mostly."
  (cond
    ((equalp movefunction #'try-2->6) "2->6")
    ((equalp movefunction #'try-2->3) "2->3")
    ((equalp movefunction #'try-4->4) "4->4")
    ((equalp movefunction #'try-3->2) "3->2")
    ((equalp movefunction #'try-6->2) "6->2")
    (t *default-move-ext*)))
   

(defvar *move-functions* (list #'try-2->6 #'try-2->3 #'try-4->4 
			       #'try-3->2 #'try-6->2)
  "A list of the different try move functions. 
   Used for certain loops below.")
;;;;-------------------------------------------------------------------------


;;;; Data Taking Functions
;;;;-------------------------------------------------------------------------

;;;; These functions take data on what moves are topologically
;;;; acceptable. They can be fed into output functions or plotting
;;;; functions below.

(defun map-move-in-sandwich (move-function t0)
  "Find out how many simplices with t-low t0 (i.e., in the t0,t0+1 sandwich)
   are topologically acceptable for a given move."
  (let* ((attempts 0)
	 (simplices (get-simplices-in-sandwich t0 (1+ t0)))
	 (num-simplices 0)) ; not necessary, but faster than (length)
    (loop for simplex in simplices do
	 (incf num-simplices)
	 (when (funcall move-function simplex)
	   (incf attempts)))
    (values attempts num-simplices)))

(defun percent-moves-in-sandwich (move-function t0)
  "Find out what percentage of simplices in the t0,t0+1 sandwich are
   topologically acceptable targets for a given move."
  (* 100 (float (/ (map-move-in-sandwich move-function t0)
		   (count-simplices-in-sandwich t0 (1+ t0))))))

(defun percent-each-move-in-sandwich (t0)
  "Prints the percent of a move in the sandwich as above, but does so for 
   every move possible."
  (loop for move in *move-functions* do
       (format t "~A: ~A~%" 
	       (move-type-names move)
	       (percent-moves-in-sandwich move t0))))
	       
(defun map-move-volume-across-time (move-function)
  "For each time slice sandwich, checks to see on how many 3-simplices
   a move is topologically acceptable. Returns a pair of lists. The
   first is the total number of topologically acceptable simplices for
   each time slice. The second is the total number of 3-simplices at
   that time slice. The index is time."
  (let* ((attempt-count nil) ; The number of simplices which are
			     ; topologically acceptable for a move
	 (simplex-count nil) ; The total number of simplices at a
			     ; given time slice.
	 ; The number of time slices we iterate over. 
	 ; Open boundaries have one extra slice.
	 (t-max (if (string= BCTYPE "OPEN") (1- NUM-T) (- NUM-T 2))))
    (loop for i from 0 to t-max do
	 (multiple-value-bind (attempts simplices)
	     (map-move-in-sandwich move-function i)
	   (push attempts attempt-count)
	   (push simplices simplex-count)))
    (values (reverse attempt-count) (reverse simplex-count))))

(defun percent-moves-across-time (move-function)
  "For each time slice sandwich, checks to see on how many 3-simplices
   a move is topologically acceptable. Returns a list where each
   element is the percentage of simplices that are topologically
   acceptable at a given time. The index is time."
  (multiple-value-bind (attempts simplices)
      (map-move-volume-across-time move-function)
    (make-percent attempts simplices)))
  
;;;;-------------------------------------------------------------------------



;;;; Output functions
;;;;-------------------------------------------------------------------------

;;;; These functions take the data taking functions above and produce
;;;; output files.

(defun print-move-percent (iostream move-function)
  "This function will run percent-moves-across-time with move-function
   on the spacetime and append the results to the end of the
   iostream."
  (format-list iostream (percent-moves-across-time move-function)))

(defun make-new-move-percent-file (move-function &optional (start-sweep 1))
  "Makes a new move-percent file and prints the current move percents to it.
   start-sweep is the starting sweep for monte carlo sweeps. For this 
   function to work properly NUM-SWEEPS must be set.

   Returns the filename."
  (let ((filename (concatenate 'string 
			       (generate-filename start-sweep) 
			       (move-count-extension move-function))))
    (with-open-file (countfile filename 
			       :direction :output
			       :if-exists :supersede)
      (print-move-percent countfile move-function)
      (format countfile "~%"))
    filename))

(defun append-move-percent (filename move-function)
  "Appends the current move percents to the file with filename. Usually gets
   the filename from make-new-move-percent-file."
  (with-open-file (countfile filename
			     :direction :output
			     :if-exists :append)
    (print-move-percent countfile move-function)
    (format countfile "~%")))

(defun make-all-new-move-percent-files (&optional (start-sweep 1))
  "Like make-new-move-percent-file above, but does so for all moves in 
   *move-functions*. Returns a list of filenames."
  (let ((filenames nil))
    (loop for move in *move-functions* do
	 (push (make-new-move-percent-file move start-sweep) filenames))
    (reverse filenames)))

(defun append-to-all-move-percent-files (filename-list)
  "Like append-move-percent above, but does so for an entire list of
   moves and for all moves in *move-functions*. Relies on the ordering
   of the filename-list, so it will probably only work properly if it
   takes the output of make-all-new-move-percent-files as input."
  (for (i 0 (1- (length *move-functions*)))
    (append-move-percent (nth i filename-list) (nth i *move-functions*))))
;;;;-------------------------------------------------------------------------



;;;; Monte Carlo Output Functions
;;;;-------------------------------------------------------------------------

;;;; These functions work like output functions, but they work over an
;;;; entire monte-carlo data-taking/thermalization simulation.

(defun generate-spacetime-movie-and-count-data (&optional (start-sweep 1))
  "identical to generate-spacetime-and-movie-data, except also
   generates movie data containing move information."
  (let* ((end-sweep (+ start-sweep NUM-SWEEPS -1))
	 (datafilename (generate-filename start-sweep))
	 (move-file-names (make-all-new-move-percent-files start-sweep)))
    ;; open and close the file, for :append to work properly
    ;; record the initial data only if start-sweep = 1
    (make-movie-file datafilename start-sweep)

    ;; Sweep
    (for (ns start-sweep end-sweep)
      (sweep)
      (when (= 0 (mod ns SAVE-EVERY-N-SWEEPS))
	(make-spacetime-file datafilename)
	(make-progress-file datafilename start-sweep ns end-sweep)
	(append-to-movie-file datafilename)
	(append-to-all-move-percent-files move-file-names)))))
;;;;-------------------------------------------------------------------------
