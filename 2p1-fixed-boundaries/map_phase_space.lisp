;;;; map_phase_space.lisp
;;;; Jonah Miller (jonah.miller@colorado.edu)
;;;; Date: July 3, 2012

;;;; This is a script that should vary the coupling constants for the
;;;; 2p1-fixed-boundaries code to map out useful points in the phase
;;;; space. This is an attempt to make a program that works faster
;;;; than the equivalent python program, by manipulating the
;;;; simulation environment directly.

;; Load module
;;;--------------------------------------------------------------------------
(load "../utilities.lisp")
;;;--------------------------------------------------------------------------

;;; Utility functions
;;;--------------------------------------------------------------------------
(defun format-hour-day-month-year ()
  "Gets the hour, day, month, and year, and formats it nicely."
  (let* ((datelist (multiple-value-list (get-decoded-time)))
	 (hours (third  datelist))
	 (days  (fourth datelist))
	 (month (fifth  datelist))
	 (year  (sixth  datelist)))
    (format nil "~ah-~ad-~am-~a" hours days month year)))
;;;--------------------------------------------------------------------------

;;; Global constants
;;;--------------------------------------------------------------------------
;; k0 data
(defvar *k0-min* 0.5 "The value we start varying k0 over.")
(defvar *k0-max* 3.5 "The value we end our variation over k0 with.")
(defvar *k0-increment* 0.1 "How much we increase k0 as we loop over it.")

;; k3 data
(defvar *k3-min* 0.6 "The minimum value of k3 for each k0.")
(defvar *k3-max* 3.5 "The maximum value of k3 for each k0.")
(defvar *k3-increment* 0.001 "How much we increase k3 as we loop over it.")

;; Simulation parameters
(defvar *target-volume* (* 8 1024)
  "The volume we want our simulation to thermalize to.") 
(defvar *acceptable-range* (map 'list 
				#'(lambda (x) (* *target-volume* x)) 
				'(0.1 1.1))
  "The minimum and maximum values of volume we find acceptable as we vary.")
(defvar *num-time-slices* 64 "Number of proper time slices")
(defvar *spatial-topology* "S2" "spatial topology of simulation")
(defvar *boundary-condition-type* "OPEN" "Boundary condition type")
(defvar *n-sweeps* 50000 "Just to make the simulation work.")
(defvar *thermalization-sweeps* 10 
  "The number of sweeps to reach target volume, plus a little bit.")
(defvar *data-sweeps* 5 "Sweeps we do statistics on.")


;; Default initial and final geometries
(defparameter initial-geometry "tetra.txt" "Default initial geometry.")
(defparameter final-geometry   "tetra.txt" "Default final geometry.")

;; Various
(defvar *outfile* (format nil "~a-~a.phasemap" 
			  (format-hour-day-month-year) (hostname))
  "The file name for the output file")
;;;--------------------------------------------------------------------------


;;; Program-Specific Functions
;;;--------------------------------------------------------------------------
(defun start-simulation (&optional (initial-geometry "tetra.txt")
			   (final-geometry "tetra.txt"))
  "Starts the simulation running."
  (load "cdt2p1.lisp")    ;; load simulation
  (setf NUM-SWEEPS *n-sweeps*) ;; Just so the simulation works
  (initialize-t-slices-with-v-volume :num-time-slices     *num-time-slices*
				    :target-volume       *target-volume*
				    :spatial-topology    "S2"
				    :boundary-conditions "OPEN"
				    :initial-spatial-geometry initial-geometry
				    :final-spatial-geometry final-geometry))

(defun open-file-for-the-first-time (outfile)
  "Opens the output file for the fist time and 
puts in the proper information."
  (with-open-file (f outfile 
		     :direction :output
		     :if-exists :supersede
		     :if-does-not-exist :create)
    (format f "# k0~Tk3~TN3-mean~Tsigma~Tdelta~%")))

(defun add-data-to-file (outfile k0 k3 n3-mean n3-sigma n3-delta)
  "Adds the requisite information to the output file."
  (with-open-file (f outfile
		     :direction :output
		     :if-exists :append
		     :if-does-not-exist :create)
    (format f "~3D~T~3D~T~3D~T~3D~T~3D~%" k0 k3 n3-mean n3-sigma n3-delta)))

(defun collect-volume-data (thermalization-sweeps data-sweeps 
				&optional (start-sweep 1))
  "Similar to the collect data functions in montecarlo.lisp, 
but designed to test the phase space."
  (let ((3-volumes nil))
    (prog nil
       (for (ns start-sweep (+ start-sweep thermalization-sweeps -1))
	 (sweep))
       (for (ns start-sweep (+ start-sweep data-sweeps -1))
	 (sweep)
	 (push (+ N3-TL-31 N3-TL-22) 3-volumes)))
    3-volumes))

(defun mean (list-of-values)
  "Calculate the mean of a list of values"
  (/ (reduce #'+ list-of-values) (length list-of-values)))

(defun standard-deviation (list-of-values)
  "Calculate the standard deviation of a list of values"
  (flet ((sqr (x) (* x x)))
    (let* ((n (length list-of-values))
	   (mean (mean list-of-values))
	   (deviation-list (mapcar #'(lambda (x) (sqr (- x mean))) 
				   list-of-values)))
      (sqrt (/ (reduce #'+ deviation-list) (1- n))))))

(defun get-delta (list-of-values)
  "Calculate a crude derivative for a
 list of values (as a function of list index."
  (- (nth (1- (length list-of-values)) list-of-values) 
     (first list-of-values)))

(defun vary-k0-k3 (outfile thermalization-sweeps data-sweeps
		   k0-min k0-max k0-increment
		   k3-min k3-max k3-increment
		   &optional (initial-geometry "tetra.txt")
		     (final-geometry "tetra.txt"))
  "Varies k0 and k3 to map out the phase space."
  (prog nil
     (start-simulation initial-geometry final-geometry)
     (open-file-for-the-first-time outfile)
     (loop for k0 from k0-min to k0-max by k0-increment do
	  (loop for k3 from k3-min to k3-max by k3-increment do
	       (prog ((data nil))
		  (set-k0-k3-alpha k0 k3 -1)
		  (setf data (collect-volume-data thermalization-sweeps
						   data-sweeps))
		  (add-data-to-file outfile k0 k3 
				    (float (mean data))
				    (float (standard-deviation data))
				    (float (get-delta data))))))))
;;;--------------------------------------------------------------------------


;;; The main program
;;;--------------------------------------------------------------------------
(vary-k0-k3 *outfile* *thermalization-sweeps* *data-sweeps* 
	    *k0-min* *k0-max* *k0-increment*
	    *k3-min* *k3-max* *k3-increment*)
;;;--------------------------------------------------------------------------
