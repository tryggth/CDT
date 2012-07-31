;;;; ascii_plotting_tools.lisp
;;;; Author: Jonah Miller
;;;; Date: July 30th, 2012

;;;; The goal of this module for 2p1-fixed-boundaries is to allow for
;;;; very quick and dirty plots of a spacetime during repl
;;;; sessions. However, this basic code might be useful other places.
;;;; The primary function, ascii-list-plot is in ../utilities.lisp


;;;; Additional counting functions required for plotting functions.
;;;;-----------------------------------------------------------------------
(defun count-3-simplices-of-time nil
  "Counts the number of 3 simplices connected to proper time tm-low
   as a function of tm-low. tm-low is the lower time slice a simplex is
   attached to. Returns a list. Index is tm-low."
  (let ((simplex-count nil)
	(t-max (if (string= BCTYPE "OPEN") (1- NUM-T) (- NUM-T 2))))
    (loop for i from 0 to t-max do
	 (push (count-simplices-in-sandwich i (1+ i)) simplex-count))
    (reverse simplex-count)))

(defun count-spacelike-triangles-of-time nil
  "Counts the number of spacelike triangles in time-slice tm
   as a function of tm. Returns a list. Index is tm."
  (let ((simplex-count nil)
	(t-max (if (string= BCTYPE "OPEN") NUM-T (1- NUM-T))))
    (loop for i from 0 to t-max do
	 (push (count-spacelike-triangles-at-time i) simplex-count))
    (reverse simplex-count)))
;;;;-----------------------------------------------------------------------


;;;; Plotting functions
;;;;-----------------------------------------------------------------------
(defun plot-3-simplices-of-time nil
  "Plots the number of 3 simplices connected to proper time tm-low
   as a function of tm-low. tm-low is the lower time slice a simplex is
   attached to."
  (ascii-list-plot (count-3-simplices-of-time) "3-simplex count:"))

(defun plot-spacelike-triangles-of-time nil
  "Plots the number of spacelike triangles in time-slice tm
   as a function of tm."
  (ascii-list-plot (count-spacelike-triangles-of-time) 
		   "Spacelike Triangle Count:"))
;;;;-----------------------------------------------------------------------

