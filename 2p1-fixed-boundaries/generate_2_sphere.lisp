;;;; generate_2_sphere.lisp
;;;; Jonah Miller (jonah.miller@colorado.edu)
;;;; Date: July 5, 2012

;;;; This is a lisp program that should generate an approximation of a
;;;; sphere using equilateral triangles. It then generates a list of
;;;; triples to give as input to 2p1-fixed-boundaries.

;; Load Modules
;;;--------------------------------------------------------------------------
(load "../utilities.lisp")
;;;--------------------------------------------------------------------------

;;; Global constants
;;;--------------------------------------------------------------------------
;; initial triangulation
(defvar *tetrahedron* '((1 2 3) (1 2 4) (1 3 4) (2 3 4)) 
  "The tetrahedron we start with to generate our triangulation (if we choose)."
(defvar *octahedron* '((1 2 3) (2 3 5) (2 4 5) (4 5 6)
		       (3 5 6) (1 3 6) (1 4 6) (1 2 4))
  "The octahedron we start with to generate our triangulation (if we choose).")
(defvar *icosohedron* '((1 2 3)  (2 3 4)  (3 4 5)  (4 5 6)   (5 6 7)
			(6 7 8)  (7 8 9)  (8 9 10) (9 10 1)  (1 2 10)
			(1 3 11) (3 5 11) (5 7 11) (7 9 11)  (9 1 11)
			(2 4 12) (4 6 12) (6 8 12) (8 10 12) (10 2 12))
  "The icosahedron we start with to generate our triangulation (if we choose).")

;;; Data structures (hash tables
(defparameter *triangles* ()
  "The list containing triangles. Each triangle is a tripple of point numbers.
 Position in the list matters!")
(defparameter *edges* ()
  "The list containing edges of triangles. Each edge is a double of point
 numbers. Position in the list matters!")
(defparameter *points* ()
  "The list containing points contained in the triangle. In the case of
 points, order in the list doesn't matter.")
;;;--------------------------------------------------------------------------


;;; Functions
;;;--------------------------------------------------------------------------
(defun set-equal (set1 set2)
  "Checks to see if two sets are equal."
  (and (subsetp set1 set2) (subsetp set2 set1)))

(defun gen-member (object list)
  "Checks to see if a set is in a list of sets."
  (cond 
    ((null list)                     nil)
    ((set-equal object (first list)) object)
    (t                               (gen-member object (rest list)))))


(defun generate-edges (triangle)
  "Generate the edge doubles in a triangle."
  (let ((edgelist nil))
    (do-tuples/c (x y) triangle (push (list x y) edgelist))
    edgelist))

(defun initialize-triangulation (initial-tetrahedron)
  "Initialize the first tetrahedron we use triangulate the sphere."
  ;; Loop over triangles
  (loop for triangle in initial-tetrahedron do
       ;; Add points in triangle to list (if they're not already there)
       (loop for point in triangle do
	    (when (not (member point *points*))
	      (push point *points*)))
       ;; Add edges in triangle to list (if they're not already there)
       (let ((edgelist (generate-edges triangle)))
	 (loop for edge in edgelist do
	      (when (not (member edge *edges*))
		(push edge *edges*))))
       ;; Add the triangle to the triangle list
       (when 
  ;; Add links to the link list 
  (loop for triangle in initial-tetrahedron do

;;;--------------------------------------------------------------------------
