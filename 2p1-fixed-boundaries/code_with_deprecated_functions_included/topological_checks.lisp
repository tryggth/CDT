;; topological_checks.lisp
;; Author: Jonah Miller (jonah.maxwell.miller@gmail.com)
;; Date: July 27, 2012

;; This program is a module meant to be used with the 2+1-dimensional
;; CDT simulations to check that topological identities are not
;; violated. A number of topological identities are defined as
;; functions, and checks are defined to determine that the identity
;; returns the correct thing (for instance to ensure that the Euler
;; characteristic of a sphere is 2).

;; Requires the following modules to be loaded first:
;; ../utilities.lisp
;; globals.lisp
;; generalized-hash-table-counting-functions.lisp
;; simplex.lisp



;;;; Constants
;;;;-------------------------------------------------------------------------
(defvar *chi-sphere* 2 "The Euler characteristic of a sphere.")
(defvar *chi-torus*  0 "The Euler characteristic of a torus.")
;;;;-------------------------------------------------------------------------


;;;; Topological identities
;;;;-------------------------------------------------------------------------
(defun euler-char (time-slice)
  "Returns the Euler characteristic of a given time slice t0. 
   The formula for the Euler characteristic is:
       Vertices - Edges + Faces
   of our triangulation."
  (let* ((t0 (bc-mod time-slice))
	 (vertices (count-points-at-time t0))
	 (edges    (count-spacelike-links-at-time t0))
	 (faces    (count-spacelike-triangles-at-time t0))
	 (chi (- (+ vertices faces) edges)))
    chi))

(defun triangle-edge-relation (time-slice)
  "A topological relation relating the number of edges on a given 
   time slice to the number of triangles. The relation is
   edges - (3/2) faces
   If topological identities hold, this should return zero."
  (let* ((t0 (bc-mod time-slice))
	 (edges (count-spacelike-links-at-time t0))
	 (faces (count-spacelike-triangles-at-time t0))
	 (relationship (- edges (* 3/2 faces))))
    relationship))

(defun tetrahedron-edge-relation (time-slice)
  "A topological relation relating the number of tetrahedra connecting
   to a given time slice and the number of edges in that
   slice. The relation is:
   (N3-TL-31 + N3-TL-13)(connected to slice) - (4/3) N1-SL-TL(time-slice).
   If the identities hold, this should return zero.

   Note that this identity will not necessarily hold on the initial or
   final slices in the open boundary conditions case."
  (let* ((t0 (bc-mod time-slice))
	 (13s (count-simplices-in-sandwich-of-type (1- t0) t0 1))
	 (31s (count-simplices-in-sandwich-of-type t0 (1+ t0) 3))
	 (sum-n3-31 (+ 13s 31s))
	 (edges (count-spacelike-links-at-time t0))
	 (relationship (- sum-n3-31 (* 4/3 edges))))
    relationship))

(defun timelike-face-tetrahedron-relation (time-slice)
  "A topological relation relating the number of spacelike 2-simplices 
   in a given sandwich to the number of 3-simplices in that sandwich.
   The relation comes from the fact that each timelike face must be part of 
   two 3-simplices. The relation is:
   3 N3-TL-31 + 3 N3-TL-13 + 4 N3-TL-22 = 2 N2-TL
   If the identities hold, this will return zero."
  (let* ((t0 (bc-mod time-slice))
	 (13s (count-simplices-in-sandwich-of-type t0 (1+ t0) 1))
	 (31s (count-simplices-in-sandwich-of-type t0 (1+ t0) 3))
	 (22s (count-simplices-in-sandwich-of-type t0 (1+ t0) 2))
	 (timelike-faces (count-timelike-triangles-in-sandwich t0 (1+ t0)))
	 (relationship 
	  (+ (* 3 31s) (* 3 13s) (* 4 22s) (* -2 timelike-faces))))
    relationship))
;;;;-------------------------------------------------------------------------


;;;; Topological checks
;;;;-------------------------------------------------------------------------
(defun check-topological-identity (identity expected-result)
  "The template to check a topological identity. Identity is a symbol
   function for one of the topological identities above. expected-result is
   the expected result of the relation."
  (let ((t-max (if (string= BCTYPE "OPEN") NUM-T (1- NUM-T)))
	(failed-time-slices nil)
	(failed-relations nil))
    (for (t0 0 t-max)
      (let ((result (funcall identity t0)))
	(when (/= result expected-result)
	  (push t0 failed-time-slices)
	  (push result failed-relations))))
    (list failed-time-slices failed-relations)))

(defun check-euler-characteristic ()
  "Checks the Euler characteristic of every time slice in the spacetime.
   If the Euler characteristic is not correct, prints the proper time of 
   the slice and the actual Euler characteristic."
  (let ((expected-result (if (string= STOPOLOGY "S2") 
			     *chi-sphere*
			     *chi-torus*)))
    (check-topological-identity #'euler-char expected-result)))

(defun check-triangle-edge-relation ()
  "Checks the triangle-edge-relation given above. If we have a topologically
   acceptable spacetime, then the result should be zero. If it is not, prints
   the proper time of the slice and the actual result of the relation."
  (check-topological-identity #'triangle-edge-relation 0))

(defun check-tetrahedron-edge-relation ()
  "Checks the tetrahedron-edge-relation given above. If we have a
   topologically acceptable spacetime, then the result should be
   zero (except at the boundary). If it is not, prints the proper time of
   the slice and the actual result of the relation."
  (check-topological-identity #'tetrahedron-edge-relation 0))

(defun check-timelike-face-tetrahedron-relation ()
  "Checks the timelike-face-tetrahedron-relation given above. If we have a
   topologically acceptable spacetime, then the result should be zero and the
   function returns nil. If it is not, prints the proper time of the lower 
   time slice in the sandwich and the actual result of the relation."
  (check-topological-identity #'timelike-face-tetrahedron-relation 0))
;;;;-------------------------------------------------------------------------
