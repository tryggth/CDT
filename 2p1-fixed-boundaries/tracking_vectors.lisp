;;;; tracking_vectors.lip
;;;; Authors: Jonah Miller (jonah.maxwell.miller@gmail.com)
;;;;          Rajesh Kommu

;;;; This file contains the f- and b-vectors and the functions that
;;;; manipulate them.


;;; Set this to true if you want to test the program with the
;;; boundaries identified (i.e., periodic boundary conditions) but you
;;; still want to keep track of the boundary. In this case, you also
;;; need to set the boundary condition type to "PERIODIC.
(defparameter *merge-faces* nil)


;;;; These are elements of the f-vector. The f-vector keeps track of
;;;; the number of simplices of each type in the entire simplectic
;;;; manifold. It is faster to keep track as we run the simulation
;;;; than to count up the number after every move. The b-vector, which
;;;; is below, counts the number of simplices of each type that are only on
;;;; the boundary.
(defparameter N0 0 "number of points")
(defparameter N1-SL 0 "number of spacelike links")
(defparameter N1-TL 0 "number of timelike links")
(defparameter N2-SL 0 "number of spacelike triangles")
(defparameter N2-TL 0 "number of timelike triangles")
(defparameter N3-TL-31 0 "number of (1,3) + (3,1) timelike tetrahedra")
(defparameter N3-TL-22 0 "number of (2,2) timelike tetrahedra")
(defmacro N3 ()
  "total number of timelike 3simplices (tetrahedra)"
  `(+ N3-TL-31 N3-TL-22))

;;; Book-keeping for the action terms. We input the f-vector and
;;; b-vector into the action to determine whether or not to accept a
;;; move. It is worth noting that we ALSO only accept a move depending
;;; on whether or not it changes the boundary. Changes to the boundary
;;; are not allowed.

;; Some macros to determine whether or not a simplex is on a boundary
;;first, define macros to determine helpful things about
;;the position of a particular simplex
(defmacro in-upper-sandwich (sxid) 
  "Is a 3-simplex in the sandwich (t_max-1,t_max)?"
  `(and 
    (or (string= BCTYPE "OPEN") *merge-faces*)
    (= (3sx-tmhi (get-3simplex ,sxid)) (bc-mod NUM-T))))
(defmacro in-lower-sandwich (sxid)
  "Is a 3-simplex in the sandwich (0,1)?"
  `(and 
    (or (string= BCTYPE "OPEN") *merge-faces*)
    (= (3sx-tmlo (get-3simplex ,sxid)) 0)))
(defmacro in-either-boundary-sandwich (sxid)
  `(or (in-upper-sandwich ,sxid) (in-lower-sandwich ,sxid)))

;; Necessary if we want to unfix the boundaries, but still keep track
;; of them. 

;; JM: to have non-periodic unfixed boundaries, you will need new
;; moves that act on the boundary. These are not implimented.
(defmacro has-face-on-boundary (sxid)
  "Is one of the triangular faces of a 3-simplex contained in a boundary?"
  `(let* ((sx (get-3simplex ,sxid))  
	  (ty (3sx-type sx))
	  (th (bc-mod (3sx-tmhi sx)))
	  (tl (bc-mod (3sx-tmlo sx))))
     (and (or (string= BCTYPE "OPEN") 
	      *merge-faces*)
	  (or (and (= ty 1)
		   (or (= th 0) (= th NUM-T)))
	      (and (= ty 3)
		   (or (= tl 0) (= tl NUM-T)))))))

;; The f-vector manipulation functions
(defun set-f-vector (v1 v2 v3 v4 v5 v6 v7)
  (setf N0 v1 N1-SL v2 N1-TL v3 N2-SL v4 N2-TL v5 N3-TL-31 v6 N3-TL-22 v7))

(defun update-f-vector (dv)
  (incf N0       (nth 0 dv))
  (incf N1-SL    (nth 1 dv))
  (incf N1-TL    (nth 2 dv))
  (incf N2-SL    (nth 3 dv))
  (incf N2-TL    (nth 4 dv))
  (incf N3-TL-31 (nth 5 dv))
  (incf N3-TL-22 (nth 6 dv)))

;; Changes in the f-vector depending on moves.
(defparameter DF26 '(1 3 2 2 6 4 0))
(defparameter DF62 '(-1 -3 -2 -2 -6 -4 0))
(defparameter DF44 '(0 0 0 0 0 0 0))
(defparameter DF23 '(0 0 1 0 2 0 1))
(defparameter DF32 '(0 0 -1 0 -2 0 -1))

(defun f-vector ()
  "Print the f-vector."
  (format nil "[~A ~A ~A ~A ~A ~A ~A] : ~A ~%"
	  N0
	  N1-SL
	  N1-TL
	  N2-SL
	  N2-TL
	  N3-TL-31
	  N3-TL-22
	  (N3)))


;;; In addition to the f-vector, which keeps track of total information
;;; (not bulk information), we also have the b-vector, which keeps
;;; track of boundary information only. This is useful for
;;; fixed-boundary systems. Below is all the information the system
;;; needs to keep track of the boundary (not including hash tables).

;;; The b-vector keeps track of the top boundary (at t=t_final) and
;;; the bottom boundary (at t=t_initial) seperately. This is because
;;; these terms are accounted for in the action with opposite sign.

;;; Note: the following definitions are not in the order they appear
;;; in the b-vector!

;; Define some quantities of geometrical objects on the boundaries
(defparameter *N1-SL-TOP* 0) ; Number of spacelike links on the
			     ; boundary. These are the bones used to
			     ; calculate the extrinsic curvature on
			     ; the boundary. For the top boundary.
(defparameter *N1-SL-BOT* 0) ; For the bottom boundary

(defparameter *N3-22-TOP* 0) ; Number of (2,2)-simplexes that have 2
			     ; vertexes on the boundary. Useful for
			     ; extrinsic curvature calculations. For
			     ; the top boundary.
(defparameter *N3-22-BOT* 0) ; For the bottom boundary.

(defparameter *N3-31-TOP* 0) ; Number of (3,1)- or (1,3)-simplexes
			     ; that have at least one vertex on the
			     ; boundary. We don't need them for
			     ; extrinsic curvature calculations, but
			     ; we need to compare them to N3-31
			     ; (total) to figure out the number of
			     ; (3,1)-simplices in the bulk. For the
			     ; top boundary.
(defparameter *N3-31-BOT* 0) ; For the bottom boundary.

;; The above three quantities compose the "b-vector," which keeps
;; track of boundary information. The following two functions
;; manipulate the b-vector.

;; Extract the order in the b-vector from these definitions.
(defun set-b-vector (n1 n2 n3 n4 n5 n6)
  "Input: *N1-SL-TOP*, *N1-SL-BOT*, *N3-22-TOP*, *N3-22-BOT*, 
*N3-31-TOP*, *N3-31-BOT*."
  (setf *N1-SL-TOP* n1
	*N3-22-TOP* n2
	*N3-31-TOP* n3
	*N1-SL-BOT* n4
	*N3-22-BOT* n5
	*N3-31-BOT* n6))

(defun update-b-vector (dv)
  "Input list: (DELTA-N1-SL-TOP DELTA-N1-SL-BOT
 DELTA-N3-22-TOP DELTA-N3-22-BOT DELTA-N3-31-TOP DELTA-N3-31-BOT)"
  (incf *N1-SL-TOP* (first  dv))
  (incf *N3-22-TOP* (second dv))
  (incf *N3-31-TOP* (third  dv))
  (incf *N1-SL-BOT* (fourth dv))
  (incf *N3-22-BOT* (fifth  dv))
  (incf *N3-31-BOT* (sixth  dv)))

;; As with the f-vector, the change in the b-vector depends on the
;; move that acts on the spacetime. Unlike the f-vector, the change in
;; the b-vector ALSO depends on where the move occurs. If the move
;; acts on a 2-simplex that is NOT on the the boundary, the b-vector
;; does not change. The only move allowed to act on the boundary is
;; the 23-move and its inverse, the 32-move. Thus these are the only
;; moves which have a chance of changing the b-vector.
(defmacro DB23 (sxid) ; Change in b-vector due to 23-move.
  `(cond ((in-upper-sandwich ,sxid) (list 0 1 0 0 0 0))
	 ((in-lower-sandwich ,sxid) (list 0 0 0 0 1 0))
	 (t                         (list 0 0 0 0 0 0))))

;; change in b-vector due to 32-move. Inverse of DB23
(defmacro DB32 (sxid) 
  `(cond ((in-upper-sandwich ,sxid) (list 0 -1 0 0 0  0))
	 ((in-lower-sandwich ,sxid) (list 0  0 0 0 -1 0))
	 (t                         (list 0  0 0 0  0 0))))

;; The following moves can't affect the boundary. However, I have
;; given them the proper input data anyway, in case we want to test
;; the action when the boundaries are identified.

(defparameter DB44 '(0 0 0 0 0 0)) ; Change in b-vector due to a
				   ; 44-move. 44 is its own inverse.

(defmacro DB26 (sxid) ; Change in b-vector due to a 26-move.
  `(cond ((and (has-face-on-boundary ,sxid) *merge-faces*)
	  (list 3 0 2 3 0 2))
	 (t (list 0 0 0 0 0 0))))

(defmacro DB62 (sxid) ; Change in b-vector due to a 62-move.
  `(cond ((and (has-face-on-boundary, sxid) *merge-faces*)
	  (list -3 0 -2 -3 0 -2))
	 (t (list 0 0 0 0 0 0))))

(defun b-vector()
  "Print the b-vector."
  (format t "[~A, ~A, ~A, ~A, ~A, ~A]~%"
	  *N1-SL-TOP*
	  *N3-22-TOP*
	  *N3-31-TOP*
	  *N1-SL-BOT*
	  *N3-22-BOT*
	  *N3-31-BOT*))

