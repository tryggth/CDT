;;;; parameter_tuning.lisp
;;;; Author: Jonah Miller (jonah.maxwell.miller@gmail.com)
;;;; July 19, 2012.

;;;; Given a k, this program will solve for a value of litL (little
;;;; lambda) in the Einstein-Hilbert action for 2+1-dimensional
;;;; CDT. This should work for other actions as well, but program will
;;;; need to be tweaked since the number of input terms in the action,
;;;; and the min and maximum values of the coupling constants (and
;;;; probabilities) may change. *alpha* is assumed to be 1.

;;;; Honestly, this could probably be done much more effectively and
;;;; much more generally by a lisp guru than by me, but... I'm doing
;;;; my best. Sorry, next guy!

;;;; This program is assumed to be loaded along with the rest of the
;;;; simulation.

;;;; SOME WEIRD CONCERNS
;; As a crude estimation of how many moves of each type are attempted,
;; we weight each probability of acceptance based on the type of move
;; attempted in each case. This is the PATTEMPTED variable for us. As
;; a result, run these algorithms only after initialization. It's
;; important to note that the boundary moves have to be counted
;; seperately, because the action for these moves is different.

;;; Global constants
;;;--------------------------------------------------------------------------
;; Values of little lambda where the probability of acceptace for the
;; action should be completely negative or completely positive
;; respectively. Because I can't figure out a good way to do it, no
;; concern is taken to make sure the action is monotone. Be careful!
(defvar min-litl 0.001  "The minimum value little lambda can take.")
(defvar max-litl 30 "The maximum value little lambda can take.")
(defvar db-no-boundary '(0 0 0 0 0 0) "We assume no change to the boundary, 
because the probability of a change to the boundary should be minimal.")
;;;--------------------------------------------------------------------------

;;; Function Definitions
;;;--------------------------------------------------------------------------
(defun delta-action (action-exposed alpha k litL df db)
  "Given an action and coupling constants, calculate the change 
in the action for a given change in the f-vector (df) and a given b-vector, 
db."
  (let* (;; BULK
	 ;(d-n0        (nth 0 df))  ; Number of points
	 (d-n1-sl     (nth 1 df))  ; Number of space-like links
	 (d-n1-tl     (nth 2 df))  ; Number of time-like links
	 ;(d-n2-sl     (nth 3 df))  ; Number of space-like triangles
	 ;(d-n2-tl     (nth 4 df))  ; Number of time-like triangles
	 (d-n3-tl-31  (nth 5 df))  ; Number of (3,1)- and (1,3)-simplices
	 (d-n3-tl-22  (nth 6 df))  ; Number of (2,2)-simplices
	 ;; BOUNDARY
	 (d-n1-sl-top (nth 0 db))  ; Number of space-like links in
				   ; upper boundary
	 (d-n3-22-top (nth 1 db))  ; Number of (2,2)-simplices
				   ; connected to the upper boundary
	 (d-n3-31-top (nth 2 db))  ; Number of (1,3)-simplices
				   ; connected to the upper boundary
	 (d-n1-sl-bot (nth 3 db))  ; Number of space-like links in the
				   ; lower boundary
	 (d-n3-22-bot (nth 4 db))  ; Number of (2,2)-simplices
				   ; connected to the lower boundary
	 (d-n3-31-bot (nth 5 db))) ; Number of (3,1)-simplices
				   ; connected to the lower boundary.
    (funcall action-exposed 
	     d-n1-sl d-n1-tl d-n3-tl-31 d-n3-tl-22
	     d-n1-sl-top d-n3-22-top d-n3-31-top
	     d-n1-sl-bot d-n3-22-bot d-n3-31-bot
	     alpha k litL)))

(defun p-accepted (action-exposed alpha k litL df db)
  "Given an action and coupling constants, generate the probability a given 
move will be accepted. Information for the move is in the f and b vectors."
  (exp (realpart (* *i* (delta-action action-exposed alpha k litL df db)))))


(defun p-v-change (action-exposed alpha k litL)
  "Given an action and coupling constants, generate the pseudo-probability 
that the volume of the spacetime will change. A positive value means the
 probability is for increased volume. A negative value means the probability 
is for decreased volume."
  (- ;; probability of increased volume minus probability of decrased volume 
   (+ ;; P of increased volume
    (p-accepted action-exposed alpha k litL DF23 db-no-boundary)
    (* 4 (p-accepted action-exposed alpha k litL DF26 db-no-boundary)))
   (+ ;; P of decreased volume
    (p-accepted action-exposed alpha k litL DF32 db-no-boundary)
    (* 4 (p-accepted action-exposed alpha k litL DF62 db-no-boundary)))))
;;;--------------------------------------------------------------------------
