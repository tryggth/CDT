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
;; seperately, because the action for these moves is different. This
;; is only counted during initialization.

;; WARNING: this program is not finished. Although the program runs,
;; we don't get correctly optimized lambda values. I'm not sure what's
;; wrong.


;;; Global constants
;;;--------------------------------------------------------------------------
;; Values of little lambda where the probability of acceptace for the
;; action should be completely negative or completely positive
;; respectively. Because I can't figure out a good way to do it, no
;; concern is taken to make sure the action is monotone. Be careful!
(defvar *min-litl* 0.001  "The minimum value little lambda can take.")
(defvar *max-litl* 30 "The maximum value little lambda can take.")

(defvar *p-tester-num-attempts* 1000000 "The number of moves attempted 
 to test probabilities for each move.")

;; Although these constants are not used in the actual simulation
;; (instead we have macros), they're useful here for the action
;; calculation.
(defvar DB23-TOP '(0 1 0 0 0 0) 
  "Change in b-vector for 2->3 move on top boundary")
(defvar DB23-BOT '(0 0 0 0 1 0)
  "Change in b-vector for 2->3 move on bottom boundary")
(defvar DB23-BULK '(0 0 0 0 0 0)
  "Change in b-vector for 2->3 move in bulk.")
(defvar DB32-TOP '(0 -1 0 0 0 0)
  "Change in b-vector for 3->2 move on top boundary.")
(defvar DB32-BOT '(0 0 0 0 -1 0)
  "Change in b-vector for 3->2 move on bottom boundary.")
(defvar DB32-BULK '(0 0 0 0 0 0)
  "Change in b-vector for 3->2 move in bulk.")
;;;--------------------------------------------------------------------------

;;; Function Definitions
;;;--------------------------------------------------------------------------
(defun progress-bar (i i-max)
  "A progress bar for use with incrementing i."
  (let ((percent-done (* 100 (/ i i-max))))
    (when (equalp (mod percent-done 10) 0)
      (format t "-~d" percent-done))))

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
  (min 1 (exp (realpart (* *i* (delta-action action-exposed alpha k litL df db))))))

(defun test-move-probabilities (num-tries)
  "This function randomly attempts (but never applies) moves to the spacetime 
   to see the probability that a given move will be allowed. It then returns
   a list of probability of move acceptance as a function of move. In order,
   the elements are:
   P(2->6 move accepted in bulk)
   P(2->3 move accepted in bulk)
   P(4->4 move accepted in bulk)
   P(3->2 move accepted in bulk)
   P(6->2 move accepted in bulk)
   P(2->3 move accepted in upper boundary sandwich)
   P(3->2 move accepted in upper boundary sandwich)
   P(2->3 move accepted in lower boundary sandwich)
   P(3->2 move accepted in lower boundary sandwich)
   Note that the only moves that are topologically acceptable in the boundary
   sandwich are 2->3 and 3->2.

   The function attempts moves num-tries times."
  (format t "Measuring move probabilities. ~A attempts.~%" num-tries)
  (let ((p-move-attempted (list 0 0 0 0 0 0 0 0 0))
	(successful-attempts 0)
	(simplex-id-list (list-all-3-simplices))
	(probability-list nil))
    (format t "Percent Complete:")
    (for (i 0 num-tries)
      (let* ((move-chooser (random 5))
	     (simplex-chooser (random-element simplex-id-list))
	     (movedata (try-move simplex-chooser move-chooser)))
	(progress-bar i num-tries)
	(when movedata
	  (incf successful-attempts)
	  (cond 
	    ((in-upper-sandwich simplex-chooser)
	     (if (equalp move-chooser 23MTYPE)
		 (incf (nth 5 p-move-attempted))
		 (incf (nth 6 p-move-attempted))))
	    ((in-lower-sandwich simplex-chooser)
	     (if (equalp move-chooser 23MTYPE)
		 (incf (nth 7 p-move-attempted))
		 (incf (nth 8 p-move-attempted))))
	    (t (incf (nth move-chooser p-move-attempted)))))))
    (setf probability-list 
	  (mapcar #'(lambda (x) (/ x num-tries)) p-move-attempted))
    (format t "~%Probabilities are: ~A~%" probability-list)
    probability-list))

(defun p-v-change (action-exposed alpha k litL p-move-attempts)
  "Given an action and coupling constants, generate the pseudo-probability 
that the volume of the spacetime will change. A positive value means the
 probability is for increased volume. A negative value means the probability 
is for decreased volume."
  (let ((p26-bulk (nth 0 p-move-attempts)) ; Bulk
	(p23-bulk (nth 1 p-move-attempts))
	;(p44-bulk (nth 2 p-move-attempts))
	(p32-bulk (nth 3 p-move-attempts))
	(p62-bulk (nth 4 p-move-attempts))
	(p23-top  (nth 5 p-move-attempts)) ; Top boundary sandwich
	(p32-top  (nth 6 p-move-attempts))
	(p23-bot  (nth 7 p-move-attempts)) ; Bottom boundary sandwich
	(p32-bot  (nth 8 p-move-attempts)))
    (- ;; probability of increased volume minus probability of
       ;; decrased volume
     (+ ;; P of increased volume
      (* p23-bulk 
	 (p-accepted action-exposed alpha k litL DF23 DB23-BULK))
      (* p23-top
	 (p-accepted action-exposed alpha k litL DF23 DB23-TOP))
      (* p23-bot
	 (p-accepted action-exposed alpha k litL DF23 DB23-BOT))
      (* 4 p26-bulk
	 (p-accepted action-exposed alpha k litL DF26 DB26)))
     (+ ;; P of decreased volume
      (* p32-bulk 
	 (p-accepted action-exposed alpha k litL DF32 DB32-BULK))
      (* p32-top
	 (p-accepted action-exposed alpha k litL DF32 DB32-TOP))
      (* p32-bot 
	 (p-accepted action-exposed alpha k litL DF32 DB32-BOT))
      (* 4 p62-bulk
	 (p-accepted action-exposed alpha k litL DF62 DB62))))))

(defun tune-litL (action-exposed alpha k min-litL max-litL)
  "Tunes litL to put the simulation on the critical surface (ideally)."
  (format t "Tuning litL. ~3f<=litL<=~3f.~%" max-litL min-litL)
  (let ((p-attempted (test-move-probabilities *p-tester-num-attempts*)))
    (flet ((func (x) (p-v-change action-exposed alpha k x p-attempted)))
      (false-position-root #'func max-litL min-litL (expt 10 -5)))))
;;;--------------------------------------------------------------------------
