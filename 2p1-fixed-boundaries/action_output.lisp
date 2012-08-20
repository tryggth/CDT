;;;; action_output.lisp
;;;; Jonah Miller (jonah.maxwell.miller@gmail.com)

;;;; This script is designed to be used with extract_action.py. It
;;;; allows the simulation to print the action of a given spacetime to
;;;; file or to the console.

;;;; Requires pretty much the entire simulation to be loaded before use

(defun current-action-wick-rotated nil
  "Returns the current value of the action for a given spacetime.
   Because it is Wick-rotated, it is completely complex."
  (let ((N3 (+ N3-TL-31 N3-TL-22)))
    (action N1-SL N1-TL N3-TL-31 N3-TL-22
	    *N1-SL-TOP* *N3-22-TOP* *N3-31-TOP*
	    *N1-SL-BOT* *N3-22-BOT* *N3-31-BOT*)))

(defun current-euclidean-action nil
  "Returns the current value of the Euclidean action, which is completely real 
   and negative."
  (* *i* (current-action-wick-rotated)))

(defun current-un-normalized-probability nil
  "Returns the the un-normalized probability of the current spacetime. 
   NOT between 0 and 1."
  (exp (current-euclidean-action)))

(defun output-action (stream)
  (format stream "~A~%" (current-un-normalized-probability)))
  
(defun measure-spacetime-action (filename)
  (with-open-file (f filename)
    (load-spacetime-from-file f))
  (output-action t))
