;;;; action_output.lisp
;;;; Jonah Miller (jonah.maxwell.miller@gmail.com)

;;;; This script is designed to be used with extract_action.py. It
;;;; allows the simulation to print the action of a given spacetime to
;;;; file or to the console.

;;;; Requires pretty much the entire simulation to be loaded before use

(defun current-action-wick-rotated nil
  "Returns the current value of the action for a given spacetime.
   Because it is Wick-rotated, it is completely complex."
  (action N1-SL N1-TL N3-TL-31 N3-TL-22
	  *N1-SL-TOP* *N3-22-TOP* *N3-31-TOP*
	  *N1-SL-BOT* *N3-22-BOT* *N3-31-BOT*))

(defun current-euclidean-action nil
  "Returns the current value of the Euclidean action, which is completely real 
   and negative."
  (let ((euclidean-action (* *i* (current-action-wick-rotated))))
    (assert (zerop (imagpart euclidean-action)))
    (realpart euclidean-action)))

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

(defun mean-and-std-dev-probability-amplitude (filename-list stream)
  "Runs measure-spacetime-action on every file in filename-list and prints 
   the mean and standard deviation."
  (let ((probability-amplitudes nil)
	(mean nil)
	(std nil))
    (loop for file in filename-list do
	 (reset-spacetime-fast)
	 (with-open-file (f file)
	   (load-spacetime-from-file f))
	 (push (* -1 (current-euclidean-action)) probability-amplitudes))
    (setf mean (mean probability-amplitudes))
    (setf std  (standard-deviation probability-amplitudes))
    (format stream "Mean: ~E Std: ~E~%" mean std)
    (values mean std)))
	 
