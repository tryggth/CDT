;;;; tuning_utilities.lisp
;;;; Jonah Miller (jonah.maxwell.miller@gmail.com)

;;;; This script is designed to be used with map_phase_space.py. It
;;;; helps the simulation take data on volume variation.

(defun run-thermalization-sweeps (sweeps)
  "Does several sweeps without taking any data on them.
   Allows for the system volume to thermalize, which takes time."
  (for (i 0 sweeps) 
    (sweep)))

(defun collect-volume-data (sweeps)
  "Similar to the collect data functions in montecarlo.lisp, 
but designed to test the phase space."
  (let ((3-volumes nil))
    (for (ns 0 sweeps)
	 (sweep)
	 (push (+ N3-TL-31 N3-TL-22) 3-volumes))
    3-volumes))

(defun get-delta (list-of-values)
  "Calculate a crude derivative for a
 list of values (as a function of list index."
  (* -1 (- (nth (1- (length list-of-values)) list-of-values) 
	   (first list-of-values))))

(defun test-volume-variation (sweeps &optional (k-and-litl nil))
  "Runs collect-volume-data and then does statistics on it, and prints it to 
  the terminal."
  (let* ((volumes (collect-volume-data sweeps))
	 (avg (mean volumes))
	 (std (standard-deviation volumes))
	 (deltav (get-delta volumes))
	 (datalist (list avg std deltav)))
    (if k-and-litl 
	(format t "~A ~A ~A ~A ~A~%" *k* *litl* avg std deltav)
	(format t "~A ~A ~A ~A ~A~%" *k0* *k3* avg std deltav))
    datalist))
	 
(defun add-data-to-file (outfile n3-mean n3-sigma n3-delta 
			 &optional (k-and-litl nil))
  "Adds the requisite information to the output file."
  (with-open-file (f outfile
		     :direction :output
		     :if-exists :append
		     :if-does-not-exist :create)
    (if k-and-litl
	(format f "~3D~T~3D~T~3D~T~3D~T~3D~%" *k* *litl* n3-mean n3-sigma 
		n3-delta)
	(format f "~3D~T~3D~T~3D~T~3D~T~3D~%" *k0* *k3* n3-mean n3-sigma 
		n3-delta))))

(defun take-data-and-print-to-file (outfile sweeps)
  "Take volume data and then print it to the given file."
    (let* ((volume-data (test-volume-variation sweeps))
	   (avg (nth 0 volume-data))
	   (std (nth 1 volume-data))
	   (deltav (nth 2 volume-data)))
      (add-data-to-file outfile avg std deltav)))


(defun test-random-number-generator (func max iterations)
  "Tests random number generator for comparison to analytic solutions."
  (let ((numbers nil)
	(mean-std (list 0 0)))
    (for (i 0 iterations)
      (push (funcall func max) numbers))
    (setf (first mean-std) (mean numbers))
    (setf (second mean-std) (standard-deviation numbers))
    (format t "Mean: ~A~%Standard: ~A~%." (first mean-std) (second mean-std))
    mean-std))
