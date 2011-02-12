(load "cdt2p1.lisp")

(defconstant SIGMA-MIN 50)
(defconstant SIGMA-MAX 60)
(defparameter NUM-TRIALS 10000)
(defparameter SPECTRAL-FILENAME "")

;;(defun get-thick-slice-with-maximal-three-volume ()
;;  (let* ((tvps (make-array NUM-T :initial-element 0)) (slice-num -1))
;;    (maphash #'(lambda (sid sx)
;;		 (incf (svref tvps (- (second sx) 1/2))))
;;	     *spacetime*)
;;
;;    (do ((ts 0 (1+ ts)) (v3 0)) ((>= ts NUM-T) slice-num)
;;      (when (>= (svref tvps ts) v3)
;;        (setf slice-num ts)
;;        (setf v3 (svref tvps ts))))))

(defmacro select-neighbor-at-random (sx)
  `(nth (random 4) (3sx-sx3ids ,sx)))

;; random-walk walks sigma steps from sx and returns the simplex at which the walk ends. The returned
;; simplex is guaranteed to be not null
(defun random-walk (sxid sigma)
  (let ((currid sxid)
	(nextid 0))
    (for (s 1 sigma)
	 (while (= 0 nextid)
	   (setf nextid (select-neighbor-at-random (get-3simplex currid))))
	 (setf currid nextid)
	 (setf nextid 0))
    currid))

;; in the approximate version, start at i0 and take NUM-TRIALS random walks of 
;; sigma steps. returns the probability of coming back to i0
(defun KT-approx (i0 sigma)
  (do ((walk-num 0 (1+ walk-num)) (num-ret 0))
      ((>= walk-num NUM-TRIALS) (float (/ num-ret NUM-TRIALS)))
    (if (= i0 (random-walk i0 sigma))
	(incf num-ret))))

(defun spacetime-spectral-dimension (start-tlo start-thi)
  (let* ((idlist (get-simplices-in-sandwich start-tlo start-thi))
         (numids (length idlist))
	 (kt 0.0))
    (for (sigma SIGMA-MIN SIGMA-MAX)
	 (for (count 1 numids)
	      (incf kt (KT-approx (nth (random numids) idlist) sigma)))
	 (with-open-file 
	     (specfile SPECTRAL-FILENAME :direction :output :if-exists :append)
	   (format specfile "~a ~a~%" sigma (/ kt numids)))
	 (setf kt 0.0))))
