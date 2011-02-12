(load "cdt3p1.lisp")
(defun setup-33-v2-223 ()
  ;;        4 5 6
  ;; t = 1 ----------
  ;;        1 2 3
  ;; t = 0 ----------
  (let ((ids (list (make-4simplex 2 0 1 1 2 4 5 6)
		   (make-4simplex 2 0 1 1 3 4 5 6)
		   (make-4simplex 3 0 1 1 2 3 4 5))))
    (connect-4simplices-within-list ids))
  (set-last-used-pt 7))

(defun setup-33-v2-332 ()
  ;;        1 2 3
  ;; t = 1 ----------
  ;;        4 5 6
  ;; t = 0 ----------
  (let ((ids (list (make-4simplex 3 0 1 4 5 6 1 2)
		   (make-4simplex 3 0 1 4 5 6 1 3)
		   (make-4simplex 2 0 1 4 5 1 2 3))))
    (connect-4simplices-within-list ids))
  (set-last-used-pt 7))
