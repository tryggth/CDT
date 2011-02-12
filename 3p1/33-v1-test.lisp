(load "cdt3p1.lisp")
(defun setup-33-v1-122 ()
  ;;        3 4 5 6
  ;; t = 1 ----------
  ;;        1 2 7
  ;; t = 0 ----------
  (let ((ids (list (make-4simplex 1 0 1 1 3 4 5 6)
		   (make-4simplex 2 0 1 1 2 3 4 5)
		   (make-4simplex 2 0 1 1 2 3 4 6)
		   (make-4simplex 3 0 1 1 2 7 4 5))))
    (connect-4simplices-within-list ids))
  (set-last-used-pt 7))

(defun setup-33-v1-433 ()
  ;;        1 2 7
  ;; t = 1 ----------
  ;;        3 4 5 6
  ;; t = 0 ----------
  (let ((ids (list (make-4simplex 4 0 1 3 4 5 6 1)
		   (make-4simplex 3 0 1 3 4 5 1 2)
		   (make-4simplex 3 0 1 3 4 6 1 2)
		   (make-4simplex 2 0 1 4 5 1 2 7))))
    (connect-4simplices-within-list ids))
  (set-last-used-pt 7))
