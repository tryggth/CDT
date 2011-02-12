(load "cdt3p1.lisp")


(let ((ids (list (make-4simplex 1 0 1 1 2 3 4 5)
		 (make-4simplex 1 0 1 1 3 4 5 6)
		 (make-4simplex 4 1 2 2 3 4 5 7)
		 (make-4simplex 4 1 2 3 4 5 6 7)
		 (make-4simplex 2 0 1 1 8 2 3 5)
		 (make-4simplex 2 0 1 1 9 2 3 4))))
  (connect-4simplices-within-list ids))
(set-last-used-pt 9)