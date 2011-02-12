(load "cdt3p1.lisp")
;; t = 1 3,4,5,6,7
;; t = 0 1,2
(let ((ids (list (make-4simplex 2 0 1 1 2 6 4 5)
		 (make-4simplex 3 0 1 3 1 2 5 4)
		 (make-4simplex 1 0 1 1 4 5 6 7)
		 (make-4simplex 4 0 1 1 2 3 8 5))))
  (connect-4simplices-within-list ids))
(set-last-used-pt 8)
