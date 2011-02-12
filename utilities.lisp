;+---------------------------------------------------------------------------------------------------------+
;| cdt-utilities.lisp --- various utility functions and macros, mostly from Paul Graham's On Lisp
;| Only dimension independent utilites are in this file
;+---------------------------------------------------------------------------------------------------------+
(defmacro circular-nth (n lst)
  `(nth (mod ,n (length ,lst)) ,lst))

(defmacro circular-subseq (seq start num)
  "return num elements from seq starting at 0-based index start"
  `(unless (> ,num (length ,seq))
    (if (<= (+ ,start ,num) (length ,seq))
	(subseq ,seq ,start (+ ,start ,num))
	(let ((part1 (subseq ,seq ,start))
	      (part2 (subseq ,seq 0 (- ,num (- (length ,seq) ,start)))))
	  (append part1 part2)))))

(defmacro take (num lst)
  `(subseq ,lst 0 ,num))

(defmacro drop (num lst)
  `(subseq ,lst ,num))

;; these return the specified element, but in the form of a list
(defmacro firstl (lst)
  `(subseq ,lst 0 1))
(defmacro secondl (lst)
  `(subseq ,lst 1 2))
(defmacro thirdl (lst)
  `(subseq ,lst 2 3))
(defmacro fourthl (lst)
  `(subseq ,lst 3 4))
(defmacro fifthl (lst)
  `(subseq ,lst 4 5))
(defmacro sixthl (lst)
  `(subseq ,lst 5 6))

(defun flatten (x)
  (labels ((rec (x acc)
	     (cond ((null x) acc)
		   ((atom x) (cons x acc))
		   (t (rec (car x) (rec (cdr x) acc))))))
    (rec x nil)))

(defun mapa-b (fn a b &optional (step 1))
  (do ((i a (+ i step))
       (result nil))
      ((> i b) (nreverse result))
    (push (funcall fn i) result)))

(defun map0-n (fn n)
  (mapa-b fn 0 n))

(defun map1-n (fn n)
  (mapa-b fn 1 n))

(defun map-> (fn start test-fn succ-fn)
  (do ((i start (funcall succ-fn i))
       (result nil))
      ((funcall test-fn i) (nreverse result))
    (push (funcall fn i) result)))

(defun mappend (fn &rest lsts)
  (apply #'append (apply #'mapcar fn lsts)))

(defun mapcars (fn &rest lsts)
  (let ((result nil))
    (dolist (lst lsts)
      (dolist (obj lst)
	(push (funcall fn obj) result)))
    (nreverse result)))

(defun rmapcar (fn &rest args)
  (if (some #'atom args)
      (apply fn args)
      (apply #'mapcar
	     #'(lambda (&rest args)
		 (apply #'rmapcar fn args))
	     args)))

(defmacro with-gensyms (syms &body body)
  `(let ,(mapcar #'(lambda (s)
		     `(,s (gensym)))
		 syms)
     ,@body))

(defmacro while (test &body body)
  `(do ()
       ((not ,test))
     ,@body))

(defmacro till (test &body body)
  `(do ()
       (, test)
     ,@body))

(defmacro repeat-until (test &body body)
  `(do ()
       (, test)
     ,@body))

(defmacro until (test &body body)
  `(do ()
       (, test)
     ,@body))

(defmacro for ((var start stop &optional (by 1)) &body body)
  (let ((gstop (gensym)))
    `(do ((,var ,start (incf ,var ,by))
	  (,gstop ,stop))
	 ((> ,var ,gstop))
       ,@body)))

(defmacro do-tuples/o (parms source &body body)
  (if parms
      (let ((src (gensym)))
	`(prog ((,src ,source))
	    (mapc #'(lambda ,parms ,@body)
		  ,@(map0-n #'(lambda (n)
				`(nthcdr ,n ,src))
			    (1- (length parms))))))))

(defun dt-args (len rest src)
  (map0-n #'(lambda (m)
	      (map1-n #'(lambda (n)
			  (let ((x (+ m n)))
			    (if (>= x len)
				`(nth ,(- x len) ,src)
				`(nth ,(1- x) ,rest))))
		      len))
	  (- len 2)))

(defmacro do-tuples/c (parms source &body body)
  (if parms
      (with-gensyms (src rest bodfn)
	(let ((len (length parms)))
	  `(let ((,src ,source))
	     (when (nthcdr ,(1- len) ,src)
	       (labels ((,bodfn ,parms ,@body))
		 (do ((,rest ,src (cdr ,rest)))
		     ((not (nthcdr ,(1- len) ,rest))
		      ,@(mapcar #'(lambda (args)
				    `(,bodfn ,@args))
				(dt-args len rest src))
		      nil)
		   (,bodfn ,@(map1-n #'(lambda (n)
					 `(nth ,(1- n)
					       ,rest))
				     len))))))))))

;; (unions '(1 2 3) '(2 3 4) '(3 4 1)) => (1 2 3 4)
(defun unions (&rest args)
  (reduce #'union args))

;; (intersections '(1 2 3) '(2 3 4) '(3 4 1)) => (3)
(defun intersections (&rest args)
  (reduce #'intersection args))

;; (differences '(1 2 3 4 5 6) '(2 3) '(4 5)) => (1 6) 
(defun differences (&rest args)
  (reduce #'set-difference args))

;; (set-equal? '(2 3 4 1) '(4 1 2 3)) => T
(defun set-equal? (set1 set2)
  (and (eql (set-difference set1 set2) nil) (eql (set-difference set2 set1) nil)))

(defmacro yyyy (lst) `(nth 5 ,lst))
(defmacro mm (lst) `(nth 4 ,lst))
(defmacro dd (lst) `(nth 3 ,lst))
(defmacro hh (lst) `(nth 2 ,lst))
(defmacro mi (lst) `(nth 1 ,lst))
(defmacro ss (lst) `(nth 0 ,lst))

(defun cdt-now-str ()
  (let ((nowlst (multiple-value-list (get-decoded-time))))
    (format nil "~4,'0d-~2,'0d-~2,'0d ~2,'0d:~2,'0d:~2,'0d" 
	    (yyyy nowlst) (mm nowlst) (dd nowlst) (hh nowlst) (mi nowlst) (ss nowlst))))