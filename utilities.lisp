(declaim (optimize (speed 3)
		   (compilation-speed 0)
		   (debug 0)
		   (safety 0)))
;+-----------------------------------------------------------------------------+
;| cdt-utilities.lisp --- various utility functions and macros, mostly from 
;| Paul Graham's "On Lisp". Only dimension independent utilites are in this file
;+-----------------------------------------------------------------------------+
(defmacro sum (lst)
  "sums the elements of lst"
  `(apply #'+ ,lst))

(defmacro filter (func lst)
  `(mapcan (lambda (y) (when (,func y) (list y))) ,lst)
)

(defmacro random-element (lst)
  "returns a random element from lst"
  `(nth (random (length ,lst)) ,lst))

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

(defmacro 2+ (num)
  `(+ ,num 2))

(defmacro nthl (n lst)
  "same as nth, but element is returned as (val) rather than val"
  `(subseq ,lst ,n (1+ ,n)))
(defmacro firstl (lst)
  `(nthl 0 ,lst))
(defmacro secondl (lst)
  `(nthl 1 ,lst))
(defmacro thirdl (lst)
  `(nthl 2 ,lst))
(defmacro fourthl (lst)
  `(nthl 3 ,lst))
(defmacro fifthl (lst)
  `(nthl 4 ,lst))
(defmacro sixthl (lst)
  `(nthl 5 ,lst))

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

;; These next two methods retrieve useful sublists of lists. 
;; They return _brand-new_ lists.
(defun pairs (input)
  "Get all pairs of distinct elements from list."
  (loop for i on input append (loop for j in (cdr i) for k = (car i) collect (list k j))))

(defun triples-decompose (triple)
  "Carry (p0 p1 p2) to ((p0 (p1 p2)) (p1 (p0 p2)) (p2 (p0 p1))). Order matters."
  (let ((p0 (elt triple 0))
	(p1 (elt triple 1))
	(p2 (elt triple 2)))
    (list (list p0 (list p1 p2))
	  (list p1 (list p0 p2))
	  (list p2 (list p0 p1)))))

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
  (and (eql (set-difference set1 set2) nil) 
       (eql (set-difference set2 set1) nil)))

(defun cdt-now-str ()
  "returns the current time in the YYYY-MM-DD-hh-mm-ss format"
  (multiple-value-bind (ss mi hh dd mm yyyy) (get-decoded-time)
    (format nil "~4,'0d-~2,'0d-~2,'0d-~2,'0d-~2,'0d-~2,'0d" 
	    yyyy mm dd hh mi ss)))

(defun hostname ()
  (let* ((machinst (machine-instance))
	 (dotpos (position #\. machinst)))
    (if dotpos
	(subseq machinst 0 dotpos)
	machinst)))

(defun change-file-suffix (filename newsuffix)
  "change-file-suffix replaces the old suffix of filename with the newsuffix.
Suffix is defined as the string following the LAST . in filename"
  (let* ((dotpos (position #\. filename :from-end t))
	 (prefix (subseq filename 0 dotpos))
	 (cpy (concatenate 'string 
			   prefix
			   "."
			   (make-string (length newsuffix) 
					:initial-element #\space))))
    (replace cpy newsuffix :start1 (1+ dotpos))))

(defun nthhash (n h)
  "returns element number n of hashtable h"
  (maphash #'(lambda (k v) 
	       (when (zerop n) 
		 (if (null v)
		     (format t "returning a nuuuuuuuuuul~%"))
		 (return-from nthhash (values v k)))
	       (decf n)) 
	   h))

(defun rndhash (h)
  "returns a random element from the hashtable h"
  (if (= 0 (hash-table-count h))
      (format t "we have a problem~%"))
  (nthhash 0 h))
  ;;(nthhash (random (hash-table-count h)) h))

;; Access table using key and try to push val to the list
;; found there. Make a list (val) if necessary
(defun gethash-push (key table val)
  (if (gethash key table) 
      (push val (gethash key table))
      (setf (gethash key table) (list val))))

;; Access table using key and remove val from the list
;; found there. Delete the list if it is left empty.
(defun gethash-delete (key table val)
  (unless (setf (gethash key table) (delete val (gethash key table)))
    (remhash key table)))

(defun gethash-update (key table vals &optional (index nil))
  "Index lets us update a field of the target, rather than the whole target."
  (let ((cur (gethash key table))
	(out 0))
    (when index
      (setf cur (elt cur index)))
    (if cur
	(setq out (mapcar (lambda (x y) (+ x y)) cur vals))
	(setq out vals)
	)
    (if (> (length (filter (lambda (x) (/= x 0)) out)) 0)
	(if index
	    (setf (elt (gethash key table) index) out)
	    (setf (gethash key table) out))
	(unless index
	  (remhash key table)))))

(defmacro zero-or (coupling term)
  `(if (/= ,coupling 0) ,term 0))

;; Numerical root finding

(defun false-position-root (func pos neg precision &optional fpos fneg)
  (unless fpos
    (setf fpos (funcall func pos)))
  (unless fneg
    (setf fneg (funcall func neg)))
  (format t "Pos: ~A, Neg: ~A ~%" fpos fneg)
  (if (< (- fneg) precision)
      neg
      (if (< fpos precision)
	  pos
	  (let* ((new (/ (- (* fpos neg) (* fneg pos))
			 (- fpos fneg)))
		 (fnew (funcall func new)))
	    (if (= 0 fnew)
		new
		(if (> fnew 0)
		    (false-position-root func new neg precision fnew fneg)
		    (false-position-root func pos new precision fpos fnew)))))))

(defun false-position-root-near (func guess precision)
  (let* ((fguess (funcall func guess))
	 (desired (if (> fguess 0) #'< #'>))
	 (other-del 0)
	 (fworking 0))
    (when (< (abs fguess) precision)
      (return-from false-position-root-near guess))
    (setf other-del (if (funcall desired (funcall func (+ guess 1)) fguess) 1 -1))
    (until (funcall desired (setf fworking (funcall func (+ guess other-del))) 0)
	   (setf other-del (* other-del 2)))
    (if (> fguess 0)
	(progn
	  (assert (< fworking 0))
	  (false-position-root func guess (+ guess other-del) precision fguess fworking))
	(progn
	  (assert (> fworking 0))
	  (false-position-root func (+ guess other-del) guess precision fworking fguess)))))
	 