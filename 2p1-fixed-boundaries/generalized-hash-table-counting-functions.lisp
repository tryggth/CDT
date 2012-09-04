;;;; generalized-hash-table-counting-functions.lisp
;;;; Author: Jonah MIller (jonah.maxwell.miller@gmail.com)


;;; What follows are a number of functions to retrieve information
;;; about simplexes, by listing them or counting them.

;; A general function to list simplex IDs with some trait. 
;; the "trait" input is a function that returns a boolean value. Thus you can
;; check to see if some property of the element of the list that is the key
;; for the simplex returns true.
;;
;; For instance, if you want the list of all space-like 2-simplices at
;; time zero:
;;
;; (list-keys-with-trait #'(lambda (x) (= x 0)) *SL2SIMPLEX->ID* 0)
;;
;; This checks each key in the hashtable to see if the zeroth 
;; element of the key is equal to zero.
;; 
;; If you want to merely count all keys in the hash table, let key-subindex=0
;; and let trait=#'(lambda (x) 0). The lambda function actually just needs to
;; return a non-nil value.
(defun list-keys-with-trait (trait hashtable key-subindex)
  (let ((keylist nil)
	(vallist nil))
    (flet ((discriminator (hkey hval)
	     (when (funcall trait (nth key-subindex hkey))
	       (push hkey keylist)
	       (push hval vallist))))
      (maphash #'discriminator hashtable)
      keylist)))

;; A general function with similar functionality to
;; list-keys-with-trait. Instead of listing the simplexes, however,
;; this one simply counts them, and returns a value.
(defun count-keys-with-trait (trait hashtable key-subindex)
  (let ((count 0)
	(vallist nil))
    (flet ((discriminator (hkey hval)
	     (when (funcall trait (nth key-subindex hkey))
	       (incf count 1)
	       (push hval vallist))))
      (maphash #'discriminator hashtable)
      count)))

;; Similar to list-keys-with-trait but works for values instead
(defun list-vals-with-trait (trait hashtable key-subindex)
  (let ((keylist nil)
	(trashlist nil))
    (flet ((discriminator (hkey hval)
	     (when (funcall trait (nth key-subindex hval))
	       (push hval keylist)
	       (push hkey trashlist))))
      (maphash #'discriminator hashtable)
      keylist)))

;; Similar to count-keys-with-trait but works for values instead 
(defun count-vals-with-trait (trait hashtable key-subindex)
  (let ((count 0)
	(trashlist nil))
    (flet ((discriminator (hkey hval)
	     (when (funcall trait (nth key-subindex hval))
	       (incf count 1)
	       (push hkey trashlist))))
      (maphash #'discriminator hashtable)
      count)))


;; Function to count objects in a single time slice over the entire
;; spacetime. I know this could be written as a macro, but it works
;; passably as a function, so we'll do that. It's true that if the
;; boundary conditions are periodic, this count will look at the empty
;; NUM-T time slice, but since the slice is empty (and there's no
;; periodic boundary conditions modulo going on here), the result of
;; the count there will be zero.
(defun count-over-all-spacetime-slices (counting-function 
					&optional (count-argument nil))
  "Runs a counting function of a single time slice over all spacetime
slices and returns the sum. If the spacetime. If a count-argument is
given, this is passed to the counting-function as the second argument
after time-slice."
  (let ((count 0)
	(tmax (if (string= BCTYPE "PERIODIC")
		  (1- NUM-T)
		  NUM-T)))
    (if (not count-argument)
	(loop for time-slice from 0 to tmax do
	     (setf count (+ count (funcall counting-function time-slice))))
	(loop for time-slice from 0 to tmax do
	     (setf count (+ count (funcall 
				   counting-function 
				   time-slice 
				   count-argument)))))
    count))
       

;; Function to count objects in a sandwich over the entire
;; spacetime. I know this could be written as a macro, but it works
;; passably as a function, so we'll do that. It's true that if the
;; boundary conditions are periodic, this count will look at the empty
;; NUM-T-1, NUM-T, sandwich, but since the slice is empty (and there's no
;; periodic boundary conditions modulo going on here), the result of
;; the count there will be zero.
(defun count-over-all-spacetime-sandwiches (counting-function 
					    &optional (count-argument nil))
  "Runs a counting function of a single sandwich over all sandwiches
in the spacetime and returns the sum. If a count-argument is given,
this is passed to the counting-function as the second argument after
time-slice."
  (let ((count 0))
    (if (not count-argument)
	(loop for time-slice from 0 to (1- NUM-T) do
	     (setf count (+ count (funcall counting-function 
					   time-slice
					   (1+ time-slice)))))
	(loop for time-slice from 0 to (1- NUM-T) do
	     (setf count (+ count (funcall 
				   counting-function 
				   time-slice 
				   (1+ time-slice)
				   count-argument)))))
    count))
