(load "cdt2p1.lisp")

(defmacro yyyy (lst) `(nth 5 ,lst))
(defmacro mm (lst) `(nth 4 ,lst))
(defmacro dd (lst) `(nth 3 ,lst))
(defmacro hh (lst) `(nth 2 ,lst))
(defmacro mi (lst) `(nth 1 ,lst))
(defmacro ss (lst) `(nth 0 ,lst))

(defun now ()
  (let ((nowlst (multiple-value-list (get-decoded-time))))
    (format nil "~4,'0d-~2,'0d-~2,'0d ~2,'0d:~2,'0d:~2,'0d" 
	    (yyyy nowlst) (mm nowlst) (dd nowlst) (hh nowlst) (mi nowlst) (ss nowlst))))

(defun convert-data (filename)
  (format t "starting conversion for ~A at ~A~%" filename (now))
  (with-open-file (infile filename :direction :input)
    (load-spacetime-from-file infile))
  (format t "finished loading old file for ~A at ~A~%" filename (now))
  (with-open-file (outfile (concatenate 'string filename ".newfmt") :direction :output)
    (save-spacetime-to-file-v2 outfile))
  (format t "finished saving new file for ~A at ~A~%" filename (now)))


(convert-data "S2-PERIODIC-T008-V065608-1.0-0.78-0.02-000000001-000200000.3sx2p1")