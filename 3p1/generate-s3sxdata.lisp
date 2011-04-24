(load "../../utilities.lisp")
(load "../globals.lisp")
(load "../simplex.lisp")

(with-open-file (infiles "4sx3p1002")
  (generate-s3sx3p1-files infiles))
