(load "spacetime-spectral-dimension.lisp")

(with-open-file (infile 
		 "S2-PERIODIC-T016-V016384-1.0-0.85-0.02-000100001-000200000.3sx2p1" 
		 :direction :input)
  (load-spacetime-from-file infile))


;; sts2p1 stands for sPACEtIMEsPECTRAL in the 2pLUS1 case
;; sss2p1 stands for sPATIALsLICEsPECTRAL in the 2pLUS1 case
;; sth2p1 spacetime hausdorff
;; ssh2p1 spatial slice hausdorff

(setf SPECTRAL-FILENAME (concatenate 'string (generate-filename 100001 200000) ".sts2p1"))
; open and close the file so append works later on
(with-open-file (specfile 
		 SPECTRAL-FILENAME
		 :direction :output :if-exists :supersede))

;; somehow figure out the slice with most number of simplices -- say this is identified by tlo and thi
;; then call the following
(spacetime-spectral-dimension 15 16)