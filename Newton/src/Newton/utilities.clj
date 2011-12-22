(ns Newton.utilities
;;  (:use [clojure.test])
  (:use midje.sweet))

;; To keep same name conventions as utilities.lisp
;; In idiomatic clojure, this could be replaced by the anonymous function
;; \#(apply + %)
(defn sum
  "Sums the elements of a list or vector"
  [list]
    (apply + list))

;; Usage/test of sum[...]
(fact (sum '(1 2 3)) => 6)
(fact (sum [1 2 3 4 5]) => 15)

