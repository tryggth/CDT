(ns Newton.utilities)

;; To keep same name conventions as utilities.lisp
;; In idiomatic clojure, this could be replaced by the anonymous function
;; #(apply + %)
(defn sum
  "sums the elements of a list"
  [list]
    (apply + list))

