(ns Newton.test.core
  (:use [Newton.core])
  (:use [clojure.test])
  (:use [Newton.utilities]))

(deftest ^{:utilities true} sum-test
  (is (= (sum '(1 2 3 4 5)) 15))
  (is (= (sum [1 2 3 4 5 6]) 21)))

;;(deftest make-3simplices-test
;;  (make-3simplices []))
