(ns Newton.test.core
  (:use [Newton.core])
  (:use [midje.sweet])
  (:use [Newton.utilities]))

(facts
  (sum '(1 2 3 4 5)) => 15
  (sum '(1 2 3 4 5 6)) => 21)

;(deftest ^{:utilities true} sum-test
;  (is (= (sum '(1 2 3 4 5)) 15))
;  (is (= (sum [1 2 3 4 5 6]) 21)))

;;(deftest make-3simplices-test
;;  (make-3simplices []))
