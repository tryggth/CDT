(ns Newton.test.utilities
  (:use [midje.sweet])
  (:use [Newton.utilities]))

(facts
  (sum '(1 2 3 4 5)) => 15
  (sum '(1 2 3 4 5 6)) => 21)
