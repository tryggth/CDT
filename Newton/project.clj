(defproject Newton "1.0.0-SNAPSHOT"
  :description "Causal Dynamical Triangulations in Clojure with Newtonian Approximation"
  :url "https://github.com/ucdavis/CDT"
  :dependencies [[org.clojure/clojure "1.3.0"]]
  :test-selectors {:default (fn [v] (not (:utilities v)))
                   :utilities :utilities
                   :all (fn [_] true)})