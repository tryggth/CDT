(defproject Newton "1.0.0-SNAPSHOT"
  :description "Newtonian Approximation for Causal Dynamical Triangulations in Clojure"
  :url "https://github.com/ucdavis/CDT"
  :dependencies [[org.clojure/clojure "1.3.0"]]
  :marginalia {:javascript["mathjax/MathJax.js"]}
  :test-selectors {:default (fn [v] (not (:utilities v)))
                   :utilities :utilities
                   :all (fn [_] true)}
  :dev-dependencies [[lein-ring "0.4.5"]
                     [lein-marginalia "0.6.1"]
                     [midje "1.3.1"]
                     [lein-midje "1.0.7"]
                     [com.stuartsierra/lazytest "1.2.3"]]
  :repositories {"stuart" "http://stuartsierra.com/maven2"})