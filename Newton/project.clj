(defproject Newton "1.0.0-SNAPSHOT"
  :description "Causal Dynamical Triangulations"
  :dependencies [[org.clojure/clojure "1.3.0"]]
  :main Newton.core
  :manifest {"Class-Path" "lib/clojure-1.3.0.jar"}
  :dev-dependencies [[lein-marginalia "0.7.0"]]
  :marginalia {:javascript ["mathjax/MathJax.js"]})