(ns Newton.simplex
  (:use midje.sweet))

;; # Load Optional Resources
;; Use external Javascript and CSS in your documentation. For example:
;; To format Latex math equations, download the
;; [MathJax](http://www.mathjax.org/) Javascript library to the docs
;; directory and then add
;;
;;     :marginalia {:javascript ["mathjax/MathJax.js"]}
;;
;; to project.clj.
;;
;; When \\(a \ne 0\\), there are two solutions to \\(ax^2 + bx + c = 0\\) and they are
;; $$x = {-b \pm \sqrt{b^2-4ac} \over 2a}.$$

;; A simplex is a d-dimensional triangle. It has d+1 points and
;; $\left(\begin{array}{c}d+1\\2\end{array}\right)$ edges
;; A face is denoted by 3 edges; each edge has a point in common with another edge,
;; and the edges should show cyclic permutation. Thus, for a set of 3 edges, there should only
;; be 3 distinct points.
;; Two simplices share a common face. They are then said to be neighbors.


(unfinished type-id)

(defn simplex?
  "Tests that a simplex satisfies the following criteria:
    every vertex has an edge with every other vertex
    at least one vertex is on t-lo
    at least one vertex is on t-hi
    all faces are defined"
  [& {type-id :type t-lo :time-lo t-hi :time-high points :vertices}]



  )

(fact
  (simplex? :type type-id
    :time-lo 0
    :time-high 1
    :vertices [0 1 2 3 4]) => truthy)

;;"failing scenario too many or too few points"
;;"abstract to d-simplex?"