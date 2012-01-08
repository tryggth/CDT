(ns Newton.test.simplex
  (:use [midje.sweet])
  (:use [Newton.simplex]))

(unfinished simplex?)

"simplex-id type tlo thi vertices neighbors"
"vertex-id time-slice"
"type should be deducible from vertices"
"vertices should be sets"

(defn simplex?
  [& {id :simplex-id type-id :type t-lo :time-lo t-hi :time-high points :vertices neighbors :neighbors}]
  )

;(fact
;  (animal-set (animals)) => #{"cow" "horse"}
;  (provided
;    (animals) => ["cow" "horse"]))


"failing scenario - not a simplex"


;(unfinished type-id)
;
;(defn simplex?
;  "Tests that a simplex satisfies the following criteria:
;    every vertex has an edge with every other vertex
;    at least one vertex is on t-lo
;    at least one vertex is on t-hi
;    all faces are defined"
;  [& {type-id :type t-lo :time-lo t-hi :time-high points :vertices neighbors :neighbors}]
;
;
;
;  )
;
;(fact
;  (simplex? :type type-id
;    :time-lo 0
;    :time-high 1
;    :vertices [0 1 2 3 4]) => truthy)

;;"failing scenario too many or too few points"
;;"abstract to d-simplex?"
