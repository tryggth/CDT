(ns Newton.core
  (:gen-class))

(defn -main [msg]
  (println "Hello" msg))

;; # Loading LaTeX
;;
;; When \\(a \ne 0\\), there are two solutions to \\(ax^2 + bx + c = 0\\) and they are
;; $$x = {-b \pm \sqrt{b^2-4ac} \over 2a}.$$
