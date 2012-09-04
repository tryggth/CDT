(TeX-add-style-hook "programmers_guide"
 (lambda ()
    (LaTeX-add-labels
     "s:intro"
     "s:lisp:tricks"
     "s:lisp:tools"
     "s:algorithm"
     "s:data-structures"
     "ss:points"
     "s:links-edges"
     "sss:sl1simplex:id"
     "sss:tl1simplex:id"
     "s:f-and-b-vectors"
     "s:functions"
     "s:f:metropolis:loop"
     "s:initialization"
     "s:initialization:top-level"
     "f:itswvv"
     "s:initialization:mid-level-functions")
    (TeX-run-style-hooks
     "hyperref"
     "color"
     "listings"
     "verbatim"
     "fullpage"
     "latex2e"
     "art12"
     "article"
     "12pt")))

