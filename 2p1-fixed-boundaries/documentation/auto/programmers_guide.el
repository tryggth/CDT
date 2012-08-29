(TeX-add-style-hook "programmers_guide"
 (lambda ()
    (LaTeX-add-labels
     "s:intro"
     "s:algorithm"
     "s:data-structures"
     "ss:points"
     "sss:sl1simplex:id"
     "sss:tl1simplex:id"
     "s:functions")
    (TeX-run-style-hooks
     "color"
     "listings"
     "verbatim"
     "fullpage"
     "latex2e"
     "art12"
     "article"
     "12pt")))

