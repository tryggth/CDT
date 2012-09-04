(TeX-add-style-hook "_region_"
 (lambda ()
    (LaTeX-add-labels
     "s:intro"
     "s:algorithm"
     "s:data-structures"
     "ss:points"
     "s:links-edges"
     "sss:sl1simplex:id"
     "sss:tl1simplex:id"
     "s:f-and-b-vectors"
     "s:functions")
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

