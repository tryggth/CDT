(TeX-add-style-hook "_region_"
 (lambda ()
    (LaTeX-add-labels
     "s:intro"
     "s:algorithm"
     "s:data-structures"
     "ss:points"
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

