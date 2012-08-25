#!/usr/bin/env sh

# If you don't have root access on the computer you want to run
# simulations on, you need to follow the instructions in cdt-doc written
# by Rajesh Kommu. CDT/documentation/cdt-doc.
# 
# If you have root access on your computer, you can install emacs,
# slime, and sbcl through the package manager. On ubuntu:
# 
sudo apt-get install sbcl
sudo apt-get install emacs
sudo apt-get install slime
# 
# Find the path to sbcl with
which sbcl
# 
# Add the following to your ~/.emacs file:
# (setq inferior-lisp-program "/path/to/sbcl")
# (add-to-list 'load-path "/usr/share/emacs/site-lisp/slime/")
# (require 'slime)
# (slime-setup)
