#!/bin/bash
#
#$-cwd
#$-j y
#$-S /bin/bash

module load compilers/sbcl-1.0
hostname
time sbcl --script generate-data-v2-0.75-0.50-100000.lisp
