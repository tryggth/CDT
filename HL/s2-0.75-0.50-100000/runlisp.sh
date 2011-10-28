#!/bin/bash
#
#$-cwd
#$-j y
#$-S /bin/bash

module load gcc sbcl
hostname
time sbcl --script generate-data-v2-0.75-0.50-100000.lisp --dynamic-space-size 200
