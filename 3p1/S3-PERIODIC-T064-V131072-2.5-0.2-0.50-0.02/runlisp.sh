#!/bin/bash
#
#$-cwd
#$-j y
#$-S /bin/bash

module load compilers/sbcl-1.0
hostname
time sbcl --script S3-PERIODIC-T064-V131072-2.5-0.2-0.50-0.02.lisp
