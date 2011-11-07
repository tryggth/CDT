#!/bin/bash
#
#$-cwd
#$-j y
#$-S /bin/bash

module load gcc sbcl
hostname
time sbcl --script generate-data-v2-1.0-0.0-100000.lisp