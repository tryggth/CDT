#!/bin/bash
#
#$-cwd
#$-j y
#$-S /bin/bash

module load gcc sbcl
hostname
time sbcl --script S3-PERIODIC-T256-V131072-2.0-0.2-0.407-0.02.lisp