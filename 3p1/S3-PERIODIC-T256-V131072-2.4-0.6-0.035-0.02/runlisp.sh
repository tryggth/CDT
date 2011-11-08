#!/bin/bash
#
#$-cwd
#$-j y
#$-S /bin/bash

module load gcc sbcl
hostname
time sbcl --script S3-PERIODIC-T256-V131072-2.4-0.6-0.035-0.02.lisp