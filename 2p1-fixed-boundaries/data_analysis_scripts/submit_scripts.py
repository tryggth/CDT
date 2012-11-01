#!/usr/bin/env python

"""
Author: Jonah Miller (jonah.maxwell.miller@gmail.com)
Time-stamp: <2012-10-31 20:50:44 (jonah)>

This script is designed to be run on Janus. It takes a number of
scripts (passed in by command line) meant to be run by SBCL and writes
a shell script for each one. It then submits them all. 
