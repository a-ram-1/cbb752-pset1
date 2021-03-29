#!/usr/bin/python
__author__ = "FirstName LastName"
__email__ = "first.last@yale.edu"
__copyright__ = "Copyright 2021"
__license__ = "GPL"
__version__ = "1.0.0"

### Usage: python hw1.py -i <input file> -s <score file>
### Example: python hw1.py -i input.txt -s blosum62.txt
### Note: Smith-Waterman Algorithm

import argparse

### This is one way to read in arguments in Python.


parser = argparse.ArgumentParser(description='Smith-Waterman Algorithm')
parser.add_argument('-i', '--input', help='input file', required=True)
parser.add_argument('-s', '--score', help='score file', required=True)
parser.add_argument('-o', '--opengap', help='open gap', required=False,default=-2)
parser.add_argument('-e', '--extgap', help='extension gap', required=False,default=-1)

args = parser.parse_args()
### Implement your Smith-Waterman Algorithm

def runSW(inputFile, scoreFile, openGap, extGap):

### calculation
### write output
### Run your Smith-Waterman Algorithm

runSW(args.input, args.score, args.opengap, args.extgap)

