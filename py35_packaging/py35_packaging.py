'''
Python packaging_test - simple Python packaging test
================================================

:Author: Antonio Berlanga-Taylor
:Release: |date|
:Date: |today|


Purpose
=======

This script is a dummy python package test.

Usage and Options
=================

Usage:
       py35_packaging.py [-h | --help]
       py35_packaging.py [--version]
       py35_packaging.py [--quiet]
       py35_packaging.py [--verbose]
       py35_packaging.py [--log=<log_file> | -L <log_file>]

Options:
    -h --help                     Show this screen.
    --version                     Show version.
    --quiet                       Print less text.
    --verbose                     Print more text.
    -L FILE --log=FILE            Log file name. [default: py35_packaging.log]

Documentation
=============

'''

# Python modules:
import sys
import os

# Modules not in core library:
import docopt

# Package module:
#from py35_packaging import xxx 
#import pyPackaging

# Get package source directory in (param path) '
#src_dir = pyPackaging.getDir('..')
#print('Python packaging main dir is:', '\n', src_dir)


def main():
    ''' with docopt main() expects a dictionary with arguments from docopt()
    docopt will automatically check your docstrings for usage, set -h, etc.
    '''
    options = docopt.docopt(__doc__)
    print(str( '\n' + 'Done, py35_packaging test worked!'
             ))
    return

def doSuperTest():
    print('doSuperTest works')
    return

if __name__ == '__main__':
    # if using docopt:
    # it will check all arguments pass, if not exits with 'Usage
    # if arguments are valid, run the program:
    sys.exit(main())
