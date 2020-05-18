'''
extract_pubmed.py
=================

:Author: |author_name|
:Release: |version|
:Date: |today|


Purpose
=======

|description|


Usage and options
=================

These are based on docopt_, see examples_.

.. _docopt: https://github.com/docopt/docopt

.. _examples: https://github.com/docopt/docopt/blob/master/examples/options_example.py


Usage:
       script_name [--main-method]
       script_name [-I FILE]
       script_name [-O FILE]
       script_name [-h | --help]
       script_name [-V | --version]
       script_name [-f --force]
       script_name [-L | --log]

Options:
    -I             Input file name.
    -O             Output file name.
    -h --help      Show this screen
    -V --version   Show version
    -f --force     Force overwrite
    -L --log       Log file name.

Documentation
=============

    For more information see:

        |url|

'''
##############
# Get all the modules needed
# System:
import os
import sys
import glob

# Options and help:
from docopt import docopt

# required to make iteritems python2 and python3 compatible
from builtins import dict


# Basic function structure:
def my_func():
    '''
    Extract each element from a Pubmed email
    '''
    # Emails come from several sources: my email with links, Pubmed email, etc.
    # Get PMID, if not separate file
    # Title

    # Authors

    # Journal

    # Year

    # PMID


    # Make sure the setup is correct for this function:
    assert True == True

    # do something

    return


# Finish and exit with docopt arguments:
if __name__ == '__main__':
    arguments = docopt(__doc__, version='xxx 0.1')
    print(arguments)
    sys.exit(main())
