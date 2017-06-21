'''script_quickstart.py - setup a new python based project
===========================================================

:Author: Antonio Berlanga-Taylor
:Release:
:Date:

Purpose
-------
This script creates a
The main idea is to
We've additionally put
Once you've quickstarted your

For a pipeline quickstart based on a Ruffus and CGAT framework see also:
https://github.com/CGATOxford/CGATPipelines/blob/master/scripts/pipeline_quickstart.py
(on which this code is based on)

Usage and Options
=================
Quickstart a script
.. These are using docopt: http://docopt.org/
.. https://github.com/docopt/docopt
.. An example for loading arguments from an INI file: https://github.com/docopt/docopt/blob/master/examples/config_file_example.py

Usage:
    script_quickstart.py (--script_name=<script_name>) ...

to create a script template.

    script_quickstart.py (--script_name | -n) <script_name>
    script_quickstart.py (--language | -lang) <python | R>
    script_quickstart.py [-f | --force]
    script_quickstart.py -h | --help
    script_quickstart.py --version
    script_quickstart.py --quiet
    script_quickstart.py --verbose
    script_quickstart.py [-L | --log] <project_quickstart.log>

Options:
    -lang --language    R or Python templates available
    -f --force   Take care, forces to overwrite files and directories.
    -h --help    Show this screen.
    --version    Show version.
    --quiet      Print less text.
    --verbose    Print more text.
    -L --log     Log file name. [default: project_quickstart.log]
Documentation
-------------
.. todo::
  Add docs
  New string formatting https://pyformat.info/ '{} {}'.format('one', 'two')

Code
----
'''
from __future__ import print_function
from __future__ import unicode_literals
from __future__ import division
from __future__ import absolute_import
##############################
from future import standard_library
standard_library.install_aliases()
import sys
import re
import os
import shutil
import collections
#import CGAT.Experiment as E
from docopt import docopt

try:
    import configparser
except ImportError: # Py2 to Py3
    import configparser as configparser

# Check configuration and print to standard out
# See https://github.com/CGATOxford/CGATPipelines/blob/master/CGATPipelines/Pipeline/Parameters.py
# https://github.com/CGATOxford/cgat/blob/master/CGAT/Experiment.py

# Global variable for configuration file ('.ini'):
CONFIG = configparser.ConfigParser()


# TO DO: (see E.py)

CONFIG.read('project_quickstart.ini')
for key in CONFIG:
    print(key, CONFIG[key])
##############################


##############################
def main():
    ''' Create a Python or R script template
    '''
    # Copy files from template directory:
    def copyTemplate(source_dir, project_dir):
        '''
        Copy across template files and directories for a Python/GitHub/etc setup.
        TO DO: 'code' dir is hard coded, change to ini parameter later
        The intention is to use the 'code' dir as a GitHub/packageable directory
        '''
        copy_from = project_template
        copy_to = os.path.join(project_dir, 'code')

        if os.path.exists(copy_to) and not options['--force']:
            raise OSError(
                '''file {} already exists - not overwriting, see --help or use --force
                to overwrite.'''.format(project_name)
                         )

        shutil.copytree(copy_from, copy_to) # https://docs.python.org/3/library/shutil.html


    # Replace all instances of template with 'name' from project_'name' as
    # specified in options:
    def rename(project_dir, old_substring, project_name):
        ''' rename 'template' to 'project' from template file names '''
        for dirpath, dirname, filename in os.walk(project_dir):
            for filename in files:
                os.rename(os.path.join(project_dir, filename),
                          os.path.join(project_dir,
                              filename.replace('template',
                                  '{}'.format(project_name))))

    # Create links for the manuscript and lab_notebook 
    # templates to go into the 'manuscript' directory:
    for template_dir, project_dir in (("manuscript_template.rst", "lab_notebook_template.rst"),
                                      ("manuscript_template.rst", "lab_notebook_template.rst")):
        d = os.path.join("", project_dir)
        if os.path.exists(d) and options['--force']:
            os.unlink(d)
        os.symlink(os.path.join(project_dir, ), d)

    # Print a nice welcome message (if successful):
    print(""" Done, welcome to your {1}!

    Remember to

    Feel free to raise issues, fork or contribute at:

    https://github.com/

    Have fun!
    """.format(1, 2, 3, )
         )

if __name__ == "__main__":
   sys.exit(main())
