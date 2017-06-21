'''
Utilities for project_quickstart.py

:Author: Antonio J Berlanga-Taylor
:Date:


Boilerplate tools for quickstarting a data analysis project:

https://github.com/AntonioJBT/project_quickstart

'''
#################
import os
import sys

# Set up calling parameters from INI file:
# Modules with Py2 to 3 conflicts
try:
    import configparser
except ImportError:  # Py2 to Py3
    import ConfigParser as configparser
# Global variable for configuration file ('.ini'):
CONFIG = configparser.ConfigParser(allow_no_value = True)
#################


################# 
_ROOT = os.path.abspath(os.path.dirname(__file__))
def getINIdir(path = _ROOT):
    ''' Search for an INI file, default is where this script is executed from '''
    f_count = 0
    for f in os.listdir(path):
        if (f.endswith('.ini') and not f.startswith('tox')):
            f_count += 1
            INI_file = f
    if (f_count > 1 or f_count == 0):
        print('You have no "xxx.ini" or more than one "xxx.ini" file ',
                'in the directory:', '\n', path)
        sys.exit()

    INI_file_dir = os.path.abspath(os.path.join(path, INI_file))

    return(INI_file_dir)
################# 


#################
# Get locations of source code
    # os.path.join note: a subsequent argument with an '/' discards anything
    # before it
    # For function to search path see: 
    # http://stackoverflow.com/questions/4519127/setuptools-package-data-folder-location
# MANIFEST.in file instructs the project_quickstart/templates folder to be included in installation

_ROOT = os.path.abspath(os.path.dirname(__file__))
def getDir(path = _ROOT):
    ''' Get the absolute path to where this function resides. Useful for
    determining the user's path to a package. If a sub-directory is given it
    will be added to the path returned. Use '..' to go up directory levels. '''
   # src_top_dir = os.path.abspath(os.path.join(_ROOT, '..'))
    src_dir = _ROOT
    return(os.path.abspath(os.path.join(src_dir, path)))
#################
