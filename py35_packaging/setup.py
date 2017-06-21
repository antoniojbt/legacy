'''
setup for |project_name|

For example on setting a Python package, see:
https://packaging.python.org/en/latest/distributing.html
https://github.com/pypa/sampleproject

For Python 3.5
Before packaging or installing run:

    pip install -U pip twine check-manifest setuptools

TO DO: to add tests see https://python-packaging.readthedocs.io/en/latest/testing.html

To package, do something like this:

    check-manifest
    python setup.py check
    python setup.py sdist bdist_wheels

which will create a dist/ directory and a compressed file inside with your package.

More notes and references in:
    https://github.com/EpiCompBio/welcome

And in the Python docs of course.

Upload to PyPI after this if for general use. Register yourself and package,
then:

Test first in https://testpypi.python.org/pypi

When ready to upload run:

    twine upload dist/*

Complete instructions in:

https://packaging.python.org/distributing/#uploading-your-project-to-pypi
'''

# Always prefer setuptools over distutils
from setuptools import setup, find_packages
# To use a consistent encoding
from codecs import open

import sys
import os

# Get location to this file:
here = os.path.abspath(os.path.dirname(__file__))
print(here)


# Get info from README and version.py:
# Get Ptyhon modules required:
install_requires = []

with open(os.path.join(here, 'requirements.rst'), encoding='utf-8') as required:
    for line in (required):
        if not line.startswith('#') and not line.startswith('\n'):
            line = line.strip()
            install_requires.append(line)

print(install_requires)

# Use README as long description if desired, otherwise get it from INI file (or
# write it out in setup()):

with open(os.path.join(here, 'README.rst'), encoding='utf-8') as readme:
    description = readme.read()

# Actual setup.py instructions
# Python docs: https://docs.python.org/3.6/distutils/setupscript.html 
# Tutorial: http://python-packaging.readthedocs.io/en/latest/
print(find_packages())

setup(
      # Package metadata:
      name = 'py35_packaging',
      version = '0.1.1',
      url = 'https://github.com/AntonioJBT/py35_packaging',
      download_url = 'https://github.com/AntonioJBT/py35_packaging.git',
      author = 'Antonio J Berlanga-Taylor',
      author_email = 'a.berlanga@imperial.ac.uk',
      license = 'GPL-3',
      description = 'Simple repository to test Python packaging',
      keywords = 'python3 packaging utility',
      long_description = description,
      # Package information:
      packages = find_packages(),
      install_requires = install_requires,
#      package_dir = {'' : 'py35_packaging'},
      scripts = ['py35_main.py.py'],
#      entry_points = {'console_scripts': [ 'py35_packaging.py = py35_packaging:main' ]},
      zip_safe = False,
          )
