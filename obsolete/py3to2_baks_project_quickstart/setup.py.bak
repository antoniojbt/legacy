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

And in the Python docs.
Upload to PyPI after this if for general use.
'''

# Always prefer setuptools over distutils
from setuptools import setup, find_packages
# To use a consistent encoding
from codecs import open

# To run custom install warning use:
from setuptools.command.install import install

import sys
import os

# Get location to this file:
here = os.path.abspath(os.path.dirname(__file__))
print(here)

# Set up calling parameters from INI file:
# Modules with Py2 to 3 conflicts
try:
    import configparser
except ImportError:  # Py2 to Py3
    import ConfigParser as configparser

# Global variable for configuration file ('.ini')
# allow_no_value addition is from:
# https://github.com/docopt/docopt/blob/master/examples/config_file_example.py
# By using `allow_no_value=True` we are allowed to
# write `--force` instead of `--force=true` below.
CONFIG = configparser.ConfigParser(allow_no_value = True)

CONFIG.read(os.path.join(here, str('project_quickstart' + '.ini')))
# str(CONFIG['metadata']['project_name'] + '.ini'))) 

# Print keys (sections):
print('Values for setup.py:', '\n')
for key in CONFIG:
    for value in CONFIG[key]:
        print(key, value, CONFIG[key][value])
#################


#################
# Get version:
sys.path.insert(0, CONFIG['metadata']['project_name'])
import version

version = version.__version__
print(version)
#################


#################
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

# Give warning:
class CustomInstall(install):
    def initialize_options(self):
        if sys.version < '3.5':
            print('Error during installation: ', '\n',
                    CONFIG['metadata']['project_name'],
                    ' requires Python 3.5 or higher.',
                    'Exiting...')
            sys.exit(1)

        return install.initialize_options(self)

# See
# http://setuptools.readthedocs.io/en/latest/setuptools.html#dynamic-discovery-of-services-and-plugins
# For rst INI entry points
#################


#################
# Actual setup.py instructions
# Python docs: https://docs.python.org/3.6/distutils/setupscript.html 
# Tutorial: http://python-packaging.readthedocs.io/en/latest/
setup(
      # Package metadata:
      name = CONFIG['metadata']['project_name'],
      version = CONFIG['metadata']['version'],
      url = CONFIG['metadata']['project_url'],
      download_url = CONFIG['metadata']['download_url'],
      author = CONFIG['metadata']['author_name'],
      author_email = CONFIG['metadata']['author_email'],
      license = CONFIG['metadata']['license'],
#      classifiers = [CONFIG['metadata']['classifiers_setup']], 
       # needs to be passed as list
       # gives many errors when registering manually in pip
      description = CONFIG['metadata']['project_short_description'],
      keywords = CONFIG['metadata']['keywords'],
      #long_description = CONFIG['metadata']['long_description'],
      long_description = description,
      # Package information:
      packages = find_packages(),
      install_requires = install_requires,
      # If there are data files to include with installation, specify here
      # (they should be in the main src dir):
      # Including them in MANIFEST.in is not favoured.
#      package_data = {'sample': ['package_data.dat']},
      # 'package_data' is preferred but if data files are outside the package's
      # main dir then use:
      # http://docs.python.org/3.4/distutils/setupscript.html#installing-additional-files # noqa
      # 'data_file' will be installed into '<sys.prefix>/my_data'
#      data_files=[('my_data', ['data/data_file'])],
      #package_dir = {CONFIG['metadata']['project_name']: CONFIG['metadata']['project_name']},
      scripts = [str(CONFIG['metadata']['project_name'] + '/main.py')],
      #entry_points = {'console_scripts': ['project_quickstart.py = project_quickstart.project_quickstart:main']},
      cmdclass = {'install': CustomInstall},
      zip_safe = False,
      #test_suite = "tests"
          )
#################
