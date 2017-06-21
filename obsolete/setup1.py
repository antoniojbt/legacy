#################

# https://python-packaging.readthedocs.io/en/latest/minimal.html
# For a fuller example see: https://github.com/CGATOxford/UMI-tools/blob/master/setup.py
# Or: https://github.com/CGATOxford/cgat/blob/master/setup.py

# TO DO: update with further options such as include README.rst and others when ready

# TO DO: to add tests see https://python-packaging.readthedocs.io/en/latest/testing.html

# See also this example: https://github.com/pypa/sampleproject/blob/master/setup.py

# This may be a better way, based on Py3: http://www.diveintopython3.net/packaging.html

# To package, check setup.py first:
# python setup.py check
# Then create a source distribution:
# python setup.py sdist
# which will create a dist/ directory and a compressed file inside with your package.
# For uploading to PyPi see http://www.diveintopython3.net/packaging.html#pypi

#################
import sys
import os

try:
    from setuptools import setup, find_packages
except ImportError:
    from ez_setup import use_setuptools
    use_setuptools("10.0")
    from setuptools import setup, find_packages

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

CONFIG.read('project_quickstart.ini')
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
#################


#################
# Get info from README and version.py:
# Get Ptyhon modules required:
install_requires = []

with open('requirements.rst') as required:
    for line in (required):
        if not line.startswith('#') and not line.startswith('\n'):
            line = line.strip()
            install_requires.append(line)

#print(install_requires)

# Use README as long description if desired, otherwise get it from INI file (or
# write it out in setup()):

#with open('README.rst', 'rt') as readme:
#    description = readme.read()


# Give warning:
class CustomInstall(install):
    def initialize_options(self):
        if sys.version < '3.6':
            print('Error during installation: ', '\n',
                    CONFIG['metadata']['project_name'],
                    ' requires Python 3.6 or higher.',
                    'Exiting...')
            sys.exit(1)

        return install.initialize_options(self)

entry_points = str('{"console_scripts"'
                  + ' : '
                  + ' [ '
                  + '"'
                  + CONFIG['metadata']['project_name']
                  + '.py'
                  + ' = '
                  + CONFIG['metadata']['project_name']
                  + '.'
                  + CONFIG['metadata']['project_name']
                  + ' : main'
                  + '" '
                  + ']}'
                  )
print(entry_points)
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
      long_description = CONFIG['metadata']['long_description'], #long_description = description,
      # Package information:
#      packages = find_packages(CONFIG['metadata']['project_name']),
      packages = find_packages(),
      #[CONFIG['metadata']['packages_setup']], # needs to be passed
                                                         # as list
      install_requires = install_requires,
#      include_package_data = True,
      package_dir = {CONFIG['metadata']['project_name']: CONFIG['metadata']['project_name']},
#      entry_points = entry_points,
       entry_points = {'console_scripts': [
           'project_quickstart.py = project_quickstart.project_quickstart:main'
           ]},
      cmdclass = {'install': CustomInstall},
      zip_safe = False,
      test_suite = "tests"
          )
#################
