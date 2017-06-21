# Probably good to use MAJOR.MINOR.PATCH ?

#    MAJOR version (when you make incompatible API changes), or something arbitrarily import?
#    MINOR version when you add functionality in a backwards-compatible manner, and
#    PATCH version when you make backwards-compatible bug fixes.

# http://semver.org/


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

#for key in CONFIG:
#    for value in CONFIG[key]:
#        print(CONFIG[key][value])

__version__ = CONFIG['metadata']['version']

print(CONFIG['metadata']['project_name'], 'version', __version__)
