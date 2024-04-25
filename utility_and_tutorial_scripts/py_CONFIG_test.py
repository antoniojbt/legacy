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

# Read INI file:
CONFIG.read('project_quickstart.ini')

# Print kyes (sections):
for key in CONFIG:
    for value in CONFIG[key]:
        print(key, value, '\n', CONFIG[key][value])

CONFIG['metadata']['project_name']
#print(CONFIG['metadata']['project_name'])
my_val = CONFIG['metadata']['project_name']
#my_val
#print(my_val)
print([CONFIG['metadata']['project_name']])
print([CONFIG['metadata']['packages_setup']])
