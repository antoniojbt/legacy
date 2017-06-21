#################
# Load arguments for docopt from an INI file 
# Modified from:
# https://github.com/docopt/docopt/blob/master/examples/config_file_example.py
# See also (which includes reading from json):
# http://stackoverflow.com/questions/8884188/how-to-read-and-write-ini-file-with-python
# https://wiki.python.org/moin/ConfigParserExamples
# CGAT tools do a better job, using argparser.

def load_ini_config():
    ''' Loads an *.ini file if present in the current directory and prepares it
    for use with docopt'''

    try:  # Python 2
        from ConfigParser import ConfigParser
        from StringIO import StringIO
    except ImportError:  # Python 3
        from configparser import ConfigParser
        from io import StringIO

    import sys
    import os
    import glob

    # By using `allow_no_value=True` we are allowed to
    # write `--force` instead of `--force=true` below.
    config = ConfigParser(allow_no_value=True)

    # Check if there's an INI file present:
    if len(glob.glob('*.ini')) == 0:
        return 'No configuration file (e.g. xxx.ini) present in the directory'
    elif len(glob.glob('*.ini')) > 1:
        return ''' More than one configuration file present
        (several files ending in ".ini" files") present in the directory')
        '''
    else:
        INI_file = glob.glob('*.ini')

    # ConfigParser requires a file-like object and
    # no leading whitespace.
    with open(INI_file, 'r') as f:
        config_file = f.StringIO('\n'.join(INI_file.split()))
        config_file = config.read_file(config_file)
        config_file = config.read(config_file)

        for section, key, value in config_file:
            print(section, key, value)

    # ConfigParsers sets keys which have no value
    # (like `--force`) to `None`. Thus we
    # need to substitute all `None` with `True` for docopt:
    return dict((key, True if value is None else value)
                for key, value in config_file)


    ini_config = load_ini_config()

    # Arguments take priority over INI, INI takes priority over JSON:
    result = arguments, ini_config

    from pprint import pprint

    print('\nINI config:')
    pprint(ini_config)
    print('\nResult:')

    pprint(result)

    return(result)
#################
