'''
Wrapper to execute main script function
'''

import sys
import os

#from py35_packaging.py35_packaging import py35_packaging
#from py35_packaging import py35_packaging
import py35_packaging # works python -m py35_main and python py35_main.py

_ROOT = os.path.abspath(os.path.dirname(__file__))

sys.path.append(_ROOT)

try:
    print(dir(py35_packaging))

    print(py35_packaging.__loader__)
    py35_packaging.__loader__

    print(py35_packaging.load_entry_point)
    py35_packaging.load_entry_point

except AttributeError:
    print('AttributeError')

print('Hello, this works')

py35_packaging.doSuperTest()

py35_packaging.main()

