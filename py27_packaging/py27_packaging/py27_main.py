'''
Wrapper to execute main script function
'''
from __future__ import print_function
from __future__ import unicode_literals
from __future__ import division
from __future__ import absolute_import

from future import standard_library
standard_library.install_aliases()
import sys
import os

import py27_packaging.py27_packaging 

_ROOT = os.path.abspath(os.path.dirname(__file__))

print(_ROOT)

sys.path.append(_ROOT)

try:
    print(dir(py27_packaging.py27_packaging))
    print(dir(py27_packaging))
    for i in dir(py27_packaging):
        print(i)
        print(py27_packaging.i)
 
    print(py27_packaging.__loader__)
    py27_packaging.__loader__

    print(py27_packaging.load_entry_point)
    py27_packaging.load_entry_point

except AttributeError:
    print('AttributeError')

except NameError:
    print('NameError')


print('Hello, this works')


#py27_packaging.doSuperTest()

try:
    doSuperTest()

except NameError:
    print('NameError')

py27_packaging.py27_packaging.main()

