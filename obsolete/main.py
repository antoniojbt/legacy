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

import project_quickstart.project_quickstart
import project_quickstart.projectQuickstart


_ROOT = os.path.abspath(os.path.dirname(__file__))
print(_ROOT)
sys.path.append(_ROOT)


print(dir(project_quickstart.project_quickstart))


print('Hello, this works')

project_quickstart.project_quickstart.doSuperTest()

project_quickstart.project_quickstart.main()

