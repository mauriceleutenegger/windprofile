import sys

# for some reason, this behavior is different for python 2 and 3
if sys.version_info.major < 3 :
    from PyWindProfile import *
else:
    from PyWindProfile.PyWindProfile import *

