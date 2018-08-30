#!usr/bin/env python
__author__ = "Min Sung Ahn"
__email__ = "aminsung@gmail.com"
__copyright__ = "Copyright 2016 RoMeLa"
__date__ = "January 1, 1999"

__version__ = "0.0.1"
__status__ = "Prototype"

import collections

"""
Only edit if you know what you are doing.
"""

"""
class NestedDict():
    def __getitem__(self, item):
        try:
            return dict.__getitem__(self, item)
        except KeyError:
            value = self[item] = type(self)()
            return value
#"""

def create():
    return collections.defaultdict(lambda: collections.defaultdict(int))