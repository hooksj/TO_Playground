#!/usr/bin/env python

__author__  	= "Joshua Hooks"
__email__   	= "hooksjrose@gmail.com"
__copyright__ 	= "Copyright 2017 RoMeLa"
__date__ 		= "Aug. 29th, 2018"

__version__ 	= "0.1.0"
__status__ 		= "Prototype"

"""
Useful high level classes to organize optimization problem. These classes are designed to be inherited by specific
constraints and costs.
"""


class Constraint(object):

    def __init__(self):

        self.ncon = 0
        self.nnjz = 0

    def evaluate(self):
        pass

    def derivative_structure(self):
        pass

    def derivative(self):
        pass

    def boundaries(self):
        pass


class Cost(object):

    def __init__(self):
        pass

    def evaluate(self):
        pass

    def derivative(self):
        pass
