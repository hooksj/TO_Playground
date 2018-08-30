#!usr/bin/env python
__author__ = "Joshua Hooks"
__email__ = "hooksjrose@gmail.com"
__copyright__ = "Copyright 2016 RoMeLa"
__date__ = "June 12, 2017"

__version__ = "0.0.1"
__status__ = "Prototype"

"""
    Used to pose the ALPHRED robot.
"""

import time
import Library.dxl_manager.dcm_controller as dcmc
import sys
import select
import os

# Connect to motors and return all motor ids between 0-76
dcmc.main_import()

def controller():
    cmd = ""
    while cmd != "q":
        cmd = str(raw_input("Input Command: "))

        if cmd == "1 on":
            dcmc.set_torque_enable((1,1), (2,1), (3,1))
        elif cmd == "2 on":
            dcmc.set_torque_enable((4,1), (5,1), (6,1))
        elif cmd == "3 on":
            dcmc.set_torque_enable((7,1), (8,1), (9,1))
        elif cmd == "4 on":
            dcmc.set_torque_enable((10,1), (11,1), (12,1))
        elif cmd == "1 off":
            dcmc.set_torque_enable((1,0), (2,0), (3,0))
        elif cmd == "2 off":
            dcmc.set_torque_enable((4,0), (5,0), (6,0))
        elif cmd == "3 off":
            dcmc.set_torque_enable((7,0), (8,0), (9,0))
        elif cmd == "4 off":
            dcmc.set_torque_enable((10,0), (11,0), (12,0))
        elif cmd == "q":
            break
        else:
            print ""


controller()
