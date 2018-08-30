#!usr/bin/env python
__author__ = "Joshua Hooks"
__email__ = "hooksjrose@gmail.com"
__copyright__ = "Copyright 2016 RoMeLa"
__date__ = "June 12, 2017"

__version__ = "0.0.1"
__status__ = "Prototype"

"""
    Used to read from all motors with ids between 0-76.
"""

import time
import Library.dxl_manager.dcm_controller as dcmc
import sys
import select
import os

# Connect to motors and return all motor ids between 0-76
motor_ids = dcmc.main_import(True)

stop = False
N = 0.0
total_time = 0.0

deltas = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] 
previous = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] 
first = True

# Continously read motor positions until "enter" is pressed
while not stop:

    # read from all motors that are connected
    t0 = time.time()
    pos = dcmc.get_current_position(*motor_ids)
    current_time = time.time() - t0

    if first:
	previous = pos
	first = False

    # Calculate average read time
    total_time += current_time
    N += 1.0
    avg_time = total_time/N

    # print out joint positions
    os.system('clear')  # clear the terminal window before printing out data
    for i in range(0, len(motor_ids)):
        print "Joint %s: %s" % (motor_ids[i], round(pos[i],5))
	if abs(previous[i] - round(pos[i],5)) > deltas[i]:
		deltas[i] = abs(previous[i] - round(pos[i],5))
	previous[i] = round(pos[i],5)
    print
    print "Avg read time: %s" % (round(avg_time,5))
    print "Avg read in Hz: %s" % (round((1.0/avg_time),5))
    print
    print "Press enter to exit."

    time.sleep(0.1)

    # check if the user has pressed enter
    if sys.stdin in select.select([sys.stdin], [], [], 0)[0]:
        stop = True
	print deltas

