#!usr/bin/env python
__author__ = "Min Sung Ahn"
__email__ = "aminsung@gmail.com"
__copyright__ = "Copyright 2017 RoMeLa"
__date__ = "November 5, 2017"

__version__ = "0.0.1"
__status__ = "Prototype"

'''

'''

import pdb
import sys
import tty
import termios

def getch():
    fd = sys.stdin.fileno()
    old = termios.tcgetattr(fd)
    try:
        tty.setraw(fd)
        return sys.stdin.read(1)
    finally:
        termios.tcsetattr(fd, termios.TCSADRAIN, old)

if __name__ == '__main__':
    while True:
        cmd = getch()
        print ("The character you typed is: {}".format(cmd))

    print("Done!")