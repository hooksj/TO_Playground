#!usr/bin/env python
__author__ = "Joshua Hooks"
__email__ = "hooksjrose@gmail.com"
__copyright__ = "Copyright 2018 RoMeLa"
__date__ = "August 9th, 2018"

__version__ = "0.0.1"
__status__ = "Prototype"

"""
Sets up a UDP server/client socket in order to send commands to the robot.
"""

import socket
if __name__ == '__main__':
    from pynput.keyboard import Controller, Listener
import pdb

frequency = 750.0

class server():

    def __init__(self):
        self.udp_id = "192.168.123.14"
        #self.udp_id = "127.0.0.1"
        self.udp_port = 4456
        self.fd = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
        self.fd.bind((self.udp_id, self.udp_port))
        self.fd.settimeout(0.01)
        self.client_addresses = []
        self.connection = False
        self.keyboard = Controller()
        self.key = None

    def monitor(self):
        listener = Listener(on_press=self.on_press, on_release=self.on_release)
        listener.start()

        while True:
            try:
                message, address = self.fd.recvfrom(1000) #Look for the FSM trying to connect.
                print "Connected, ready to send commands"
                print
                if message == 'connect':
                    self.connection = True
                    self.client_addresses.append(address)
            except socket.timeout:
                pass

            if self.connection:
                if self.key is not None:
                    for address in self.client_addresses:
                        self.fd.sendto(self.key, address)

    # Callback function for listener thread
    def on_press(self, key):
        self.key = format(key)

    # Call back function for listener thread
    def on_release(self, key):
        self.key = None


class client():

    def __init__(self):
        self.udp_id = "192.168.123.14"
        #self.udp_id = "127.0.0.1"
        self.udp_port = 4456
        self.fd = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
        self.fd.settimeout(0.01)
        self.fd.sendto('connect', (self.udp_id, self.udp_port))

    def get_input(self):
        try:
            key = self.fd.recvfrom(1000)
            if key is not None:
                if len(key[0]) == 4:
                    return key[0][2]
        except socket.timeout:
            pass


    def disconnect(self):
        self.fd.sendto('disconnect', (self.udp_id, self.udp_port))



if __name__ == '__main__':
    s = server()
    s.monitor()
