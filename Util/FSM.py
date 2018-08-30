#!usr/bin/env python
__author__ = "Min Sung Ahn"
__email__ = "aminsung@gmail.com"
__copyright__ = "Copyright 2016 RoMeLa"
__date__ = "July 1, 2016"

__version__ = "0.1.0"
__status__ = "Production"

"""
FSM class for individual FSM's
"""

import re
import pdb

from collections import deque

class FSM(object):
    def __init__(self, state1, state_scripts, name=None):
        self.sm_name = name
        self.tag = "[FSM | " + self.sm_name + "]  "
        self.cs = state1
        self.state_scripts = state_scripts
        self.state_list = {}
        self.state_list[self.cs] = getattr(self.state_scripts, self.cs)
        self.triggers = {}
        self.triggers[self.cs] = {}
        self.current_state = self.state_list[self.cs]
        self.next_state = None

        self.trigger = deque() # Note it is different to self.triggers. You can trigger transitions through this attribute.

        self.state_stack = []

    def add_state(self, new_state):
        self.state_list[new_state] = getattr(self.state_scripts, new_state)
        self.triggers[new_state] = {}

    def set_transition(self, from_state, trigger, to_state):
        if (not from_state in self.triggers):
            self.triggers[from_state] = {}
        self.triggers[from_state][trigger] = to_state

    def set_state(self, state_to_set):
        """
        Forces a certain state to be the current state. Good for getch!
        """
        self.cs = state_to_set
        self.current_state = self.state_list[self.cs]
        self.entry()

    def _parse_state_name(self, state_name):
        '''Function to turn <this_type_of_name> into <this type of name>'''
        split_names = re.split('_', state_name)[1:]
        return " ".join(split_names).title()

    def start(self):
        while self.cs != "exit":
            # Call current state entry function
            print self.tag + "Entering state: {}".format(self._parse_state_name(self.cs))
            self.current_state.entry()

            # Call current state update function
            print self.tag + "Updating state: {}".format(self._parse_state_name(self.cs))
            self.trigger.append(self.current_state.update())

            # Based on the returned trigger, the next step is decided
            self.next_state = self.triggers[self.cs][self.trigger.popleft()]

            # TODO: Create an exception for non-defined triggers
            while self.next_state == self.cs:
                self.trigger.append(self.current_state.update())
                self.next_state = self.triggers[self.cs][self.trigger.popleft()]

            self.trigger.clear()
            print(self.tag + "A next/different step has been triggered!")

            print self.tag + "Exiting state: {}".format(self._parse_state_name(self.cs))
            self.current_state.exit()
            self.cs = self.next_state
            self.current_state = getattr(self.state_scripts, self.next_state)

        print self.tag+"Shutting down."

    def go_to(self, val_trigger):
        self.trigger.append(val_trigger)
        print("Trigger set to: {}".format(self.trigger))
