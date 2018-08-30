#!usr/bin/env python
__author__ = "Hosik Chae"
__email__ = "CKMagenta@gmail.com"
__copyright__ = "Copyright 2017 RoMeLa"
__date__ = "10/30/2017"

import time
class SyncTimer(object):
    t0 = 0
    t1 = 0

    def __init__(self, dt):
        self.dt = dt

    def start(self, t0 = None):
        if t0 == None:
            self.t0 = time.time()
            self.t1 = self.t0
        else:
            self.t0 = t0
            self.t1 = t0

        return self.t0

    def set_dt(self, dt):
        self.dt = dt

    def set_next_timestep(self, t1):
        self.t1 = t1 + self.t0

    def to_next_timestep(self, n_step = 1):
        self.t1 += self.dt * n_step

    def get_time_elapsed(self):
        return time.time() - self.t0

    def wait_until_next_timestep(self, to_next_n_timestep=None):
        if time.time() < self.t1:
            return True
        else:
            if to_next_n_timestep != None:
                self.to_next_timestep(to_next_n_timestep)
            return False

    @classmethod
    def sleep(cls, t, dt, simulator_step_advance=None):
        t0 = time.time()
        if simulator_step_advance == None:
            while (time.time() - t0 < t + dt) :
                pass
        else:
            for idx in range(int(t//dt)):
                simulator_step_advance()