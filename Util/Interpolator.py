#!usr/bin/env python
__author__ = "Hosik Chae"
__email__ = "CKMagenta@gmail.com"
__copyright__ = "Copyright 2016 RoMeLa"

import numpy as np
import pylab

class TrapezoidalVelocity(object):
    STATE_PAUSE = 0
    STATE_PRE   = 1
    STATE_MAIN  = 2
    STATE_POST  = 3

    ZERO_FINISH = 0
    NON_ZERO_FINISH = 1

    def __init__(self, dt, w_max, n_joints, q0, w0, T_acc_default=0.2):
        self.dt = dt
        self.t_zero = max(0.1 * dt, np.finfo(float).eps)
        self.w_max = w_max
        self.n_joints = n_joints

        ### Vectors of size (joint.NUMBER_OF_JOINTS,)
        # position vectors
        self.q0 = np.zeros((self.n_joints,))  # when entering this STEP or current position
        self.q1 = np.zeros((self.n_joints,))  # when  leaving this STEP or target position in this STEP
        self.dq = np.zeros((self.n_joints,))  # the position difference that should be moved in this STEP, i.e. = q1-q0
        self.q  = np.zeros((self.n_joints,))  # (return value) interpolated, current target position at the current time t

        # velocity vectors of size (joint.NUMBER_OF_JOINTS,)
        self.w0     = np.zeros((self.n_joints,))  # when entering PRE_SECTION of this STEP
        self.w_main = np.zeros((self.n_joints,))  # when  leaving PRE_SECTION or max velocity in this STEP
        self.w      = np.zeros((self.n_joints,))  # (return value) interpolated, current target velocity at the current time t

        # The amount of position to move in each SECTION
        self.dq_main = np.zeros((self.n_joints,))  # in MAIN_SECTION
        self.dq_pre  = np.zeros((self.n_joints,))  # in PRE_SECTION
        self.dq_post = np.zeros((self.n_joints,))  # in POST_SECTION

        # Flag that indicates if each motors should have zero velocity or not when this STEP finish
        self.finish_type = np.zeros((self.n_joints,))  # ZERO_FINISH or NON_ZERO_FINISH

        ### MotionPlanner's internal state
        self._just_started = np.array([True,]*self.n_joints, dtype=bool)
        self._finished = np.array([False,]*self.n_joints, dtype=bool)
        self.section = np.array([None,]*self.n_joints)  # Indicator which section currently it is in

        ### Time step values ###
        self.T_acc_default = T_acc_default
        self.t_global = 0.0
        self.T  = 0.0 * np.ones((self.n_joints,), dtype=float)  # Total time for this STEP (including PRE/MAIN/POST/PAUSE_SECTION)
        self.t  = 0.0 * np.ones((self.n_joints,), dtype=float)  # Current time in each section
        self.nT = 0   * np.ones((self.n_joints,), dtype=int)

        self.T_section = 0.0 * np.ones((self.n_joints,), dtype=float)
        self.T_pre     = 0.0 * np.ones((self.n_joints,), dtype=float)
        self.T_main    = 0.0 * np.ones((self.n_joints,), dtype=float)
        self.T_post    = 0.0 * np.ones((self.n_joints,), dtype=float)

        if q0 != None:
            self.q0 = q0

        if w0 != None:
            self.w0 = w0


    def set_next_point(self, T, q1, q2=None):

        q1 = [q1v if q1[idx] != None else self.q1[idx] for idx, q1v in enumerate(q1)]
        self.q1 = q1
        input_idx = [idx for idx, val in enumerate(q1) if val != None]

        if q2 == None:
            q2 = np.array([None, ] * self.n_joints)

        self.q2 = [q2v if q2[idx] != None else self.q1[idx] for idx, q2v in enumerate(q2)]

        self.T[input_idx] = T
        self._just_started[input_idx] = [True, ] * len(input_idx)
        self.section[input_idx] = TrapezoidalVelocity.STATE_PAUSE

        for idx in input_idx:
            self._process(idx)

    def take_step(self):

        just_started_idx = (self._just_started == True)

        self._just_started[just_started_idx] = False
        self._finished[just_started_idx] = False

        self.t[just_started_idx] = 0.0

        self.q[just_started_idx] = self.q0[just_started_idx]
        self.w[just_started_idx] = self.w0[just_started_idx]

        self.section[just_started_idx] = [TrapezoidalVelocity.STATE_PRE,]*just_started_idx.sum()
        self.T_section[just_started_idx] = self.T_pre[just_started_idx]
        #print "PAUSE -> PRE"


        idx = [idx for idx, val in enumerate(self.q1) if val != None]
        self.t_global += self.dt
        self.t[idx] += self.dt

        section_string = ["PAUSE", "PRE", "MAIN", "POST"]
        #print self._t, "\t", section_string[self.section], "\t", self.t, "\t / ", self.T_section

        dq = np.zeros((self.n_joints,))
        dw = np.zeros((self.n_joints,))

        idx_PRE = (self.section == TrapezoidalVelocity.STATE_PRE)
        dw[idx_PRE] = (self.w_main[idx_PRE] - self.w0[idx_PRE]) * self.t[idx_PRE] / self.T_section[idx_PRE]
        self.w[idx_PRE] = self.w0[idx_PRE] + dw[idx_PRE]
        dq[idx_PRE] = 0.5 * (self.w0[idx_PRE] + self.w[idx_PRE]) * self.t[idx_PRE]
        self.q[idx_PRE] = self.q0[idx_PRE] + dq[idx_PRE]

        idx_MAIN = (self.section == TrapezoidalVelocity.STATE_MAIN)
        dq[idx_MAIN] = self.dq_main[idx_MAIN] * self.t[idx_MAIN] / self.T_main[idx_MAIN]
        self.q[idx_MAIN] = self.q0[idx_MAIN] + dq[idx_MAIN]
        self.w[idx_MAIN] = self.w_main[idx_MAIN]

        idx_POST = (self.section == TrapezoidalVelocity.STATE_POST)
        dw[idx_POST] = np.where(self.finish_type[idx_POST] == TrapezoidalVelocity.NON_ZERO_FINISH,
                      0, -self.w_main[idx_POST] * self.t[idx_POST] / self.T_post[idx_POST]) # TODO
        self.w[idx_POST] = self.w0[idx_POST] + dw[idx_POST]
        dq[idx_POST] = 0.5 * (self.w_main[idx_POST] + self.w[idx_POST]) * self.t[idx_POST]
        self.q[idx_POST] = self.q0[idx_POST] + dq[idx_POST]

        """
        if self.section == TrapezoidalVelocity.STATE_PRE:
            dw = (self.w_main - self.w0) * self.t / self.T_section
            self.w = self.w0 + dw

            dq = 0.5 * (self.w0 + self.w) * self.t
            self.q = self.q0 + dq

        elif self.section == TrapezoidalVelocity.STATE_MAIN:
            dq = self.dq_main * self.t / self.T_main
            self.q = self.q0 + dq

            self.w = self.w_main

        elif self.section == TrapezoidalVelocity.STATE_POST:
            # dw = -self.w_main * self.t / self.T_post
            # self.w = self.w0 + dw
            # dq = 0.5 * (self.w_main + self.w) * self.t
            # self.q = self.q0 + dq

            dw = np.where(self.finish_type == TrapezoidalVelocity.NON_ZERO_FINISH,
                          0, -self.w_main * self.t / self.T_post)
            self.w = self.w0 + dw

            dq = 0.5 * (self.w_main + self.w) * self.t
            self.q = self.q0 + dq

        else: # TrapezoidalVelocity.PAUSE_SECTION
            pass
        """

        #print self.t_global, "\t", section_string[self.section], "\t", self.q, "\t", self.w, "\t", self._finished

        advance_idx = [idx for idx, val in enumerate(self.t) if self.t[idx] >= self.T_section[idx] - self.t_zero]
        self.advance_section(advance_idx)
        """
        if self.t >= self.T_section - self.t_zero:
            print "advance"
            self.advance_section()
        else: #self.t < self.T_section
            pass
        """
        return (self.q, self.w, self._finished)

    def advance_section(self, idx):

        self.t[idx] = 0.0

        # Re-assign starting point
        self.q0[idx] = self.q[idx]
        self.w0[idx] = self.w[idx]

        # Section Advance
        for iidx in idx:
            if self.section[iidx] == TrapezoidalVelocity.STATE_PRE:
                self.section[iidx] = TrapezoidalVelocity.STATE_MAIN
                self.T_section[iidx] = self.T_main[iidx]
                # print "PRE -> MAIN"

            elif self.section[iidx] == TrapezoidalVelocity.STATE_MAIN:
                self.section[iidx] = TrapezoidalVelocity.STATE_POST
                self.T_section[iidx] = self.T_post[iidx]
                # print "MAIN -> POST"

            elif self.section[iidx] == TrapezoidalVelocity.STATE_POST:
                self.section[iidx] = TrapezoidalVelocity.STATE_PAUSE
                self.T_section[iidx] = None
                self._finished[iidx] = True
                # print "POST -> PAUSE"
            else:
                self.section[iidx] = TrapezoidalVelocity.STATE_PRE
                self.T_section[iidx] = self.T_pre[iidx]
                # print "PAUSE -> PRE"

    def _process(self, idx):
        T = self.T[idx]
        dq = self.q1[idx] - self.q0[idx]
        nT = np.ceil(T / self.dt).astype(int)

        # Find direction change
        direction_changes = ((self.q0[idx] - self.q1[idx]) * (self.q2[idx] - self.q1[idx]) >= 0)

        # Find finish type
        if direction_changes == True:
            finish_type = TrapezoidalVelocity.ZERO_FINISH
        else:
            finish_type = TrapezoidalVelocity.NON_ZERO_FINISH

        if T > (2 * self.T_acc_default):
            T_pre  = self.T_acc_default
            T_post = self.T_acc_default
        elif nT >= 3:
            # adjusting acceleration time to be as long as possible
            T_pre  = (T // 2) * self.dt
            T_post = (T // 2) * self.dt
        else:
            # infeasible, ignoring the STEP
            T = 0.0
            T_pre = 0.0
            T_post = 0.0

        T_main = T - T_pre - T_post


        w_num = 0.0
        w_den = 1.0 # different when Zero_Finish /-\, Non-Zero Finish /-


        if self.w0[idx] * dq >= 0:  # TODO : if w0 ~= 0, dq == 0
            w_num = dq - 0.5 * T_pre * self.w0[idx]

            if finish_type == TrapezoidalVelocity.NON_ZERO_FINISH:
                w_den = 0.5 * T_pre + (T_main + T_post)
            else:
                w_den = 0.5 * T_pre + T_main + 0.5 * T_post  #


        else:   # w0 * w_main < 0
            w_num = dq - 0.5 * T_pre * self.w0[idx]
            if finish_type == TrapezoidalVelocity.NON_ZERO_FINISH:
                w_den = T_pre + (T_main + T_post)
            else:
                w_den = 0.5* T_pre + T_main + 0.5 * T_post

        # TODO: when T_pre, T_post == 0 -> den == 0
        w_main = np.clip( w_num / w_den, -self.w_max, self.w_max)

        if self.w0[idx] * dq >= 0:
            dq_pre = 0.5 * T_pre * (self.w0[idx] + w_main)
        else: # w0 * w_main < 0
            dq_pre = T_pre * (self.w0[idx] + w_main)

        if finish_type == TrapezoidalVelocity.NON_ZERO_FINISH:
            dq_main = dq-dq_pre   # NON_ZERO : No decel
            dq_post = 0.0
        else:
            dq_main = T_main * w_main
            dq_post = 0.5 * T_post * w_main

        self.dq[idx] = dq
        self.dq_pre[idx] = dq_pre
        self.dq_main[idx] = dq_main
        self.dq_post[idx] = dq_post

        self.w_main[idx] = w_main
        self.nT[idx] = nT
        self.T[idx] = T
        self.T_pre[idx] = T_pre
        self.T_post[idx] = T_post
        self.T_main[idx] = T_main
        self.finish_type[idx] = finish_type


        print "T=", self.T, ", T_pre=", self.T_pre, ", T_main=", self.T_main, ", T_post=", self.T_post
        print "dq=", self.dq, ", dq_pre=", self.dq_main, ", dq_main=", self.dq_main, ", dq_post=", self.dq_post, ", w_main=", self.w_main


    def _is_position_valid(self, data):
        return np.ones((self.n_joints,))

    @staticmethod
    def plot_result(result):

        ### plotting position values
        f, axarr = pylab.subplots(4, 5)
        for pi in range(0, 20):
            axarr[pi / 5][pi % 5].plot(result['time'], result['given_value'][pi, :], 'r*', markersize=15.0)
            axarr[pi / 5][pi % 5].plot(0.008 * np.array(range(0, result['count'])), result['value'][pi, :])  # 0.00025
            axarr[pi / 5][pi % 5].set_ylim([0, 4095])

        ### plotting velocity values
        f2, axarr2 = pylab.subplots(4, 5)
        for pi in range(0, 20):
            axarr2[pi / 5][pi % 5].plot(0.008 * np.array(range(0, resNult['count'])), result['velocity'][pi, :])  # 0.00025
            axarr2[pi / 5][pi % 5].set_ylim([-1200, 1200])

        pylab.show()

    @classmethod
    def test():
        idx = 0

        dt = 0.05
        w_max = 2000
        n_joints = 2
        q0 = np.array([40.0, -10.0])
        w0 = np.array([-10.0, -20.0])
        T_acc_default = 0.25
        q_target = np.array([[100, None, 50], [50, 100, 200]])
        T = [1, 1, 1]

        intp = TrapezoidalVelocity(dt, w_max, n_joints, q0, w0, T_acc_default=T_acc_default)
        intp.set_next_point(T[0], q_target[:, 0], q_target[:, 1])

        qs = np.multiply(intp.q0.reshape(len(intp.q0), 1), np.ones((intp.n_joints, 1)))
        ws = np.multiply(intp.w0.reshape(len(intp.w0), 1), np.ones((intp.n_joints, 1)))
        while True:
            (q, w, finished) = intp.take_step()
            # print q, w, finished

            qs = np.append(qs, q.reshape((intp.n_joints, 1)), axis=1)
            ws = np.append(ws, w.reshape((intp.n_joints, 1)), axis=1)

            if finished.any() == True:
                idx += 1
                if idx > q_target.shape[1] - 1:
                    break
                next_point = q_target[:, idx]
                next_next_point = q_target[:, idx] if idx == q_target.shape[1] - 1 else None
                intp.set_next_point(T[idx], next_point, next_next_point)

        ts = np.arange(0, intp.T[0] * idx + intp.dt, intp.dt)
        plot_time = ts.reshape((len(ts), 1))
        plot_qdata = qs.transpose()
        plot_wdata = ws.transpose()
        # print plot_time.shape, plot_data.shape


        ax1 = pylab.subplot(211)
        for qdx in range(plot_qdata.shape[1]):
            pylab.plot(plot_time, plot_qdata[:, qdx], label="joint %i" % qdx)

        t_bound = ax1.get_xlim()
        q_bound = ax1.get_ylim()
        for qt in q_target.flatten():
            pylab.plot(t_bound, [qt, qt], 'k--')
        for tdx in range(len(T)):
            pylab.plot([tdx, tdx], q_bound, 'k--')

        ax2 = pylab.subplot(212)
        pylab.plot(plot_time, plot_wdata)
        w_bound = ax2.get_ylim()
        for tdx in range(len(T)):
            pylab.plot([tdx, tdx], w_bound, 'k--')

        pylab.show()
        print "done!"