#!/usr/bin/env python
"""
A class used to solve for a trajectory using ipopt.
     - This serves as the structure needed to pass into ipopt
"""

__author__  	= "Joshua Hooks"
__email__   	= "hooksjrose@gmail.com"
__copyright__ 	= "Copyright 2017 RoMeLa"
__date__ 		= "Aug. 4th, 2018"

__version__ 	= "0.1.0"
__status__ 		= "Prototype"

"""
Step planner for teleoperated moves
"""

import numpy as np
import Util.MathFcn as MF
import pdb


def teleop_creep(start, initial_foot_pos, direction, orientation, step_dis, step_time, step_height, height, right):
    num_steps = 8
    dx = step_dis * np.cos(direction)
    dy = step_dis * np.sin(direction)
    delta = np.array([dx,dy,0.0])
    cm_x = 0.0
    cm_y = 0.0
    for i in range(4):
        cm_x += initial_foot_pos[i][0]/4.0
        cm_y += initial_foot_pos[i][1]/4.0
    CM = np.array([cm_x, cm_y, 0.0])

    if start:
        num_phases = 11
    else:
        num_phases = 10

    normal = True
    if orientation > 0.0:
        normal = False

    Tk = np.zeros(num_phases)
    T_H = np.zeros(num_phases)
    C1 = np.ones(num_phases)
    C2 = np.ones(num_phases)
    C3 = np.ones(num_phases)
    C4 = np.ones(num_phases)
    sh1 = np.zeros(num_phases)
    sh2 = np.zeros(num_phases)
    sh3 = np.zeros(num_phases)
    sh4 = np.zeros(num_phases)
    H = np.zeros(num_phases)
    F1 = np.zeros((num_phases,3))
    F2 = np.zeros((num_phases,3))
    F3 = np.zeros((num_phases,3))
    F4 = np.zeros((num_phases,3))

    if start:
        F1[0] = initial_foot_pos[0]
        F2[0] = initial_foot_pos[1]
        F3[0] = initial_foot_pos[2]
        F4[0] = initial_foot_pos[3]
        H[0] = height

        # double support
        Tk[1] = round(Tk[0] + 0.5,6)
        T_H[1] = round(T_H[0] + 0.5,6)
        F1[1] = initial_foot_pos[0]
        F2[1] = initial_foot_pos[1]
        F3[1] = initial_foot_pos[2]
        F4[1] = initial_foot_pos[3]
        H[1] = height
        index = 2
    else:
        #if right:
        #    C4[0] = 0.0
        #else:
        #    C3[0] = 0.0
        F1[0] = initial_foot_pos[0]
        F2[0] = initial_foot_pos[1]
        F3[0] = initial_foot_pos[2]
        F4[0] = initial_foot_pos[3]
        H[0] = height
        index = 1


    for i in range(2):
        CM += delta
        for j in range(4):
            Tk[index] = round(Tk[index-1] + step_time,6)
            T_H[index] = round(T_H[index-1] + step_time,6)
            F1[index] = F1[index-1]
            F2[index] = F2[index-1]
            F3[index] = F3[index-1]
            F4[index] = F4[index-1]
            H[index] = height
            if j == 0:
                if normal:
                    C1[index] = 0.0
                    sh1[index] = step_height
                    if i == 0:
                        f1_ni = F1[index] + delta
                        f1_nr = f1_ni - CM
                        f1_rr = np.dot(MF.yaw_rot_mat(orientation),f1_nr.reshape(3,1)).reshape(3,)
                        F1[index] = f1_rr + CM
                    else:
                        F1[index] += delta
                else:
                    C2[index] = 0.0
                    sh2[index] = step_height
                    if i == 0:
                        f2_ni = F2[index] + delta
                        f2_nr = f2_ni - CM
                        f2_rr = np.dot(MF.yaw_rot_mat(orientation),f2_nr.reshape(3,1)).reshape(3,)
                        F2[index] = f2_rr + CM
                    else:
                        F2[index] += delta
                index += 1
            elif j == 1:
                if normal:
                    C2[index] = 0.0
                    sh2[index] = step_height
                    if i == 0:
                        f2_ni = F2[index] + delta
                        f2_nr = f2_ni - CM
                        f2_rr = np.dot(MF.yaw_rot_mat(orientation),f2_nr.reshape(3,1)).reshape(3,)
                        F2[index] = f2_rr + CM
                    else:
                        F2[index] += delta
                else:
                    C1[index] = 0.0
                    sh1[index] = step_height
                    if i == 0:
                        f1_ni = F1[index] + delta
                        f1_nr = f1_ni - CM
                        f1_rr = np.dot(MF.yaw_rot_mat(orientation),f1_nr.reshape(3,1)).reshape(3,)
                        F1[index] = f1_rr + CM
                    else:
                        F1[index] += delta
                index += 1
            elif j == 2:
                if normal:
                    C3[index] = 0.0
                    sh3[index] = step_height
                    if i == 0:
                        f3_ni = F3[index] + delta
                        f3_nr = f3_ni - CM
                        f3_rr = np.dot(MF.yaw_rot_mat(orientation),f3_nr.reshape(3,1)).reshape(3,)
                        F3[index] = f3_rr + CM
                    else:
                        F3[index] += delta
                else:
                    C4[index] = 0.0
                    sh4[index] = step_height
                    if i == 0:
                        f4_ni = F4[index] + delta
                        f4_nr = f4_ni - CM
                        f4_rr = np.dot(MF.yaw_rot_mat(orientation),f4_nr.reshape(3,1)).reshape(3,)
                        F4[index] = f4_rr + CM
                    else:
                        F4[index] += delta
                index += 1
            elif j == 3:
                if normal:
                    C4[index] = 0.0
                    sh4[index] = step_height
                    if i == 0:
                        f4_ni = F4[index] + delta
                        f4_nr = f4_ni - CM
                        f4_rr = np.dot(MF.yaw_rot_mat(orientation),f4_nr.reshape(3,1)).reshape(3,)
                        F4[index] = f4_rr + CM
                    else:
                        F4[index] += delta
                else:
                    C3[index] = 0.0
                    sh3[index] = step_height
                    if i == 0:
                        f3_ni = F3[index] + delta
                        f3_nr = f3_ni - CM
                        f3_rr = np.dot(MF.yaw_rot_mat(orientation),f3_nr.reshape(3,1)).reshape(3,)
                        F3[index] = f3_rr + CM
                    else:
                        F3[index] += delta
                index += 1

            """
            print orientation
            print direction
            print delta
            print F1[index-1]
            print F2[index-1]
            print F3[index-1]
            print F4[index-1]
            """


    # stop
    Tk[index] = round(Tk[index-1] + 0.5,6)
    T_H[index] = round(T_H[index-1] + 0.5,6)
    F1[index] = F1[index-1]
    F2[index] = F2[index-1]
    F3[index] = F3[index-1]
    F4[index] = F4[index-1]
    H[index] = height

    return T_H, Tk, Tk, Tk, Tk, sh1, sh2, sh3, sh4, F1, F2, F3, F4, C1, C2, C3, C4, H, num_steps
