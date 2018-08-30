#!/usr/bin/env python
"""
A class used to solve for a trajectory using ipopt.
     - This serves as the structure needed to pass into ipopt
"""

__author__  	= "Joshua Hooks"
__email__   	= "hooksjrose@gmail.com"
__copyright__ 	= "Copyright 2017 RoMeLa"
__date__ 		= "Aug. 6th, 2018"

__version__ 	= "0.1.0"
__status__ 		= "Prototype"

"""
Classes to setup teleoperated trajectories. All trajectories are meant to be stitched together.
"""

from numpy import *
import matplotlib.pyplot as plt
import ipopt_problem as ipopt_problem
import footstep_planner as fp
from opt_structures import Constraint, Cost
#import Memory.MemoryManager as MM
import pdb

g = 9.81

class Structure(object):
    def __init__(self):

        #### Ipopt parameters ####
        self.nvar = 0 # number of variables to optimize for
        self.ncon = 0 # number of constraint equations
        self.nnjz = 0 # number of non-zero values in the gradient matrix of the constraint equations
        self.nnzh = 0 # number of non=zero values in the hessian matrix (Currently not used)
        self.x_L = 0  # Lower bound on variables
        self.x_U = 0  # Upper bound on variables
        self.g_L = 0  # Lower bound on constraint functions
        self.g_U = 0  # Upper bound on constraint functions

        #### Footstep schedule parameters ####
        self.H = 0    # List of heights at the beginning of each step
        self.num_steps = 0 # Number of steps taken in given trajectory
        self.T_H = 0  # List of times when the robot height should change
        self.T_F1 = 0   # Foot 1 list: Times when changes occur
        self.T_F2 = 0   # Foot 2 list:
        self.T_F3 = 0   # Foot 3 list
        self.T_F4 = 0   # Foot 4 list
        self.F1 = 0    # Foot 1 list: gives the foot position for a given time
        self.F2 = 0    # Foot 2 list:
        self.F3 = 0    # Foot 3 list:
        self.F4 = 0    # Foot 4 list:

        #### Initial and final conditions for center of mass ####
        self.CM_x_init = 0.0
        self.CM_y_init = 0.0
        self.dCM_x_init = 0.0
        self.dCM_y_init = 0.0
        self.CM_x_final = 0.0
        self.CM_y_final = 0.0
        self.dCM_x_final = 0.0
        self.dCM_y_final = 0.0

        #### Parameters for foot positions ####
        self.f1_x_init = 0.0    # Foot 1's initial x position
        self.f1_x_nom = 0.0     # Foot 1's nominal x distance from the center of mass, used for kinematic constraint
        self.f2_x_init = 0.0
        self.f2_x_nom = 0.0
        self.f3_x_init = 0.0
        self.f3_x_nom = 0.0
        self.f4_x_init = 0.0
        self.f4_x_nom = 0.0
        self.f1_y_init = 0.0    # Foot 1's initial y position
        self.f1_y_nom = 0.0     # Foot 1's nominal y distance from the center of mass, used for kinematic constraint
        self.f2_y_init = 0.0
        self.f2_y_nom = 0.0
        self.f3_y_init = 0.0
        self.f3_y_nom = 0.0
        self.f4_y_init = 0.0
        self.f4_y_nom = 0.0
        self.final_foot_pos = None

        #### Optimized variables ####
        self.x = 0.0    # array of optimized parameters [wc,wu] wc = cm polynomials, wu = load distribution
        self.wu_index = 0.0   # Starting index of wu
        self.wf1_index = 0     # Starting index of f1 positions
        self.wf2_index = 0
        self.wf3_index = 0
        self.wf4_index = 0

        #### Time parameters ####
        self.T = 0.0    # Total time for trajectory
        self.tp = 0.0   # polynomial time
        self.np = 0.0   # number of polynomials
        self.step_time = 0.0
        self.check_time = 0.0  # The time when the main thread should check if another trajectory has been generated
        self.heading = 0.0

        #### Constraints ####
        self.SCon = None
        self.DCon = None
        self.KCon = None
        self.LCon = None

        #### Costs ####
        self.ACost = None
        self.LCost = None

        self.h = 0.0    # height of the center of mass
        self.step_height = 0.0  # height of step

        self.new = 1.0

        self.solved = False

    """
    @classmethod
    def load_from_mem(cls):
        traj = cls()
        data = MM.TRAJECTORY_STATE.get()
        nvar = int(data['nvar'][0][0].copy())
        num_phases = int(data['num_phases'][0][0].copy())
        traj.check_time = data['transition_time'][0][0].copy()
        traj.wu_index = int(data['wu_index'][0][0].copy())
        traj.tp = data['tp'][0][0].copy()
        traj.T = data['T'][0][0].copy()
        traj.new = data['new'][0][0].copy()
        traj.heading = data['heading'][0][0].copy()

        traj.T_H = data['T_H'][:num_phases].reshape(num_phases).copy()
        traj.T_F1 = data['T_F1'][:num_phases].reshape(num_phases).copy()
        traj.T_F2 = data['T_F2'][:num_phases].reshape(num_phases).copy()
        traj.T_F3 = data['T_F3'][:num_phases].reshape(num_phases).copy()
        traj.T_F4 = data['T_F4'][:num_phases].reshape(num_phases).copy()
        traj.C_F1 = data['C_F1'][:num_phases].reshape(num_phases).copy()
        traj.C_F2 = data['C_F2'][:num_phases].reshape(num_phases).copy()
        traj.C_F3 = data['C_F3'][:num_phases].reshape(num_phases).copy()
        traj.C_F4 = data['C_F4'][:num_phases].reshape(num_phases).copy()
        traj.sh_F1 = data['sh_F1'][:num_phases].reshape(num_phases).copy()
        traj.sh_F2 = data['sh_F2'][:num_phases].reshape(num_phases).copy()
        traj.sh_F3 = data['sh_F3'][:num_phases].reshape(num_phases).copy()
        traj.sh_F4 = data['sh_F4'][:num_phases].reshape(num_phases).copy()
        traj.H = data['H'][:num_phases].reshape(num_phases).copy()
        traj.F1 = data['F1'][:num_phases][:].copy()
        traj.F2 = data['F2'][:num_phases][:].copy()
        traj.F3 = data['F3'][:num_phases][:].copy()
        traj.F4 = data['F4'][:num_phases][:].copy()
        traj.x = data['x'][:nvar].reshape(nvar).copy()
        return traj
    """
        
        
    def curr_a(self,t):
        """
        Get the CM polynomial coefficients for the specific time.
    
        :param x: Optimized values
        :param t: Current time
        :return: X and Y polynomial coefficients and their indices in the optimized value array.
        """
        index_c = int(round(t/self.tp,6))*10 # determine which polynomial to be using
    
        if t == self.T:
            index_c = index_c - 10
    
        a_x = self.x[index_c:index_c+5]
        a_y = self.x[index_c+5:index_c+10]
        da_x = arange(index_c,index_c+5)
        da_y = arange(index_c+5,index_c+10)
    
        return a_x, a_y, da_x, da_y


    def curr_lam(self,t):
        """
        Get load distribution values (lambda) for the specific time
    
        :param x: Optimized values
        :param t: Current time
        :return: A vector of lambda values for each foot and their indices in the optimized value array
        """
        index_u = int(round(t/(self.tp/2.0),6))*4 # determine which load distribution values to use
    
        lam = self.x[self.wu_index+index_u:self.wu_index+index_u+4]
        dlam = arange(self.wu_index+index_u,self.wu_index+index_u+4)
    
        return lam, dlam


    def curr_f(self,t,shift=True):
        """
        Get the feet position and if they are on the ground or not for the specific time
    
        :param x: Optimized values
        :param t: Current time
        :return: 2 objects for each foot containing: the foot position [x,y,z] and
                 whether the foot is on the ground or not (1.0 if on, 0.0 if off)
        """
        t = round(t,8)

        # Limb 1 trajectory
        index_t1 = searchsorted(self.T_F1,t)
        if index_t1 % 2 == 0:
            c_f1 = 1.0
        else:
            c_f1 = 0.0
        index_f1 = index_t1/2
        if c_f1:
            f1 = self.F1[index_f1]
        else:
            s = (t-self.T_F1[index_t1-1])*2*pi/(self.T_F1[index_t1] - self.T_F1[index_t1-1])
            f1_x = self.F1[index_f1+1][0] - self.F1[index_f1][0]
            f1_y = self.F1[index_f1+1][1] - self.F1[index_f1][1]
            f1_z = self.F1[index_f1+1][2] - self.F1[index_f1][2]
            x_pos = self.F1[index_f1][0] + f1_x*(s-sin(s))/(2*pi)
            y_pos = self.F1[index_f1][1] + f1_y*(s-sin(s))/(2*pi)
            step_h1 = self.step_height
            if s < pi:
                z_pos = self.F1[index_f1][2] + step_h1*(1-cos(s))/2.0
            else:
                z_pos = self.F1[index_f1][2] + f1_z + (step_h1-f1_z)*(1-cos(s))/2.0
            f1 = array([x_pos, y_pos, z_pos])
        df1 = array([self.wf1_index+2*index_f1, self.wf1_index+2*index_f1+1])

        # Limb 2 trajectory
        index_t2 = searchsorted(self.T_F2,t)
        if index_t2 % 2 == 0:
            c_f2 = 1.0
        else:
            c_f2 = 0.0
        index_f2 = index_t2/2
        if c_f2:
            f2 = self.F2[index_f2]
        else:
            s = (t-self.T_F2[index_t2-1])*2*pi/(self.T_F2[index_t2] - self.T_F2[index_t2-1])
            f2_x = self.F2[index_f2+1][0] - self.F2[index_f2][0]
            f2_y = self.F2[index_f2+1][1] - self.F2[index_f2][1]
            f2_z = self.F2[index_f2+1][2] - self.F2[index_f2][2]
            x_pos = self.F2[index_f2][0] + f2_x*(s-sin(s))/(2*pi)
            y_pos = self.F2[index_f2][1] + f2_y*(s-sin(s))/(2*pi)
            step_h2 = self.step_height
            if s < pi:
                z_pos = self.F2[index_f2][2] + step_h2*(1-cos(s))/2.0
            else:
                z_pos = self.F2[index_f2][2] + f2_z + (step_h2-f2_z)*(1-cos(s))/2.0
            f2 = array([x_pos, y_pos, z_pos])
        df2 = array([self.wf2_index+2*index_f2, self.wf2_index+2*index_f2+1])

        # Limb 3 trajectory
        index_t3 = searchsorted(self.T_F3,t)
        if index_t3 % 2 == 0:
            c_f3 = 1.0
        else:
            c_f3 = 0.0
        index_f3 = index_t3/2
        if c_f3:
            f3 = self.F3[index_f3]
        else:
            s = (t-self.T_F3[index_t3-1])*2*pi/(self.T_F3[index_t3] - self.T_F3[index_t3-1])
            f3_x = self.F3[index_f3+1][0] - self.F3[index_f3][0]
            f3_y = self.F3[index_f3+1][1] - self.F3[index_f3][1]
            f3_z = self.F3[index_f3+1][2] - self.F3[index_f3][2]
            x_pos = self.F3[index_f3][0] + f3_x*(s-sin(s))/(2*pi)
            y_pos = self.F3[index_f3][1] + f3_y*(s-sin(s))/(2*pi)
            step_h3 = self.step_height
            if s < pi:
                z_pos = self.F3[index_f3][2] + step_h3*(1-cos(s))/2.0
            else:
                z_pos = self.F3[index_f3][2] + f3_z + (step_h3-f3_z)*(1-cos(s))/2.0
            f3 = array([x_pos, y_pos, z_pos])
        df3 = array([self.wf3_index+2*index_f3, self.wf3_index+2*index_f3+1])

        # Limb 1 trajectory
        index_t4 = searchsorted(self.T_F4,t)
        if index_t4 % 2 == 0:
            c_f4 = 1.0
        else:
            c_f4 = 0.0
        index_f4 = index_t4/2
        if c_f4:
            f4 = self.F4[index_f4]
        else:
            s = (t-self.T_F4[index_t4-1])*2*pi/(self.T_F4[index_t4] - self.T_F4[index_t4-1])
            f4_x = self.F4[index_f4+1][0] - self.F4[index_f4][0]
            f4_y = self.F4[index_f4+1][1] - self.F4[index_f4][1]
            f4_z = self.F4[index_f4+1][2] - self.F4[index_f4][2]
            x_pos = self.F4[index_f4][0] + f4_x*(s-sin(s))/(2*pi)
            y_pos = self.F4[index_f4][1] + f4_y*(s-sin(s))/(2*pi)
            step_h4 = self.step_height
            if s < pi:
                z_pos = self.F4[index_f4][2] + step_h4*(1-cos(s))/2.0
            else:
                z_pos = self.F4[index_f4][2] + f4_z + (step_h4-f4_z)*(1-cos(s))/2.0
            f4 = array([x_pos, y_pos, z_pos])
        df4 = array([self.wf4_index+2*index_f4, self.wf4_index+2*index_f4+1])

        F = array([f1,f2,f3,f4])
        C = array([c_f1,c_f2,c_f3,c_f4])
        #Df = array([df1,df2,df3,df4])
    
        return F, C


    def curr_h(self,t):

        t = round(t,8)
        """
        if t == self.T_H[-1]:
            h = self.H[-1]
            ddh = 0.0
        else:
            index = searchsorted(self.T_H,t)
            H1 = self.H[index-1]
            H2 = self.H[index]
            T1 = self.T_H[index-1]
            T2 = self.T_H[index]

            s = 2*pi*(t-T1)/(T2-T1)
            h = H1 + (H2-H1)*(s-sin(s))/(2*pi)
            ddh = 2*pi*(H2-H1)*sin(s)/(T2-T1)**2
        """
        h = self.H[0]
        ddh = 0.0

        return h, ddh


    def solve(self, initial_CM, initial_dCM, final_CM, final_dCM, TP, step_height, foot_constraint, warm_start=None):

        self.step_height = step_height
    
        self.f1_x_nom = self.F1[0][0] - initial_CM[0]
        self.f1_y_nom = self.F1[0][1] - initial_CM[1]
    
        self.f2_x_nom = self.F2[0][0] - initial_CM[0]
        self.f2_y_nom = self.F2[0][1] - initial_CM[1]
    
        self.f3_x_nom = self.F3[0][0] - initial_CM[0]
        self.f3_y_nom = self.F3[0][1] - initial_CM[1]
    
        self.f4_x_nom = self.F4[0][0] - initial_CM[0]
        self.f4_y_nom = self.F4[0][1] - initial_CM[1]

        foot_nominals = array([[self.f1_x_nom, self.f1_y_nom],
                               [self.f2_x_nom, self.f2_y_nom],
                               [self.f3_x_nom, self.f3_y_nom],
                               [self.f4_x_nom, self.f4_y_nom]])
    
        self.CM_x_init = initial_CM[0]
        self.CM_y_init = initial_CM[1]
        self.dCM_x_init = initial_dCM[0]
        self.dCM_y_init = initial_dCM[1]

        self.CM_x_final = final_CM[0]
        self.CM_y_final = final_CM[1]
        self.dCM_x_final = final_dCM[0]
        self.dCM_y_final = final_dCM[1]

        self.tp = TP # Duration of polynomials
        self.np = int(round(self.T/self.tp,3))+1 # number of polynomials
    
        wc = ones(5*2*self.np)*0.1
        wu = ones(4*(2*self.np+1))*0.25

        if warm_start != None:
            wc[:len(warm_start[0])] = warm_start[0]
            wu[:len(warm_start[1])] = warm_start[1]
    
        self.wu_index = len(wc)
    
        x0 = concatenate((wc,wu), axis=0)
        
        self.x = x0

        self.SCon = StitchingConstraint(self.tp)
        self.DCon = DynamicConstraint()
        self.KCon = KinematicConstraint(foot_nominals, foot_constraint)
        self.LCon = LoadConstraint()

        acceleration_gain = 0.0
        load_gain = 1.0
        self.ACost = AccelerationCost(acceleration_gain)
        self.LCost = LoadCost(load_gain)
    
        self.nnjz = (self.SCon.nnjz + 2*self.DCon.nnjz + self.KCon.nnjz + 2*self.LCon.nnjz)*(self.np-1)  # Number of non-zero values in the jacobian matrix
        self.nvar = len(x0)  # Number of variables to solve for
        self.nnzh = 0  # Not used because no Hessian was derive
        self.ncon = (self.SCon.ncon + 2*self.DCon.ncon + self.KCon.ncon + 2*self.LCon.ncon)*(self.np-1)  # number of constraints
    
        # Upper and lower bounds on variables
        self.x_L = ones((self.nvar), dtype=float_) * -2.0*pow(10.0, 19)
        self.x_U = ones((self.nvar), dtype=float_) * 2.0*pow(10.0, 19)
    
        # Constrain load distribution between 0 and 1
        self.x_L[self.wu_index:] = zeros(len(wu), dtype=float_)
        self.x_U[self.wu_index:] = ones(len(wu), dtype=float_)

        # Force the initial postion and velocity to be correct
        self.x_L[0] = self.x_U[0] = initial_CM[0]
        self.x_L[1] = self.x_U[1] = initial_dCM[0]
        self.x_L[5] = self.x_U[5] = initial_CM[1]
        self.x_L[6] = self.x_U[6] = initial_dCM[1]
        self.x_L[-10] = self.x_U[-10] = final_CM[0]
        self.x_L[-9] = self.x_U[-9] = final_dCM[0]
        self.x_L[-5] = self.x_U[-5] = final_CM[1]
        self.x_L[-4] = self.x_U[-4] = final_dCM[1]
    
        # Upper and lower bounds on constraints
        self.g_L = zeros(self.ncon)
        self.g_U = zeros(self.ncon)

        t = 0.0
        index = 0
        for j in range(self.np-1):
            _, C = self.curr_f(t)
            self.g_L[index:index+self.SCon.ncon], self.g_U[index:index+self.SCon.ncon] = self.SCon.boundaries()
            index += self.SCon.ncon
            self.g_L[index:index+self.DCon.ncon], self.g_U[index:index+self.DCon.ncon] = self.DCon.boundaries()
            index += self.DCon.ncon
            self.g_L[index:index+self.KCon.ncon], self.g_U[index:index+self.KCon.ncon] = self.KCon.boundaries()
            index += self.KCon.ncon
            self.g_L[index:index+self.LCon.ncon], self.g_U[index:index+self.LCon.ncon] = self.LCon.boundaries(C)
            index += self.LCon.ncon

            # Perform 2 constraint checks in the middle of the time period
            t_bar = round(self.tp/2.0,6)
            t = round(t+t_bar,6)
            _, C = self.curr_f(t)
            self.g_L[index:index+self.DCon.ncon], self.g_U[index:index+self.DCon.ncon] = self.DCon.boundaries()
            index += self.DCon.ncon
            self.g_L[index:index+self.LCon.ncon], self.g_U[index:index+self.LCon.ncon] = self.LCon.boundaries(C)
            index += self.LCon.ncon

            t = round(t+t_bar,6)
    
        ipopt_problem.solve(self)


    def plot_solution(self, trajectories = None):

        frames = int(self.T/0.005)
        t = linspace(0.0,self.T,frames)
        tk = 0.0

        lam, dlam = self.curr_lam(t[0])
        F, C = self.curr_f(t[0])
        a_x, a_y, da_x, da_y = self.curr_a(t[0])

        plt.ion()
        fig1 = plt.figure(1)
        ax = fig1.add_subplot(1,1,1)
        ax.set_ylim((-0.7,1.3))
        ax.set_xlim((-0.7,1.3))
        CM_X = []
        CM_Y = []
        line, = ax.plot(CM_X, CM_Y)
        point1, = ax.plot(F[0][0], F[0][1], 'x')
        point2, = ax.plot(F[1][0], F[1][1], 'x')
        point3, = ax.plot(F[2][0], F[2][1], 'x')
        point4, = ax.plot(F[3][0], F[3][1], 'x')

        #ux = lam[0]*f1[0] + lam[1]*f2[0] + lam[2]*f3[0] + lam[3]*f4[0]
        #uy = lam[0]*f1[1] + lam[1]*f2[1] + lam[2]*f3[1] + lam[3]*f4[1]

        n_trajectories = 0
        if trajectories != None:
            n_trajectories = len(trajectories)

        first = True
        for w in range(1+n_trajectories):
            if not first:
                traj = trajectories[w-1]
                frames = int(traj.T/0.005)
                t = linspace(0.0,traj.T,frames)
                tk = 0.0

            for i in range(0,len(t)):
                if n_trajectories > 0 and w < n_trajectories:
                    if first:
                        if t[i] > self.T_F1[5]:
                            break
                    else:
                        if t[i] > trajectories[w].T_F1[4]:
                            break
                if first:
                    lam, dlam = self.curr_lam(t[i])
                    F, C = self.curr_f(t[i])
                    a_x, a_y, da_x, da_y = self.curr_a(t[i])
                    if t[i] >= tk + self.tp:
                        tk = tk + self.tp
                else:
                    lam, dlam = traj.curr_lam(t[i])
                    F, C, _ = traj.curr_f(t[i])
                    a_x, a_y, da_x, da_y = traj.curr_a(t[i])
                    if t[i] >= tk + traj.tp:
                        tk = tk + traj.tp


                CM_x = a_x[0]
                CM_y = a_y[0]
                dCM_x = 0.0
                dCM_y = 0.0
                for j in range(1,5):
                    CM_x = CM_x + a_x[j]*(t[i]-tk)**j
                    CM_y = CM_y + a_y[j]*(t[i]-tk)**j
                    dCM_x = dCM_x + j*a_x[j]*(t[i]-tk)**(j-1)
                    dCM_y = dCM_y + j*a_y[j]*(t[i]-tk)**(j-1)

                CM_X.append(CM_x)
                CM_Y.append(CM_y)
                line.set_xdata(CM_X)
                line.set_ydata(CM_Y)


                #ux = lam[0]*f1[0] + lam[1]*f2[0] + lam[2]*f3[0] + lam[3]*f4[0]
                #uy = lam[0]*f1[1] + lam[1]*f2[1] + lam[2]*f3[1] + lam[3]*f4[1]

                #ax.plot(ux,uy,'.')

                if C[0] == 1.0:
                    point1.set_xdata(F[0][0])
                    point1.set_ydata(F[0][1])
                else:
                    point1.set_xdata('nan')
                    point1.set_ydata('nan')

                if C[1] == 1.0:
                    point2.set_xdata(F[1][0])
                    point2.set_ydata(F[1][1])
                else:
                    point2.set_xdata('nan')
                    point2.set_ydata('nan')

                if C[2] == 1.0:
                    point3.set_xdata(F[2][0])
                    point3.set_ydata(F[2][1])
                else:
                    point3.set_xdata('nan')
                    point3.set_ydata('nan')

                if C[3] == 1.0:
                    point4.set_xdata(F[3][0])
                    point4.set_ydata(F[3][1])
                else:
                    point4.set_xdata('nan')
                    point4.set_ydata('nan')

                fig1.canvas.draw()

            first = False


        CM_x = a_x[0]
        CM_y = a_y[0]
        for j in range(1,5):
            CM_x = CM_x + a_x[j]*(t[i]-tk)**j
            CM_y = CM_y + a_y[j]*(t[i]-tk)**j

        plt.ioff()
        ax.plot(CM_x, CM_y, 'o')
        #ax.set_xticks(arange(-0.3, 0.3, 0.02))
        #ax.set_yticks(arange(-0.1, 0.3, 0.02))
        #plt.grid()
        plt.show()


    def plot3D(self):

        frames = int(self.T/0.01)
        t = linspace(0.0,self.T,frames)
        tk = 0.0

        lam, dlam = self.curr_lam(t[0])
        f1, f2, f3, f4, c_f1, c_f2, c_f3, c_f4 = self.curr_f(t[0],False)
        a_x, a_y, da_x, da_y = self.curr_a(t[0])

        plt.ion()
        fig1 = plt.figure(1)
        ax = fig1.add_subplot(111, projection='3d')
        ax.set_xlim3d((-0.4,1.2))
        ax.set_ylim3d((-0.4,1.2))
        ax.set_zlim3d((0.0,1.0))
        CM_X = []
        CM_Y = []
        CM_Z = []
        line, = ax.plot(CM_X, CM_Y, CM_Z)

        FX_x = [f1[0], f2[0], f3[0], f4[0]]
        FY_x = [f1[1], f2[1], f3[1], f4[1]]
        FZ_x = [f1[2], f2[2], f3[2], f4[2]]
        FX_o = []
        FY_o = []
        FZ_o = []
        points, = ax.plot(FX_x,FY_x,FZ_x,'x')

        for i in range(0,len(t),5):
            lam, dlam = self.curr_lam(t[i])
            f1, f2, f3, f4, c_f1, c_f2, c_f3, c_f4 = self.curr_f(t[i],False)
            h, dhh = self.curr_h(t[i])
            a_x, a_y, da_x, da_y = self.curr_a(t[i])
            if t[i] >= tk + self.tp:
                tk = tk + self.tp


            CM_x = a_x[0]
            CM_y = a_y[0]
            for j in range(1,5):
                CM_x = CM_x + a_x[j]*(t[i]-tk)**j
                CM_y = CM_y + a_y[j]*(t[i]-tk)**j

            CM_X.append(CM_x)
            CM_Y.append(CM_y)
            CM_Z.append(h)
            line.set_xdata(CM_X)
            line.set_ydata(CM_Y)
            line.set_3d_properties(CM_Z)

            FX_x.append(f1[0])
            FY_x.append(f1[1])
            FZ_x.append(f1[2])

            FX_x.append(f2[0])
            FY_x.append(f2[1])
            FZ_x.append(f2[2])

            FX_x.append(f3[0])
            FY_x.append(f3[1])
            FZ_x.append(f3[2])

            FX_x.append(f4[0])
            FY_x.append(f4[1])
            FZ_x.append(f4[2])


            points.set_xdata(FX_x)
            points.set_ydata(FY_x)
            points.set_3d_properties(FZ_x)

            fig1.canvas.draw()


        CM_x = a_x[0]
        CM_y = a_y[0]
        for j in range(1,5):
            CM_x = CM_x + a_x[j]*(t[i]-tk)**j
            CM_y = CM_y + a_y[j]*(t[i]-tk)**j

        CM_z = h

        plt.ioff()
        #ax.plot(CM_x, CM_y, CM_z, 'o')
        #ax.set_xticks(arange(-0.3, 0.3, 0.02))
        #ax.set_yticks(arange(-0.1, 0.3, 0.02))
        #plt.grid()
        plt.show()


    """
    def write2mem(self):
        data = {}
        data['nvar'] = array([self.nvar])
        data['num_phases'] = array(len(self.C_F1))
        data['transition_time'] = array(self.check_time)
        data['wu_index'] = array(self.wu_index)
        data['tp'] = array(self.tp)
        data['T'] = array(self.T)
        data['heading'] = array(self.heading)

        self.F1 = zeros(15-len(self.C_F1))
        self.F2 = zeros((15-len(self.C_F1),3))
        self.F3 = zeros(1500-len(self.x))

        data['T_H'] = concatenate((self.T_H, self.F1))
        data['T_F1'] = concatenate((self.T_F1, self.F1))
        data['T_F2'] = concatenate((self.T_F2, self.F1))
        data['T_F3'] = concatenate((self.T_F3, self.F1))
        data['T_F4'] = concatenate((self.T_F4, self.F1))
        data['C_F1'] = concatenate((self.C_F1, self.F1))
        data['C_F2'] = concatenate((self.C_F2, self.F1))
        data['C_F3'] = concatenate((self.C_F3, self.F1))
        data['C_F4'] = concatenate((self.C_F4, self.F1))
        data['sh_F1'] = concatenate((self.sh_F1, self.F1))
        data['sh_F2'] = concatenate((self.sh_F2, self.F1))
        data['sh_F3'] = concatenate((self.sh_F3, self.F1))
        data['sh_F4'] = concatenate((self.sh_F4, self.F1))
        data['H'] = concatenate((self.H, self.F1))
        data['F1'] = concatenate((self.F1, self.F2), axis=0)
        data['F2'] = concatenate((self.F2, self.F2), axis=0)
        data['F3'] = concatenate((self.F3, self.F2), axis=0)
        data['F4'] = concatenate((self.F4, self.F2), axis=0)
        data['x'] = concatenate((self.x, self.F3))
        data['new'] = ones((1,1))

        MM.TRAJECTORY_STATE.set(data)
    """


class StitchingConstraint(Constraint):

    def __init__(self, tp):

        self.ncon = 4
        self.nnjz = 22
        self.T = array([1.0, tp, tp**2.0, tp**3.0, tp**4.0])
        self.dT = array([0.0, 1.0, 2.0*tp, 3.0*tp**2.0, 4.0*tp**3.0])

    def evaluate(self, a_x0, a_y0, a_x1, a_y1):
        """
        :param a_x0: x coefficients for current time
        :param a_y0: y coefficients for current time
        :param a_x1: x coefficients for next polynomial
        :param a_y1: y coefficients for next polynomial
        :return: comparison between last value in current polynomial to first value in next polynomial
        """
        return array([dot(a_x0,self.T) - a_x1[0], dot(a_y0,self.T) - a_y1[0],
                      dot(a_x0,self.dT) - a_x1[1], dot(a_y0,self.dT) - a_y1[1],
                      ])

    def derivative(self):
        return concatenate((self.T, array([-1.0]), self.T, array([-1.0]),
                            self.dT[1:], array([-1.0]), self.dT[1:], array([-1.0])))

    def derivative_structure(self, da_x0, da_y0, da_x1, da_y1):
        col1 = concatenate((da_x0, array([da_x1[0]])))
        col2 = concatenate((da_y0, array([da_y1[0]])))
        col3 = concatenate((da_x0[1:], array([da_x1[1]])))
        col4 = concatenate((da_y0[1:], array([da_y1[1]])))
        return array([col1, col2, col3, col4])

    def boundaries(self):
        return zeros(self.ncon), zeros(self.ncon)


class DynamicConstraint(Constraint):

    def __init__(self):

        self.ncon = 2
        self.nnjz = 18

    def evaluate(self, t_bar, a_x, a_y, F, C, lam, h, ddh):
        T = array([1.0, t_bar, t_bar**2.0, t_bar**3.0, t_bar**4.0])
        ddT = array([0.0, 0.0, 2.0, 6.0*t_bar, 12.0*t_bar**2.0])

        u_x = dot(lam,F[:,0]*C)
        u_y = dot(lam,F[:,1]*C)

        cm_x = dot(a_x,T)
        cm_ddx = dot(a_x,ddT)
        cm_y = dot(a_y,T)
        cm_ddy = dot(a_y,ddT)
        return array([cm_ddx - (g+ddh)*(cm_x - u_x)/h, cm_ddy - (g+ddh)*(cm_y - u_y)/h])

    def derivative(self, t_bar, F, C, lam, h, ddh):
        T = array([1.0, t_bar, t_bar**2.0, t_bar**3.0, t_bar**4.0])
        ddT = array([0.0, 0.0, 2.0, 6.0*t_bar, 12.0*t_bar**2.0])
        return array([ddT[0]-T[0]*(g+ddh)/h,
                      ddT[1]-T[1]*(g+ddh)/h,
                      ddT[2]-T[2]*(g+ddh)/h,
                      ddT[3]-T[3]*(g+ddh)/h,
                      ddT[4]-T[4]*(g+ddh)/h,
                      C[0]*F[0][0]*(g+ddh)/h,
                      C[1]*F[1][0]*(g+ddh)/h,
                      C[2]*F[2][0]*(g+ddh)/h,
                      C[3]*F[3][0]*(g+ddh)/h,
                      ddT[0]-T[0]*(g+ddh)/h,
                      ddT[1]-T[1]*(g+ddh)/h,
                      ddT[2]-T[2]*(g+ddh)/h,
                      ddT[3]-T[3]*(g+ddh)/h,
                      ddT[4]-T[4]*(g+ddh)/h,
                      C[0]*F[0][1]*(g+ddh)/h,
                      C[1]*F[1][1]*(g+ddh)/h,
                      C[2]*F[2][1]*(g+ddh)/h,
                      C[3]*F[3][1]*(g+ddh)/h
                      ])

    def derivative_structure(self, da_x, da_y, dlam):
        col1 = concatenate((da_x, dlam))
        col2 = concatenate((da_y, dlam))
        return array([col1, col2])

    def boundaries(self):
        return zeros(2), zeros(2)


class KinematicConstraint(Constraint):

    def __init__(self, f_initial, foot_boundary):

        self.ncon = 8
        self.nnjz = 8
        self.f1_x_nom = f_initial[0][0]
        self.f1_y_nom = f_initial[0][1]
        self.f2_x_nom = f_initial[1][0]
        self.f2_y_nom = f_initial[1][1]
        self.f3_x_nom = f_initial[2][0]
        self.f3_y_nom = f_initial[2][1]
        self.f4_x_nom = f_initial[3][0]
        self.f4_y_nom = f_initial[3][1]
        self.foot_boundary = foot_boundary

    def evaluate(self, a_x, a_y, F, C):
        return array([C[0]*(F[0][0] - a_x[0] - self.f1_x_nom),
                      C[0]*(F[0][1] - a_y[0] - self.f1_y_nom),
                      C[1]*(F[1][0] - a_x[0] - self.f2_x_nom),
                      C[1]*(F[1][1] - a_y[0] - self.f2_y_nom),
                      C[2]*(F[2][0] - a_x[0] - self.f3_x_nom),
                      C[2]*(F[2][1] - a_y[0] - self.f3_y_nom),
                      C[3]*(F[3][0] - a_x[0] - self.f4_x_nom),
                      C[3]*(F[3][1] - a_y[0] - self.f4_y_nom),
                      ])

    def derivative(self, C):
        return array([-C[0], -C[0], -C[1], -C[1], -C[2], -C[2], -C[3], -C[3]])

    def derivative_structure(self, da_x, da_y):
        col1 = array([da_x[0]])
        col2 = array([da_y[0]])
        col3 = array([da_x[0]])
        col4 = array([da_y[0]])
        col5 = array([da_x[0]])
        col6 = array([da_y[0]])
        col7 = array([da_x[0]])
        col8 = array([da_y[0]])
        return array([col1, col2, col3, col4, col5, col6, col7, col8])


    def boundaries(self):
        return ones(self.ncon)*-self.foot_boundary, ones(self.ncon)*self.foot_boundary


class LoadConstraint(Constraint):

    def __init__(self):
        self.ncon = 5
        self.nnjz = 8

    def evaluate(self, lam):
        return array([lam.sum(), lam[0], lam[1], lam[2], lam[3]])

    def derivative(self):
        return ones(8)

    def derivative_structure(self, dlam):
        col1 = array([dlam[0], dlam[1], dlam[2], dlam[3]])
        col2 = array([dlam[0]])
        col3 = array([dlam[1]])
        col4 = array([dlam[2]])
        col5 = array([dlam[3]])
        return array([col1, col2, col3, col4, col5])

    def boundaries(self, C):
        return array([1.0, 0.0, 0.0, 0.0, 0.0]), array([1.0, C[0], C[1], C[2], C[3]])


class AccelerationCost(Cost):

    def __init__(self, gain):
        self.gain = gain

    def evaluate(self, t_bar, a_x, a_y):
        cm_acc_x_squared = 4.0*a_x[2]**2.0 + 36.0*a_x[3]**2.0*t_bar**2.0 + 144.0*a_x[4]**2.0*t_bar**4.0 + 24.0*a_x[2]*a_x[3]*t_bar + \
                           48.0*a_x[2]*a_x[4]*t_bar**2.0 + 72.0*a_x[3]*a_x[4]*t_bar**3.0

        cm_acc_y_squared = 4.0*a_y[2]**2.0 + 36.0*a_y[3]**2.0*t_bar**2.0 + 144.0*a_y[4]**2.0*t_bar**4.0 + 24.0*a_y[2]*a_y[3]*t_bar + \
                           48.0*a_y[2]*a_y[4]*t_bar**2.0 + 72.0*a_y[3]*a_y[4]*t_bar**3.0

        return self.gain*(cm_acc_x_squared + cm_acc_y_squared)

    def derivative(self, t_bar, a_x, a_y, da_x, da_y):
        index = array([da_x[2], da_x[3], da_x[4], da_y[2], da_y[3], da_y[4],])
        values = array([self.gain*(8.0*a_x[2] + 24.0*a_x[3]*t_bar + 48.0*a_x[4]*t_bar**2.0),
                        self.gain*(72.0*a_x[3]*t_bar**2.0 + 24.0*a_x[2]*t_bar + 72.0*a_x[4]*t_bar**3.0),
                        self.gain*(288*a_x[4]*t_bar**4.0 + 48.0*a_x[2]*t_bar**2.0 + 72.0*a_x[3]*t_bar**3.0),
                        self.gain*(8.0*a_y[2] + 24.0*a_y[3]*t_bar + 48.0*a_y[4]*t_bar**2.0),
                        self.gain*(72.0*a_y[3]*t_bar**2.0 + 24.0*a_y[2]*t_bar + 72.0*a_y[4]*t_bar**3.0),
                        self.gain*(288*a_y[4]*t_bar**4.0 + 48.0*a_y[2]*t_bar**2.0 + 72.0*a_y[3]*t_bar**3.0)])

        return index, values


class LoadCost(Cost):

    def __init__(self, gain):
        self.gain = gain

    def evaluate(self, lam, C):
        n = C.sum()
        lam_star = array([C[0]/n, C[1]/n, C[2]/n, C[3]/n])
        error = lam - lam_star
        return self.gain*dot(error, error)

    def derivative(self, lam, dlam, C):
        index = array([dlam[0], dlam[1], dlam[2], dlam[3]])
        n = C.sum()
        lam_star = array([C[0]/n, C[1]/n, C[2]/n, C[3]/n])
        error = lam - lam_star
        values = self.gain*(2.0*error)

        return index, values



class Creep(Structure):

    def __init__(self, initial_CM, initial_dCM, final_CM, final_dCM, initial_foot_pos, step_dis, step_time,
                 step_height, height, start, direction, orientation, right, warm_start=None):

        self.T_H, self.T_F1, self.T_F2, self.T_F3, self.T_F4, self.sh_F1, self.sh_F2, self.sh_F3, self.sh_F4, self.F1, \
        self.F2, self.F3, self.F4, self.C_F1, self.C_F2, self.C_F3, self.C_F4,  self.H, self.num_steps \
        = fp.teleop_creep(start, initial_foot_pos, direction, orientation, step_dis, step_time, step_height, height, right)

        self.h = height

        TP = 0.05
        foot_constraint = 1.0

        if start:
            self.check_time = self.T_F1[5]
        else:
            self.check_time = self.T_F1[4]

        self.solve(initial_CM, initial_dCM, final_CM, final_dCM, initial_foot_pos, height, TP, step_height, foot_constraint, warm_start)


if __name__ == '__main__':
    test = Structure()

    test.T_F1 = array([0.5, 0.7, 1.3, 1.5])
    test.F1 = array([[0.5, 0.0, 0.0],
                     [0.6, 0.0, 0.0],
                     [0.7, 0.0, 0.0]])

    test.T_F2 = array([0.7, 0.9, 1.5, 1.7])
    test.F2 = array([[0.5, 0.4, 0.0],
                     [0.6, 0.4, 0.0],
                     [0.7, 0.4, 0.0]])

    test.T_F3 = array([0.9, 1.1, 1.7, 1.9])
    test.F3 = array([[0.0, 0.4, 0.0],
                     [0.1, 0.4, 0.0],
                     [0.2, 0.4, 0.0]])

    test.T_F4 = array([1.1, 1.3, 1.9, 2.1])
    test.F4 = array([[0.0, 0.0, 0.0],
                     [0.1, 0.0, 0.0],
                     [0.2, 0.0, 0.0]])
    test.T = 2.3

    test.T_H = array([0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 1.7])
    test.H =   array([0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5])
    test.step_height = 0.1

    i_CM = array([0.25, 0.2])
    i_dCM = array([0.0, 0.0])
    f_CM = array([0.45, 0.2])
    f_dCM = array([0.0, 0.0])

    TP = 0.05
    step_height = 0.1

    test.solve(i_CM, i_dCM, f_CM, f_dCM, TP, step_height, 0.5)
    test.plot_solution()

