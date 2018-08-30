#!/usr/bin/env python
"""
Contains all function to use the pyiopt libraries.
 - Needs to have the trajectory.structure class passed in
"""

__author__  	= "Joshua Hooks"
__email__   	= "hooksjrose@gmail.com"
__copyright__ 	= "Copyright 2017 RoMeLa"
__date__ 		= "Jan. 8th, 2018"

__version__ 	= "0.1.0"
__status__ 		= "Prototype"

#!/usr/bin/env python
"""
Solves for a CM trajectory that satisfies ZMP constraints given
initial pose and final position. Uses IPOpt to minimize the acceleration of the
CM and satisfies all constraints.

Method derived from:
" Fast P.Trajectory Optimization for Legged Robots using Vertex-based ZMP Constraints"


 *** Note the global variables in this can be confusing. Not sure how to do it with a class though
     becasue IPOPP.T needs the class handles passed in with only one variable.
"""
__author__  	= "Joshua Hooks"
__email__   	= "hooksjrose@gmail.com"
__copyright__ 	= "Copyright 2017 RoMeLa"
__date__ 		= "Sept. 7, 2017"

__version__ 	= "0.1.0"
__status__ 		= "Prototype"

import pyipopt
import matplotlib.pyplot as plt
from numpy import *
import pdb
import time

# Problem parameters
P = None

g = 9.81

def solve (p):
    global P
    P = p
    t0 = time.time()

    # Create the problem for Ipopt to solve
    nlp = pyipopt.create(P.nvar, P.x_L, P.x_U, P.ncon, P.g_L, P.g_U, P.nnjz, P.nnzh, eval_f, eval_grad_f, eval_g, eval_jac_g)

    # Set Ipopt solve options
    nlp.int_option("print_level",5)
    nlp.num_option("acceptable_tol",1000.0)
    nlp.num_option("tol",1000.0)
    nlp.num_option("dual_inf_tol",10000.0)
    nlp.num_option("compl_inf_tol",10000.0)
    nlp.int_option("max_iter", 1000)
    #nlp.str_option("derivative_test","first-order")

    stitching_constraint(0.0, 0)

    print "Calling solve"

    x, zl, zu, constrain_mu, obj, status = nlp.solve(P.x) # Solve the NLP problem with Ipopt
    nlp.close()
    if status == 0:
        solve_time = time.time() - t0
        print "Solution found in: ",solve_time," (sec)"
        P.x = x
        P.solved = True
    else:
        print "Failed to find solution!"
        P.x = x
        P.solved = False


def eval_f(x, user_data = None):
    """
    Cost function: the sum of the L2 norm of the difference between the load distribution vector and the
                   optimal load distribution vector for all time

    :param x: variable array
    :param user_data:
    :return:  Evaluation of the cost function
    """

    f = 0.0
    t = 0.0
    P.x = x
    tk = 0.0

    load_gain = 1.0
    acc_gain =0.0
    
    """
    used for stepping over
    load_gain = 1.0
    acc_gain = 0.0
    """

    for i in range(2*P.np+1):
        _, C, _ = P.curr_f(t)
        a_x, a_y, da_x, da_y = P.curr_a(t)
        lam, dlam = P.curr_lam(t)

        if t >= tk + P.tp:
                tk = tk + P.tp

        tb = t-tk

        CM_acc_X_squared = 4.0*a_x[2]**2.0 + 36.0*a_x[3]**2.0*tb**2.0 + 144.0*a_x[4]**2.0*tb**4.0 + 24.0*a_x[2]*a_x[3]*tb + \
                           48.0*a_x[2]*a_x[4]*tb**2.0 + 72.0*a_x[3]*a_x[4]*tb**3.0

        CM_acc_Y_squared = 4.0*a_y[2]**2.0 + 36.0*a_y[3]**2.0*tb**2.0 + 144.0*a_y[4]**2.0*tb**4.0 + 24.0*a_y[2]*a_y[3]*tb + \
                           48.0*a_y[2]*a_y[4]*tb**2.0 + 72.0*a_y[3]*a_y[4]*tb**3.0

        n = C.sum()
        lam_star = array([C[0]/n, C[1]/n, C[2]/n, C[3]/n])
        diff = lam - lam_star
        f = f + load_gain*dot(diff,diff) + acc_gain*(CM_acc_X_squared + CM_acc_Y_squared)
        t = round(t+P.tp/2.0,8)

    return f


def eval_grad_f(x, user_data = None):
    """
    Gradient of the Cost function

    :param x: variable array
    :param user_data:
    :return: gradient of the cost function
    """
    P.x = x
    grad_f = zeros(len(x))
    t = 0.0
    tk = 0.0

    load_gain = 1.0
    acc_gain = 0.0

    """
    used for stepping over
    load_gain = 1.0
    acc_gain = 0.0
    """

    for i in range(2*P.np+1):
        _, C, _ = P.curr_f(t)
        a_x, a_y, da_x, da_y = P.curr_a(t)
        lam, dlam = P.curr_lam(t)

        if t >= tk + P.tp:
                tk = tk + P.tp

        tb = t-tk

        grad_f[da_x[2]] = acc_gain*(8.0*a_x[2] + 24.0*a_x[3]*tb + 48.0*a_x[4]*tb**2.0)
        grad_f[da_x[3]] = acc_gain*(72.0*a_x[3]*tb**2.0 + 24.0*a_x[2]*tb + 72.0*a_x[4]*tb**3.0)
        grad_f[da_x[4]] = acc_gain*(288*a_x[4]*tb**4.0 + 48.0*a_x[2]*tb**2.0 + 72.0*a_x[3]*tb**3.0)

        grad_f[da_y[2]] = acc_gain*(8.0*a_y[2] + 24.0*a_y[3]*tb + 48.0*a_y[4]*tb**2.0)
        grad_f[da_y[3]] = acc_gain*(72.0*a_y[3]*tb**2.0 + 24.0*a_y[2]*tb + 72.0*a_y[4]*tb**3.0)
        grad_f[da_y[4]] = acc_gain*(288*a_y[4]*tb**4.0 + 48.0*a_y[2]*tb**2.0 + 72.0*a_y[3]*tb**3.0)

        n = C.sum()
        lam_star = array([C[0]/n, C[1]/n, C[2]/n, C[3]/n])
        diff = lam - lam_star
        t = round(t + P.tp/2.0,8)
        grad_f[dlam[0]:dlam[-1]+1] = load_gain*(2.0*diff)

    return grad_f


def stitching_constraint(t, derivative, flag=False):

    T = array([1.0, P.tp, P.tp**2.0, P.tp**3.0, P.tp**4.0])
    dT = array([0.0, 1.0, 2.0*P.tp, 3.0*P.tp**2.0, 4.0*P.tp**3.0])

    if derivative == 0:
        a_x0, a_y0, da_x0, da_y0 = P.curr_a(t)
        a_x1, a_y1, da_x1, da_y1 = P.curr_a(t+P.tp)
        return array([dot(a_x0,T) - a_x1[0], dot(a_y0,T) - a_y1[0],
                      dot(a_x0,dT) - a_x1[1], dot(a_y0,dT) - a_y1[1],
                      ])
    elif derivative == 1:
        if flag:
            a_x0, a_y0, da_x0, da_y0 = P.curr_a(t)
            a_x1, a_y1, da_x1, da_y1 = P.curr_a(t+P.tp)

            col1 = concatenate((da_x0, array([da_x1[0]])))
            col2 = concatenate((da_y0, array([da_y1[0]])))
            col3 = concatenate((da_x0[1:], array([da_x1[1]])))
            col4 = concatenate((da_y0[1:], array([da_y1[1]])))
            return array([col1, col2, col3, col4])
        else:
            return concatenate((T, array([-1.0]),
                                T, array([-1.0]),
                                dT[1:], array([-1.0]),
                                dT[1:], array([-1.0])))


def dynamic_constraint(t, t_bar, derivative, flag=False):

    T = array([1.0, t_bar, t_bar**2.0, t_bar**3.0, t_bar**4.0])
    ddT = array([0.0, 0.0, 2.0, 6.0*t_bar, 12.0*t_bar**2.0])

    if derivative == 0:
        a_x, a_y, _, _ = P.curr_a(t)
        F, C, _ = P.curr_f(t)
        lam, _ = P.curr_lam(t)
        h, ddh = P.curr_h(t)

        u_x = dot(lam,F[:,0]*C)
        u_y = dot(lam,F[:,1]*C)

        cm_x = dot(a_x,T)
        cm_ddx = dot(a_x,ddT)
        cm_y = dot(a_y,T)
        cm_ddy = dot(a_y,ddT)

        return array([cm_ddx - (g+ddh)*(cm_x - u_x)/h,
                      cm_ddy - (g+ddh)*(cm_y - u_y)/h
                      ])

    if derivative == 1:
        if flag:
            a_x, a_y, da_x, da_y = P.curr_a(t)
            F, C, dF = P.curr_f(t)
            lam, dlam = P.curr_lam(t)

            col1 = concatenate((da_x, dlam, dF[:,0]))
            col2 = concatenate((da_y, dlam, dF[:,1]))

            return array([col1, col2])
        else:
            a_x, a_y, _, _ = P.curr_a(t)
            F, C, _ = P.curr_f(t)
            lam, _ = P.curr_lam(t)
            h, ddh = P.curr_h(t)

            return array([ddT[0]-T[0]*(g+ddh)/h,
                          ddT[1]-T[1]*(g+ddh)/h,
                          ddT[2]-T[2]*(g+ddh)/h,
                          ddT[3]-T[3]*(g+ddh)/h,
                          ddT[4]-T[4]*(g+ddh)/h,
                          C[0]*F[0][0]*(g+ddh)/h,
                          C[1]*F[1][0]*(g+ddh)/h,
                          C[2]*F[2][0]*(g+ddh)/h,
                          C[3]*F[3][0]*(g+ddh)/h,
                          C[0]*lam[0]*(g+ddh)/h,
                          C[1]*lam[1]*(g+ddh)/h,
                          C[2]*lam[2]*(g+ddh)/h,
                          C[3]*lam[3]*(g+ddh)/h,
                          ddT[0]-T[0]*(g+ddh)/h,
                          ddT[1]-T[1]*(g+ddh)/h,
                          ddT[2]-T[2]*(g+ddh)/h,
                          ddT[3]-T[3]*(g+ddh)/h,
                          ddT[4]-T[4]*(g+ddh)/h,
                          C[0]*F[0][1]*(g+ddh)/h,
                          C[1]*F[1][1]*(g+ddh)/h,
                          C[2]*F[2][1]*(g+ddh)/h,
                          C[3]*F[3][1]*(g+ddh)/h,
                          C[0]*lam[0]*(g+ddh)/h,
                          C[1]*lam[1]*(g+ddh)/h,
                          C[2]*lam[2]*(g+ddh)/h,
                          C[3]*lam[3]*(g+ddh)/h
                          ])



def kinematic_constraint(t, derivative, flag=False):

    F, C, dF = P.curr_f(t)

    if derivative == 0:
        a_x, a_y, da_x, da_y = P.curr_a(t)

        return array([C[0]*(F[0][0] - a_x[0] - P.f1_x_nom),
                      C[0]*(F[0][1] - a_y[0] - P.f1_y_nom),
                      C[1]*(F[1][0] - a_x[0] - P.f2_x_nom),
                      C[1]*(F[1][1] - a_y[0] - P.f2_y_nom),
                      C[2]*(F[2][0] - a_x[0] - P.f3_x_nom),
                      C[2]*(F[2][1] - a_y[0] - P.f3_y_nom),
                      C[3]*(F[3][0] - a_x[0] - P.f4_x_nom),
                      C[3]*(F[3][1] - a_y[0] - P.f4_y_nom),
                      ])

    elif derivative == 1:
        if flag:
            a_x, a_y, da_x, da_y = P.curr_a(t)

            col1 = array([dF[0][0], da_x[0]])
            col2 = array([dF[0][1], da_y[0]])
            col3 = array([dF[1][0], da_x[0]])
            col4 = array([dF[1][1], da_y[0]])
            col5 = array([dF[2][0], da_x[0]])
            col6 = array([dF[2][1], da_y[0]])
            col7 = array([dF[3][0], da_x[0]])
            col8 = array([dF[3][1], da_y[0]])
            return array([col1, col2, col3, col4, col5, col6, col7, col8])
        else:
            return array([C[0], -C[0], C[0], -C[0], C[1], -C[1], C[1], -C[1],
                          C[2], -C[2], C[2], -C[2], C[3], -C[3], C[3], -C[3]
                          ])


def load_constraint(t, derivative, flag=False):

    lam, dlam = P.curr_lam(t)
    if derivative == 0:
        return array([lam.sum(), lam[0], lam[1], lam[2], lam[3]])

    elif derivative == 1:
        if flag:
            col1 = array([dlam[0], dlam[1], dlam[2], dlam[3]])
            col2 = array([dlam[0]])
            col3 = array([dlam[1]])
            col4 = array([dlam[2]])
            col5 = array([dlam[3]])
            return array([col1, col2, col3, col4, col5])
        else:
            return ones(8)


def eval_g(x, user_data = None):
    """
    Constraint functions:
     - Motion constraint, stitch polynomials together
     - Dynamic constraint, satisfy ZMP
     - Kinematic constraint
     - Load distribution constraint

    :param x: Variable array
    :param user_data:
    :return: Evaluation of all constraints
    """

    P.x = x
    g = zeros(26*(P.np-1), dtype=float_)
    index = 0
    t = 0.0
    for j in range(P.np-1):

        t_bar = 0.0

        # Perform all 4 constraint checks at the beginning of the time period
        stitch = stitching_constraint(t,0)
        g[index:index+len(stitch)] = stitch
        index += len(stitch)

        dyn = dynamic_constraint(t,t_bar,0)
        g[index:index+len(dyn)] = dyn
        index += len(dyn)

        kin = kinematic_constraint(t,0)
        g[index:index+len(kin)] = kin
        index += len(kin)

        lc = load_constraint(t,0)
        g[index:index+len(lc)] = lc
        index += len(lc)

        # Perform 2 constraint checks in the middle of the time period
        t_bar = round(P.tp/2.0,6)
        t = round(t+t_bar,6)

        dyn = dynamic_constraint(t,t_bar,0)
        g[index:index+len(dyn)] = dyn
        index += len(dyn)

        lc = load_constraint(t,0)
        g[index:index+len(lc)] = lc
        index += len(lc)

        t = round(t+t_bar,6)

    return g



def eval_jac_g(x, flag, user_data = None):
    """
    Jacobian of constraint functions
    - "P.Triplet Format for Sparse Matrices" is used per Ipopt requirements

    :param x: Variable array
    :param flag: Determine if position of values or the values themselves are returned
    :param user_data:
    :return: Either the position of all non-zero values in the jacobian or all non-zero
             values in the jacobian.
    """
    P.x = x
    if flag:

        hrow = zeros(106*(P.np-1), dtype=int)
        hcol = zeros(106*(P.np-1), dtype=int)
        row = 0
        index = 0
        t = 0.0
        for j in range(P.np-1):

            t_bar = 0.0

            # Perform all 4 constraint checks at the beginning of the time period
            stitch = stitching_constraint(t,1,flag)
            for constraint in stitch:
                hrow[index:index+len(constraint)] = ones(len(constraint), dtype=int)*row
                hcol[index:index+len(constraint)] = constraint
                index += len(constraint)
                row += 1

            dyn = dynamic_constraint(t,t_bar,1,flag)
            for constraint in dyn:
                hrow[index:index+len(constraint)] = ones(len(constraint), dtype=int)*row
                hcol[index:index+len(constraint)] = constraint
                index += len(constraint)
                row += 1

            kin = kinematic_constraint(t,1,flag)
            for constraint in kin:
                hrow[index:index+len(constraint)] = ones(len(constraint), dtype=int)*row
                hcol[index:index+len(constraint)] = constraint
                index += len(constraint)
                row += 1

            lc = load_constraint(t,1,flag)
            for constraint in lc:
                hrow[index:index+len(constraint)] = ones(len(constraint), dtype=int)*row
                hcol[index:index+len(constraint)] = constraint
                index += len(constraint)
                row += 1

            # Perform 2 constraint checks in the middle of the time period
            t_bar = round(P.tp/2.0,6)
            t = round(t+t_bar,6)

            dyn = dynamic_constraint(t,t_bar,1,flag)
            for constraint in dyn:
                hrow[index:index+len(constraint)] = ones(len(constraint), dtype=int)*row
                hcol[index:index+len(constraint)] = constraint
                index += len(constraint)
                row += 1

            lc = load_constraint(t,1,flag)
            for constraint in lc:
                hrow[index:index+len(constraint)] = ones(len(constraint), dtype=int)*row
                hcol[index:index+len(constraint)] = constraint
                index += len(constraint)
                row += 1

            t = round(t+t_bar,6)

        return (array(hrow), array(hcol))

    else:
        values = zeros(106*(P.np-1), dtype=float_)

        index = 0
        t = 0.0
        for j in range(P.np-1):

            t_bar = 0.0

            # Perform all 4 constraint checks at the beginning of the time period
            stitch = stitching_constraint(t,1)
            values[index:index+len(stitch)] = stitch
            index += len(stitch)

            dyn = dynamic_constraint(t,t_bar,1)
            values[index:index+len(dyn)] = dyn
            index += len(dyn)

            kin = kinematic_constraint(t,1)
            values[index:index+len(kin)] = kin
            index += len(kin)

            lc = load_constraint(t,1)
            values[index:index+len(lc)] = lc
            index += len(lc)

            # Perform 2 constraint checks in the middle of the time period
            t_bar = round(P.tp/2.0,6)
            t = round(t+t_bar,6)

            dyn = dynamic_constraint(t,t_bar,1)
            values[index:index+len(dyn)] = dyn
            index += len(dyn)

            lc = load_constraint(t,1)
            values[index:index+len(lc)] = lc
            index += len(lc)

            t = round(t+t_bar,6)

        return values

