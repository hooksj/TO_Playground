#!/usr/bin/env python
"""
Contains all function to use the pyiopt libraries.
 - Needs to have the trajectory.structure class passed in
"""

__author__  	= "Joshua Hooks"
__email__   	= "hooksjrose@gmail.com"
__copyright__ 	= "Copyright 2017 RoMeLa"
__date__ 		= "Aug. 29th, 2018"

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

def solve (p, warm_start):
    global P
    P = p
    t0 = time.time()



    # Create the problem for Ipopt to solve
    P.nlp = pyipopt.create(P.nvar, P.x_L, P.x_U, P.ncon, P.g_L, P.g_U, P.nnjz, P.nnzh, eval_f, eval_grad_f, eval_g, eval_jac_g)

    # Set Ipopt solve options
    """
    P.nlp.int_option("print_level",1)
    P.nlp.num_option("acceptable_tol",10.0)
    P.nlp.num_option("tol",10.0)
    P.nlp.num_option("dual_inf_tol",1000.0)
    P.nlp.num_option("compl_inf_tol",1000.0)
    P.nlp.int_option("max_iter", 10000)
    P.nlp.str_option("linear_solver", "mumps")
    """

    P.nlp.str_option("linear_solver", "mumps")
    P.nlp.str_option("jacobian_approximation", "exact")
    P.nlp.str_option("hessian_approximation", "limited-memory")
    P.nlp.num_option("max_cpu_time", 40.0)
    P.nlp.num_option("tol", 10.0)
    P.nlp.str_option("print_timing_statistics", "no")
    P.nlp.str_option("print_user_options", "no")
    P.nlp.int_option("print_level",4)

    if warm_start:
        P.nlp.str_option('mu_strategy', 'adaptive')
        P.nlp.str_option("warm_start_init_point", "yes")
    #P.nlp.str_option("derivative_test","first-order")

    print "Calling solve"

    x, zl, zu, constrain_mu, obj, status = P.nlp.solve(P.x) # Solve the NLP problem with Ipopt
    #nlp.close()
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

    P.x = x
    f = 0.0

    for t in P.time_vector:
        t = round(t,6)
        _, C, _ = P.curr_f(t)
        a_x, a_y, da_x, da_y = P.curr_a(t)
        lam, dlam = P.curr_lam(t)
        tb = 0.0

        if P.ACost.on:
            f += P.ACost.evaluate(tb, a_x, a_y)

        if P.LCost.on:
            f += P.LCost.evaluate(lam, C)

        _, C, _ = P.curr_f(round(t+P.tp/2.0,6))
        a_x, a_y, da_x, da_y = P.curr_a(round(t+P.tp/2.0,6))
        lam, dlam = P.curr_lam(round(t+P.tp/2.0,6))
        tb = round(P.tp/2.0,6)

        if P.ACost.on:
            f += P.ACost.evaluate(tb, a_x, a_y)

        if P.LCost.on:
            f += P.LCost.evaluate(lam, C)

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

    for t in P.time_vector:
        t = round(t,6)
        _, C, _ = P.curr_f(t)
        a_x, a_y, da_x, da_y = P.curr_a(t)
        lam, dlam = P.curr_lam(t)
        tb = 0.0

        if P.ACost.on:
            index, values = P.ACost.derivative(tb, a_x, a_y, da_x, da_y)
            grad_f[index] += values

        if P.LCost.on:
            index, values = P.LCost.derivative(lam, dlam, C)
            grad_f[index] += values

        _, C, _ = P.curr_f(round(t+P.tp/2.0,6))
        a_x, a_y, da_x, da_y = P.curr_a(round(t+P.tp/2.0,6))
        lam, dlam = P.curr_lam(round(t+P.tp/2.0,6))
        tb = round(P.tp/2.0,6)

        if P.ACost.on:
            index, values = P.ACost.derivative(tb, a_x, a_y, da_x, da_y)
            grad_f[index] += values

        if P.LCost.on:
            index, values = P.LCost.derivative(lam, dlam, C)
            grad_f[index] += values

    return grad_f


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
    g = zeros(P.ncon, dtype=float_)
    index = 0

    for t in P.time_vector:
        t = round(t,6)
        F, C, _ = P.curr_f(t)
        a_x, a_y, da_x, da_y = P.curr_a(t)
        a_x1, a_y1, da_x1, da_y1 = P.curr_a(round(t+P.tp,6))
        lam, dlam = P.curr_lam(t)
        h, ddh = P.curr_h(t)
        t_bar = 0.0

        if P.SCon.on:
            g[index:index+P.SCon.ncon] = P.SCon.evaluate(a_x, a_y, a_x1, a_y1)
            index += P.SCon.ncon

        if P.DCon.on:
            g[index:index+P.DCon.ncon] = P.DCon.evaluate(t_bar, a_x, a_y, F, C, lam, h, ddh)
            index += P.DCon.ncon

        if P.KCon.on:
            g[index:index+P.KCon.ncon] = P.KCon.evaluate(a_x, a_y, F, C)
            index += P.KCon.ncon

        if P.LCon.on:
            g[index:index+P.LCon.ncon] = P.LCon.evaluate(lam)
            index += P.LCon.ncon

        t_bar = round(P.tp/2.0,6)
        a_x, a_y, da_x, da_y = P.curr_a(round(t+t_bar,6))
        F, C, _ = P.curr_f(round(t+t_bar,6))
        lam, dlam = P.curr_lam(round(t+t_bar,6))
        h, ddh = P.curr_h(round(t+t_bar,6))

        if P.DCon.on:
            g[index:index+P.DCon.ncon] = P.DCon.evaluate(t_bar, a_x, a_y, F, C, lam, h, ddh)
            index += P.DCon.ncon

        if P.LCon.on:
            g[index:index+P.LCon.ncon] = P.LCon.evaluate(lam)
            index += P.LCon.ncon

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

        hrow = zeros(P.nnjz, dtype=int)
        hcol = zeros(P.nnjz, dtype=int)
        row = 0
        index = 0

        for t in P.time_vector:
            t = round(t,6)
            a_x, a_y, da_x, da_y = P.curr_a(t)
            a_x1, a_y1, da_x1, da_y1 = P.curr_a(round(t+P.tp,6))
            lam, dlam = P.curr_lam(t)
            _, _, dF = P.curr_f(t)

            if P.SCon.on:
                stitch = P.SCon.derivative_structure(da_x, da_y, da_x1, da_y1)
                for constraint in stitch:
                    hrow[index:index+len(constraint)] = ones(len(constraint), dtype=int)*row
                    hcol[index:index+len(constraint)] = constraint
                    index += len(constraint)
                    row += 1

            if P.DCon.on:
                dyn = P.DCon.derivative_structure(da_x, da_y, dlam, dF)
                for constraint in dyn:
                    hrow[index:index+len(constraint)] = ones(len(constraint), dtype=int)*row
                    hcol[index:index+len(constraint)] = constraint
                    index += len(constraint)
                    row += 1

            if P.KCon.on:
                kin = P.KCon.derivative_structure(da_x, da_y, dF)
                for constraint in kin:
                    hrow[index:index+len(constraint)] = ones(len(constraint), dtype=int)*row
                    hcol[index:index+len(constraint)] = constraint
                    index += len(constraint)
                    row += 1

            if P.LCon.on:
                lc = P.LCon.derivative_structure(dlam)
                for constraint in lc:
                    hrow[index:index+len(constraint)] = ones(len(constraint), dtype=int)*row
                    hcol[index:index+len(constraint)] = constraint
                    index += len(constraint)
                    row += 1

            # Perform 2 constraint checks in the middle of the time period
            t_bar = round(P.tp/2.0,6)
            a_x, a_y, da_x, da_y = P.curr_a(round(t+t_bar,6))
            lam, dlam = P.curr_lam(round(t+t_bar,6))
            _, _, dF = P.curr_f(round(t+t_bar,6))

            if P.DCon.on:
                dyn = P.DCon.derivative_structure(da_x, da_y, dlam, dF)
                for constraint in dyn:
                    hrow[index:index+len(constraint)] = ones(len(constraint), dtype=int)*row
                    hcol[index:index+len(constraint)] = constraint
                    index += len(constraint)
                    row += 1

            if P.LCon.on:
                lc = P.LCon.derivative_structure(dlam)
                for constraint in lc:
                    hrow[index:index+len(constraint)] = ones(len(constraint), dtype=int)*row
                    hcol[index:index+len(constraint)] = constraint
                    index += len(constraint)
                    row += 1

        return (array(hrow), array(hcol))

    else:
        values = zeros(P.nnjz+100, dtype=float_)

        index = 0
        for t in P.time_vector:
            t = round(t,4)

            t_bar = 0.0
            F, C, dF = P.curr_f(t)
            lam, dlam = P.curr_lam(t)
            h, ddh = P.curr_h(t)

            if P.SCon.on:
                values[index:index+P.SCon.nnjz] = P.SCon.derivative()
                index += P.SCon.nnjz

            if P.DCon.on:
                dcon = P.DCon.derivative(t_bar, F, dF, C, lam, h, ddh)
                values[index:index+len(dcon)] = dcon
                index += len(dcon)

            if P.KCon.on:
                kcon = P.KCon.derivative(C, dF)
                values[index:index+len(kcon)] = kcon
                index += len(kcon)

            if P.LCon.on:
                values[index:index+P.LCon.nnjz] = P.LCon.derivative()
                index += P.LCon.nnjz

            # Perform 2 constraint checks in the middle of the time period
            t_bar = round(P.tp/2.0,4)
            F, C, dF = P.curr_f(round(t+t_bar,4))
            lam, dlam = P.curr_lam(round(t+t_bar,4))
            h, ddh = P.curr_h(round(t+t_bar,4))

            if P.DCon.on:
                dcon = P.DCon.derivative(t_bar, F, dF, C, lam, h, ddh)
                values[index:index+len(dcon)] = dcon
                index += len(dcon)

            if P.LCon.on:
                values[index:index+P.LCon.nnjz] = P.LCon.derivative()
                index += P.LCon.nnjz

        #pdb.set_trace()

        return values

