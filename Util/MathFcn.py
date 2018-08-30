#!usr/bin/env python
__author__ = "Min Sung Ahn"
__email__ = "aminsung@gmail.com"
__copyright__ = "Copyright 2016 RoMeLa"
__date__ = "June 25, 2016"

__version__ = "0.1.0"
__status__ = "Production"

import sys

import re

import numpy as np
import scipy.linalg
import pdb

import Util.NestedDict as NestedDict

def wrap_between_pi_and_neg_pi(angle):
    return (angle + np.pi) % (2*np.pi) - np.pi

# ===== Awesome math functions
def rodrigues(a_vec, th):
    # Kajita pg. 35
    # print type(a_vec)
    try:
        if a_vec.shape != (3, 1):
            raise IncorrectMatrixSize
        hat_mat = create_hat_mat(a_vec)
        # return np.eye(3) + hat_mat * np.sin(th * np.pi / 180) + (np.dot(hat_mat, hat_mat)) * (1 - np.cos(th * np.pi / 180.0))
        return np.eye(3) + hat_mat * np.sin(th) + (np.dot(hat_mat, hat_mat)) * (1 - np.cos(th))
    except IncorrectMatrixSize:
        raise IncorrectMatrixSize("Your angular velocity vector size is incorrect!")

def create_hat_mat(vec):
    try:
        if vec.shape != (3, 1):
            raise IncorrectMatrixSize
        return np.array([ [0.0, -vec[2,0], vec[1,0]], [vec[2,0], 0.0, -vec[0,0]], [-vec[1,0], vec[0,0], 0.0] ])
    except IncorrectMatrixSize:
        raise IncorrectMatrixSize("Vector you are trying to turn into a hat matrix may not be the same size.")

def roll_rot_mat(th):
    s = np.sin(th)
    c = np.cos(th)
    rot = np.array([ [1.0, 0.0, 0.0], [0.0, c, -s], [0, s, c]])
    return rot

def pitch_rot_mat(th):
    s = np.sin(th)
    c = np.cos(th)
    rot = np.array([ [c, 0.0, s], [0.0, 1.0, 0.0], [-s, 0.0, c]])
    return rot

def yaw_rot_mat(th):
    s = np.sin(th)
    c = np.cos(th)
    rot = np.array([ [c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])
    return rot

def auto_eval_spline(coef, order, t):
    """
    Automatically evaluate a spline by providing the coefficients, the order, and the time to eval at.

    :param coef: Descending list of coefficients of the spline
    :param order: Highest order of the spline
    :param t: Time value to evaulate the spline at
    :return: Evaluated spline value
    """
    val = 0
    for idx in reversed(range(order+1)):
        val += coef[order-idx]*(t**idx)
    return val

def mech_zero(arr, tol=1e-6):
    """
    Machine zero numpy arrays when they don't feel like doing it themselves.

    :param arr: numpy array with potentially very small numbers
    :return arr: Clean array with actual zeros
    """
    arr[np.abs(arr) < tol] = 0
    return arr

# ===== Awesome control functions
def dlqr(A, B, Q, R, N=0.0):
    """
    Cheap MATLAB port of linear quadratic regulator design for discrete time.

    [K, S, E] = dlqr(A, B, Q, R, N)

    K: Control gain
    S: Solution to the Algebraic Riccati Eqn.
    E: Closed-loop eigenvalues
    """

    S = scipy.linalg.solve_discrete_are(A, B, Q, R)
    K = scipy.linalg.inv(R + B.T.dot(S).dot(B)).dot(B.T).dot(S).dot(A)
    E = scipy.linalg.eigvals(A-B.dot(K))
    return K, S, E


# ===== Awesome quaternion math functions
"""
    All quaternion math is done in accordance with:
    "Indirect Kalman Filter for 3D Attitude Estimation: A Tutorial for Quaternion Algebra"
"""
def skew(q):
    """
    Creates a matrix for quaternion multiplication

    :param q: quaternion in the form q = [x, y, z, w]
    :return: 3x3 matrix
    """
    skew = np.array([[0.0,   -q[2][0],   q[1][0]],
                     [q[2][0],   0.0,     -q[0][0]],
                     [-q[1][0],  q[0][0],    0.0]])
    return skew

def quat2rotm(Q):
    """
    Converts a quaternion into a rotation matrix

    :param Q: quaternion in the form Q = [x, y, z, w]
    :return: 3x3 rotation matrix
    """
    q = Q[0:3]
    q4 = Q[3]
    C = (2.0*q4*q4 - 1.0)*np.eye(3) - 2.0*q4*skew(q) + 2.0*np.outer(q,q)
    return C

def rotm2quat(C):
    """
    Converts a 3x3 rotation matrix into a 4x1 quaternion

    :param C: 3x3 rotation matrix
    :return: 4x1 quaternion
    """
    q4 = np.sqrt(1 + np.trace(C)) / 2.0
    q1 = (C[1,2] - C[2,1]) / (4.0*q4)
    q2 = (C[2,0] - C[0,2]) / (4.0*q4)
    q3 = (C[0,1] - C[1,0]) / (4.0*q4)
    q = np.array([[q1], [q2], [q3], [q4]])
    return q

def eul2rotm(rX, rY, rZ): #
    """
    Converts Euler angles to a 3x3 rotation matrix using the ZYX convention

    :param rX: rotation about the X axis (radians)
    :param rY: rotation about the Y axis (radians)
    :param rZ: rotation about the Z axis (radians)
    :return: 3x3 rotation matrix
    """
    cA = np.cos(rZ)
    sA = np.sin(rZ)
    cG = np.cos(rX)
    sG = np.sin(rX)
    cB = np.cos(rY)
    sB = np.sin(rY)

    C = np.array([[cA*cB, cA*sB*sG-sA*cG, cA*sB*cG+sA*sG],
                 [sA*cB, sA*sB*sG+cA*cG, sA*sB*cG-cA*sG],
                 [-sB,       cB*sG,          cB*cG]])
    return C

def fix2rotm(rX, rY, rZ):
    """
    Converts Fixed angles to a 3x3 rotation matrix using the XYZ convention

    :param rX: rotation about the X axis (radians)
    :param rY: rotation about the Y axis (radians)
    :param rZ: rotation about the Z axis (radians)
    :return: 3x3 rotation matrix
    """
    cA = np.cos(rX)
    sA = np.sin(rX)
    cG = np.cos(rZ)
    sG = np.sin(rZ)
    cB = np.cos(rY)
    sB = np.sin(rY)

    C = np.array([[cA*cB, cA*sB*sG-sA*cG, cA*sB*cG+sA*sG],
                 [sA*cB, sA*sB*sG+cA*cG, sA*sB*cG-cA*sG],
                 [-sB,       cB*sG,          cB*cG]])
    return C

def rotm2fix(C):
    """
    Converts a 3x3 rotation matrix to Fixed angles

    :param: 3x3 rotation matrix
    :return rX: rotation about the X axis (radians)
    :return rY: rotation about the Y axis (radians)
    :return rZ: rotation about the Z axis (radians)
    """
    rZ = np.arctan2(C[2,1], C[2,2])
    rY = np.arctan2(-C[2,0],np.sqrt(C[2,1]**2.0 + C[2,2]**2.0))
    rX = np.arctan2(C[1,0], C[0,0])
    return rX, rY, rZ

def rotm2eul(C):
    """
    Converts a 3x3 rotation matrix to Euler angles

    :param: 3x3 rotation matrix
    :return rX: rotation about the X axis (radians)
    :return rY: rotation about the Y axis (radians)
    :return rZ: rotation about the Z axis (radians)
    """
    rX = np.arctan2(C[2,1], C[2,2])
    rY = np.arctan2(-C[2,0],np.sqrt(C[2,1]**2.0 + C[2,2]**2.0))
    rZ = np.arctan2(C[1,0], C[0,0])
    return rX, rY, rZ

def eul2quat(rX, rY, rZ):
    """
    Converts euler angles into a 4x1 quaternion

    :param rX: rotation about the X axis (radians)
    :param rY: rotation about the Y axis (radians)
    :param rZ: rotation about the Z axis (radians)
    :return: 4x1 quaternion q = [x, y, z, w]
    """
    return  rotm2quat(eul2rotm(rX,rY,rZ))

def quat2eul(q):
    """
    Converts a 4x1 quaternion into Euler angles

    :param q: 4x1 quaternion q = [x, y, z, w]
    :return rX: rotation about the X axis (radians)
    :return rY: rotation about the Y axis (radians)
    :return rZ: rotation about the Z axis (radians)
    """
    return rotm2eul(quat2rotm(q))

def fix2quat(rX, rY, rZ):
    """
    Converts fixed angles into a 4x1 quaternion

    :param rX: rotation about the X axis (radians)
    :param rY: rotation about the Y axis (radians)
    :param rZ: rotation about the Z axis (radians)
    :return: 4x1 quaternion q = [x, y, z, w]
    """
    return  rotm2quat(fix2rotm(rX,rY,rZ))

def quat2fix(q):
    """
    Converts a 4x1 quaternion into Fixed angles

    :param q: 4x1 quaternion q = [x, y, z, w]
    :return rX: rotation about the X axis (radians)
    :return rY: rotation about the Y axis (radians)
    :return rZ: rotation about the Z axis (radians)
    """
    return rotm2fix(quat2rotm(q))

def quat_cross(q, p):
    Q = np.array([[q[3][0], q[2][0],  -q[1][0], q[0][0]],
                 [-q[2][0], q[3][0],  q[0][0],  q[1][0]],
                 [q[1][0],  -q[0][0], q[3][0],  q[2][0]],
                 [-q[0][0], -q[1][0], -q[2][0], q[3][0]]])
    return Q.dot(p)

def quat_inv(q):
    q_inv = np.array([[-q[0][0]], [-q[1][0]], [-q[2][0]], [q[3][0]]])

    return q_inv


class IncorrectMatrixSize(Exception):
    pass


if __name__ == "__main__":
    rx = 1.253
    ry = 1.456
    rz = 1.563
    q = eul2quat(rx, ry, rz)
    q_inv = quat_inv(q)
    test = quat_cross(q, q_inv)
    print test
