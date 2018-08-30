__author__ = "Daniel Sun"
__email__ = "sundaniel3@gmail.com"
__copyright__ = "Copyright 2017 RoMeLa"
__date__ = "December 27 2017"

__version__ = "0.0.1"
__status__ = "Prototype"


import numpy.matlib as np
import math

"""
Matrix functions and the Transform class are located here

"""


class Transform(np.matrix):
    """
    This is a simple class that inherits from numpy's matrix class.
    You can use the regular operators on it but it also has special functions to let you access things
    and create transforms easily.

    >>> T = Transform() # basic transform
    >>> print Transform()
    [[ 1.  0.  0.  0.]
     [ 0.  1.  0.  0.]
     [ 0.  0.  1.  0.]
     [ 0.  0.  0.  1.]]

    >>> R = Rpitch(3.14159/2)
    >>> p = np.matrix('3;4;2')
    >>> T = Transform(R,p) # creates transform with rotation matrices R and p
    >>> print T
    [[  1.32679490e-06   0.00000000e+00   1.00000000e+00   3.00000000e+00]
     [  0.00000000e+00   1.00000000e+00   0.00000000e+00   4.00000000e+00]
     [ -1.00000000e+00   0.00000000e+00   1.32679490e-06   2.00000000e+00]
     [  0.00000000e+00   0.00000000e+00   0.00000000e+00   1.00000000e+00]]

    >>> print T.p # gets the position vector, 3x1 matrix
    [[ 3.]
     [ 4.]
     [ 2.]]
    >>> print T.P # gets the extended position vector, 4x1 matrix
    [[ 3.]
     [ 4.]
     [ 2.]
     [ 1.]]
    >>> print T.R # gets the 3x3 rotation matrix
    [[  1.32679490e-06   0.00000000e+00   1.00000000e+00]
     [  0.00000000e+00   1.00000000e+00   0.00000000e+00]
     [ -1.00000000e+00   0.00000000e+00   1.32679490e-06]]

    >>> print T.I # gets the inverse of the matrix using R.T and R.T*(-T.p)
    [[  1.32679490e-06   0.00000000e+00  -1.00000000e+00   1.99999602e+00]
     [  0.00000000e+00   1.00000000e+00   0.00000000e+00  -4.00000000e+00]
     [  1.00000000e+00   0.00000000e+00   1.32679490e-06  -3.00000265e+00]
     [  0.00000000e+00   0.00000000e+00   0.00000000e+00   1.00000000e+00]]
    >>> print T.rotate_y(3.1415/2) # applies an Euler pitch rotation
    [[ -9.99999999e-01   0.00000000e+00   4.76535898e-05   3.00000000e+00]
     [  0.00000000e+00   1.00000000e+00   0.00000000e+00   4.00000000e+00]
     [ -4.76535898e-05   0.00000000e+00  -9.99999999e-01   2.00000000e+00]
     [  0.00000000e+00   0.00000000e+00   0.00000000e+00   1.00000000e+00]]

    >>> print T.translate_z(1.3) # applies a vertical translation
    [[  1.32679490e-06   0.00000000e+00   1.00000000e+00   4.30000000e+00]
     [  0.00000000e+00   1.00000000e+00   0.00000000e+00   4.00000000e+00]
     [ -1.00000000e+00   0.00000000e+00   1.32679490e-06   2.00000172e+00]
     [  0.00000000e+00   0.00000000e+00   0.00000000e+00   1.00000000e+00]]

    ### NOTE ###
    Due to the way that python handles operators, right multiplying matrix by a transform results in a matrix.
    You can always view a 4x4 matrix as a Transform again by using the view method.

    M*T # creates matrix
    M*T.view(Transform) # gives view of transform
    """

    def __new__(cls, R=np.eye(3), p=np.matrix([0, 0, 0]).T):
        if R is None and p is None:
            obj = super(Transform, cls).__new__(cls, np.eye(4))
        else:
            obj = super(Transform, cls).__new__(cls, np.block([[R, p], [0, 0, 0, 1]]))
            # obj = np.block([[R, p], [0, 0, 0, 1]]) appears to be about the same
        return obj.view(cls)

    def rotate_x(self, q):
        return self * Transform(Rroll(q))

    def rotate_y(self, q):
        return self * Transform(Rpitch(q))

    def rotate_z(self, q):
        return self * Transform(Ryaw(q))

    def translate_x(self, d):
        return self * Transform(p=np.matrix([d, 0, 0]).T)

    def translate_y(self, d):
        return self * Transform(p=np.matrix([0, d, 0]).T)

    def translate_z(self, d):
        return self * Transform(p=np.matrix([0, 0, d]).T)

    def mDH(self, alpha, a, theta, d):
        return self.rotate_x(alpha).translate_x(a).rotate_z(theta).translate_z(d)

    @property
    def p(self):
        return self[:3, 3]

    @property
    def P(self):
        return self[:, 3]

    @property
    def R(self):
        return self[:3, :3]

    @property
    def I(self):
        return Transform(self.R.T, self.R.T * -self.p)


def zyz_euler_angle(a, b, g):
    """
    Creates a rotation matrix based on z-y-z rotation angles

    :param a: 1st z rotation, in radians
    :param b: 2nd y rotation, in radians
    :param g: 3rd z rotation, in radians
    :return: 3x3 rotation matrix
    """
    return np.matrix([[math.cos(a) * math.cos(b) * math.cos(g) - math.sin(a) * math.sin(g),
                       -math.cos(a) * math.cos(b) * math.sin(g) - math.sin(a) * math.cos(g), math.cos(a) * math.sin(b)],
                      [math.sin(a) * math.cos(b) * math.cos(g) + math.cos(a) * math.sin(g),
                       -math.sin(a) * math.cos(b) * math.sin(g) + math.cos(a) * math.cos(g), math.sin(a) * math.sin(b)],
                      [-math.sin(b) * math.cos(g), math.sin(b) * math.sin(g), math.cos(b)]])


def inverse_zyz_euler_angle(R, flip_beta=True):
    """
    Given a rotation matrix R, find the corresponding ZYZ Euler angle decomposition
    :param R: 3x3 rotation matrix
    :return: (yaw, math.pitch, yaw) rotation tuple
    """

    if flip_beta:
        beta = math.atan2(-math.sqrt(R[2, 0] ** 2 + R[2, 1] ** 2), R[2, 2])
    else:
        beta = math.atan2(-math.sqrt(R[2, 0] ** 2 + R[2, 1] ** 2), R[2, 2])
    if beta == 0 or beta == 2 * math.pi:
        alpha = 0
        gamma = math.atan2(-R[0, 1], R[0, 0])

    elif beta == math.pi or beta == -math.pi:
        alpha = 0
        gamma = math.atan2(R[0, 1], -R[0, 0])

    else:
        alpha = math.atan2(R[1, 2] / math.sin(beta), R[0, 2] / math.sin(beta))
        gamma = math.atan2(R[2, 1] / math.sin(beta), -R[2, 0] / math.sin(beta))
    return (alpha, beta, gamma)


def Rroll(q):
    return np.matrix([[1, 0, 0], [0, math.cos(q), -math.sin(q)], [0, math.sin(q), math.cos(q)]])


def Rpitch(q):
    return np.matrix([[math.cos(q), 0, math.sin(q)], [0, 1, 0], [-math.sin(q), 0, math.cos(q)]])


def Ryaw(q):
    return np.matrix([[math.cos(q), -math.sin(q), 0], [math.sin(q), math.cos(q), 0], [0, 0, 1]])


def interpolate_matrix(Goal, Curr, nTimes):
    nTimes = int(nTimes)
    dT = np.true_divide((Goal - Curr), nTimes)

    return [Curr + dT * i for i in xrange(nTimes + 1)]


if __name__ == "__main__":
    import unittest
    from numpy.testing import assert_allclose
    from math import pi


    class TestTransform(unittest.TestCase):
        def test_basic(self):
            # test basic construction
            T = Transform()
            assert_allclose(T, np.matrix(np.eye(4)))

            # test initializing with rotation matrix
            T = Transform(Rpitch(pi / 4))
            N = np.eye(4)
            N[:3, :3] = Rpitch(pi / 4)
            assert_allclose(T, N)

            # test initializing with position vector
            pos = np.matrix([5, 2, 4]).T
            T = Transform(p=pos)
            N = np.eye(4)
            N[:3, 3] = pos
            assert_allclose(T, N)

            # test initializing with both rotation and position vector
            R = Ryaw(5)
            p = np.matrix('3;2;4')
            T = Transform(R, p)
            N = np.eye(4)
            N[:3, :3] = R
            N[:3, 3] = p
            assert_allclose(T, N)

        def test_offsets(self):
            # test rotation
            p = np.matrix('1;0;0')
            T = Transform(p=p)
            assert_allclose(T.rotate_x(pi / 2).rotate_z(pi / 2).rotate_y(pi / 2).rotate_z(-pi / 2), T, atol=1e-7)

            # test translation
            T = Transform()
            assert_allclose(T.translate_x(5).translate_y(-2).translate_z(3).p, np.mat('5;-2;3'))

            # test mDH
            T = Transform().mDH(pi, 1, pi / 4, 2).mDH(0, 1, 0, 0)
            N = np.eye(4)
            N[:3, :3] = Rroll(pi) * Ryaw(pi / 4)
            N[:3, 3] = np.mat('1.7071;-.7071;-2')
            assert_allclose(T, N, atol=1e-4)

            # test inverse
            assert_allclose(T * T.I, Transform(), atol=1e-7)

        def test_advanced(self):
            T = Transform()
            # right multiplying matrix by a transform results in a matrix
            self.assertEqual(type((np.eye(4) * T)), np.matrix)
            # right multiplying transform by a matrix results in a transform
            self.assertEqual(type(T * np.eye(4)), Transform)
            # right multiplying transform by a non 4x4 matrix results in a matrix
            p = np.matrix('1;3;2')
            print type(T * p)
            self.assertEqual(type(T * p), type(p))


    class TestZYZEuler(unittest.TestCase):
        def testzyz_fk_ik(self):
            R = np.matrix([[1., -0., 0.],
                           [0., 1., 0.],
                           [-0., 0., 1.]])
            self.assertTrue(np.allclose(zyz_euler_angle(0, 0, 0), R))
            self.assertTrue(np.allclose(zyz_euler_angle(*inverse_zyz_euler_angle(R)), R))

            R = np.matrix([[0.70710678, -0.70710678, 0.],
                           [0.70710678, 0.70710678, 0.],
                           [-0., 0., 1.]])
            self.assertTrue(np.allclose(zyz_euler_angle(math.pi / 4, 0, 0), R))
            self.assertTrue(np.allclose(zyz_euler_angle(*inverse_zyz_euler_angle(R)), R))

            R = np.matrix([[0.5, -0., 0.8660254],
                           [0., 1., 0.],
                           [-0.8660254, 0., 0.5]])
            self.assertTrue(np.allclose(zyz_euler_angle(0, math.pi / 3, 0), R))
            self.assertTrue(np.allclose(zyz_euler_angle(*inverse_zyz_euler_angle(R)), R))

            R = np.matrix([[0.5, -0.5, 0.70710678],
                           [0.70710678, 0.70710678, 0.],
                           [-0.5, 0.5, 0.70710678]])
            self.assertTrue(np.allclose(zyz_euler_angle(0, math.pi / 4, math.pi / 4), R))
            self.assertTrue(np.allclose(zyz_euler_angle(*inverse_zyz_euler_angle(R)), R))

            R = np.matrix([[-0.55034481, -0.70029646, 0.45464871],
                           [0.70029646, -0.09064712, 0.70807342],
                           [-0.45464871, 0.70807342, 0.54030231]])
            self.assertTrue(np.allclose(zyz_euler_angle(1, 1, 1), R))
            self.assertTrue(np.allclose(zyz_euler_angle(*inverse_zyz_euler_angle(R)), R))


    unittest.main()

