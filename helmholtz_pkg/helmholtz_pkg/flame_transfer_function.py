"""
    flame transfer functions:

    n-tau model,
    state-space representation
"""

import numpy as np
import matplotlib.pyplot as plt


def n_tau(N3, tau):
    """
    :param N3: non-dimensional interaction index
    :param tau: time delay [s]
    :return: function
    """
    def inner_func(omega, k=0):
        return N3 * (1j * tau)**k * np.exp(1j * omega * tau)
    return inner_func


def state_space(A, b, c, d):
    """
    vectfit3.m is written using the expansion e^{i omega t},
    while the Helmholtz solver is written using e^{- i omega t}.
    omega in vectfit3.m is the complex conjugate of omega in the Helmholtz solver
    and the same goes for the flame transfer function.

    :param A: (N, N) array
    :param b: (N, 1) array
    :param c: (1, N) array
    :param d: (1, 1) array
    :return: function
    """
    Id = np.eye(*A.shape)

    def inner_func(omega, k=0):
        """
        :param omega: complex angular frequency [rad/s]
        :param k: k-th order derivative
        :return: transfer function/frequency response
        """
        omega = np.conj(omega)
        Mat = (- 1j) ** k * np.math.factorial(k) * \
            np.linalg.matrix_power(1j * omega * Id - A, - (k + 1))
        row = np.dot(c, Mat)
        H = np.dot(row, b)
        if k == 0:
            H += d
        return np.conj(H[0][0])
    return inner_func


if __name__ == '__main__':

    from scipy.io import loadmat

    N3_ = 1.
    tau_ = 0.003

    mat = loadmat('/Users/stefanofalco/VectorFitting/ftf/ftf.mat')

    freq = mat['freq'][0]
    f = mat['f'][0]
    fit = mat['fit'][0]

    A_ = mat['A']
    b_ = mat['b']
    c_ = mat['c']
    d_ = mat['d']


    def test_1():

        # ftf1 = n_tau(N3_, tau_)

        ftf2 = state_space(A_, b_, c_, d_)

        f_ = np.linspace(0, 1060, 1061)
        omega_ = 2 * np.pi * f_

        # z1 = ftf1(omega_)

        z2 = np.zeros_like(omega_, dtype=np.complex128)
        for i in range(1061):
            z2[i] = ftf2(omega_[i])

        plt.plot(freq, np.abs(f), 's', f_, np.abs(z2))
        plt.show()
        plt.plot(freq, np.angle(f), 's', f_, np.angle(z2))
        plt.show()


    def test_2():

        ftf2 = state_space(A_, b_, c_, d_)

        f_ = np.linspace(0, 1060, 1061)
        omega_ = 2 * np.pi * f_

        z = np.zeros(1061, dtype=np.complex128)
        for i in range(1061):
            z[i] = ftf2(omega_[i])

        dz1 = np.zeros(1060, dtype=np.complex128)
        for i in range(1060):
            dz1[i] = (z[i+1] - z[i])/(omega_[i+1] - omega_[i])

        dz2 = np.zeros(1061, dtype=np.complex128)
        for i in range(1061):
            dz2[i] = ftf2(omega_[i], k=1)

        plt.plot(omega_[:-1], np.abs(dz1))
        plt.plot(omega_, np.abs(dz2), '--')
        plt.show()

        plt.plot(omega_[:-1], np.angle(dz1))
        plt.plot(omega_, np.angle(dz2), '--')
        plt.show()


    test_2()
