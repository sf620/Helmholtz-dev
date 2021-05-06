from petsc4py import PETSc

import numpy as np

from helmholtz_pkg.petsc4py_utils import mult_complex_scalar_real_matrix as mult_C
from helmholtz_pkg.petsc4py_utils import mult_complex_scalar_complex_matrix as mult_B
from helmholtz_pkg.petsc4py_utils import vector_matrix_vector

from helmholtz_pkg.eigensolvers import eps_solver
from helmholtz_pkg import eigenvectors


def newton(operators, D,
           init, nev=2, i=0,
           tol=1e-8, maxiter=50,
           print_results=False):
    """
    The convergence strongly depends/relies on the initial value assigned to omega.
    Targeting zero in the shift-and-invert (spectral) transformation or, more in general,
    seeking for the eigenvalues nearest to zero might also be problematic.
    The implementation uses the TwoSided option to compute the adjoint eigenvector
    (IT HAS BEEN TESTED).

    helmholtz_pkg.eigensolvers and eigenvectors might have been changed in the meantime...

    :param operators:
    :param D:
    :param init: initial value assigned to omega
    :param nev:
    :param i:
    :param tol:
    :param maxiter:
    :param print_results:
    :return:
    """

    A = operators.A
    C = operators.C
    C_astuple = operators.C_astuple
    B = operators.B
    B_astuple = operators.B_astuple

    vr, vi = A.createVecs()

    # mesh and degree as instance variables of ActiveFlame
    mesh = D.mesh
    degree = D.degree

    omega = np.zeros(maxiter, dtype=complex)
    omega[0] = init

    domega = 2 * tol
    k = 0

    # formatting
    tol_ = "{:.0e}".format(tol)
    tol_ = int(tol_[-2:])
    s = "{{:+.{}f}}".format(tol_)

    while abs(domega) > tol:

        D.assemble_matrix(omega[k])
        if not B:
            L = A + \
                mult_C(omega[k] ** 2, *C_astuple) - \
                D.matrix
            dL_domega = mult_C(2 * omega[k], *C_astuple) - D.get_derivative(omega[k])
        else:
            L = A + \
                mult_B(omega[k], *B_astuple) + \
                mult_C(omega[k] ** 2, *C_astuple) - \
                D.matrix
            dL_domega = B + mult_C(2 * omega[k], *C_astuple) - D.get_derivative(omega[k])

        # solve the eigenvalue problem L(\omega) * p = \lambda * C * p
        # set the target to zero (shift-and-invert)
        E = eps_solver(L, - C, 0, nev, two_sided=True, print_results=print_results)

        eig = E.getEigenvalue(i)

        E.getEigenvector(i, vr, vi)
        p = eigenvectors.normalize_2(mesh, vr, vi, degree)

        E.getLeftEigenvector(i + 1, vr, vi)
        p_adj = eigenvectors.normalize_2(mesh, vr, vi, degree)

        # convert into PETSc.Vec type
        p_vec = p.vector().vec()
        p_adj_vec = p_adj.vector().vec()

        # numerator and denominator
        num = vector_matrix_vector(p_adj_vec, dL_domega, p_vec)
        den = vector_matrix_vector(p_adj_vec, C, p_vec)

        deig = num / den
        domega = - eig / deig

        omega[k + 1] = omega[k] + domega

        print('iter = {:2d},  omega = {}  {}j,  |domega| = {:.2e}'.format(
            k, s.format(omega[k + 1].real), s.format(omega[k + 1].imag), abs(domega)))

        k += 1

    return E
