from datetime import datetime
import logging
import numpy as np

from .petsc4py_utils import mult_complex_scalar_real_matrix as mult_c
from .petsc4py_utils import mult_complex_scalar_complex_matrix as mult_b
from .petsc4py_utils import vector_matrix_vector

from .eigensolvers import eps_solver
from .eigenvectors import normalize_2

# ________________________________________________________________________________

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

formatter = logging.Formatter('%(message)s')

filename = '{}.log'.format(datetime.now().strftime('%Y%b%d-%H%M%S'))
file_handler = logging.FileHandler(filename)
file_handler.setFormatter(formatter)

logger.addHandler(file_handler)

# ________________________________________________________________________________


def fixed_point_iteration(operators, D,
                          target, nev=2, i=0,
                          tol=1e-8, maxiter=50,
                          print_results=False,
                          two_sided=True):

    A = operators.A
    C = operators.C
    C_astuple = operators.C_astuple
    B = operators.B
    B_astuple = operators.B_astuple
    B_adj_astuple = None
    if not two_sided:
        B_adj_astuple = operators.B_adj_astuple

    E = eps_solver(A, C, target, nev, print_results=print_results)
    eig = E.getEigenvalue(i)

    domega = 2 * tol
    k = - 1

    # formatting
    s = "{:.0e}".format(tol)
    s = int(s[-2:])
    s = "{{:+.{}f}}".format(s)

    vr, vi = A.createVecs()

    mesh = D.mesh
    degree = D.degree

    f = np.zeros(maxiter, dtype=complex)
    f_adj = np.zeros(maxiter, dtype=complex)
    df_domega = np.zeros(maxiter, dtype=complex)
    alpha = np.zeros(maxiter, dtype=complex)
    omega = np.zeros(maxiter, dtype=complex)
    omega_adj = np.zeros(maxiter, dtype=complex)

    omega[0] = np.sqrt(eig)
    omega_adj[0] = np.sqrt(eig)

    print('iter = {:2d},  omega = {}  {}j'.format(
        k + 1, s.format(omega[k + 1].real), s.format(omega[k + 1].imag)
    ))

    logger.info("iteration {:2d}".format(
        k + 1
    ))

    logger.info("-------------------------------------------------------\n")

    logger.info("omega = {}  {}j\n\n".format(
        s.format(omega[k + 1].real), s.format(omega[k + 1].imag)
    ))

    E_adj = None

    while abs(domega) > tol:

        k += 1

        D.assemble_matrix(omega[k])
        if not B:
            nlinA = A - D.matrix
            dnlinA_domega = - D.get_derivative(omega[k])
        else:
            nlinA = A + mult_b(omega[k], *B_astuple) - D.matrix
            dnlinA_domega = B - D.get_derivative(omega[k])

        nlinA_adj = None
        if not two_sided:
            D.assemble_matrix(omega_adj[k], 'adjoint')
            if not B:
                nlinA_adj = A - D.adjoint_matrix
            else:
                nlinA_adj = A + mult_b(omega_adj[k], *B_adj_astuple) - D.adjoint_matrix

        E = eps_solver(nlinA, C, target, nev, two_sided=two_sided, print_results=print_results)
        eig = E.getEigenvalue(i)
        f[k] = np.sqrt(eig)
        E.getEigenvector(i, vr, vi)
        p = normalize_2(mesh, vr, vi, degree)

        if not two_sided:
            E_adj = eps_solver(nlinA_adj, C, target, nev, print_results=print_results)
            eig = E_adj.getEigenvalue(i + 1)
            f_adj[k] = np.sqrt(eig)
            E_adj.getEigenvector(i + 1, vr, vi)
        else:
            E.getLeftEigenvector(i + 1, vr, vi)
        p_adj = normalize_2(mesh, vr, vi, degree)

        p_1 = p.vector().vec()
        p_adj_1 = p_adj.vector().vec()
        num = vector_matrix_vector(p_adj_1, dnlinA_domega, p_1)
        den = vector_matrix_vector(p_adj_1, mult_c(2 * f[k], *C_astuple), p_1)
        df_domega[k] = - num / den

        alpha[k] = 1 / (1 - df_domega[k])
        omega[k+1] = alpha[k] * f[k] + (1 - alpha[k]) * omega[k]
        domega = omega[k + 1] - omega[k]

        domega_adj = None
        if not two_sided:
            omega_adj[k + 1] = alpha[k].conjugate() * f_adj[k] + (1 - alpha[k].conjugate()) * omega_adj[k]
            domega_adj = omega_adj[k + 1] - omega_adj[k]

        print('iter = {:2d},  omega = {}  {}j,  |domega| = {:.2e}'.format(
            k + 1, s.format(omega[k + 1].real), s.format(omega[k + 1].imag), abs(domega)
        ))

        logger.info("iteration {:2d}".format(
            k + 1
        ))

        logger.info("-------------------------------------------------------\n")

        logger.info("f  = {}  {}j".format(
            s.format(f[k].real), s.format(f[k].imag)
        ))  # double space after f...

        logger.info("f' = {}  {}j,  |f'| = {}".format(
            s.format(df_domega[k].real), s.format(df_domega[k].imag), s.format(abs(df_domega[k]))
        ))

        logger.info("alpha = {}  {}j\n".format(
            s.format(alpha[k].real), s.format(alpha[k].imag)
        ))

        logger.info("omega = {}  {}j,  |domega| = {:.2e}\n".format(
            s.format(omega[k + 1].real), s.format(omega[k + 1].imag), abs(domega)
        ))

        if not two_sided:

            logger.info("adjoint (not two_sided)...\n")

            logger.info("f  = {}  {}j".format(
                s.format(f_adj[k].real), s.format(f_adj[k].imag)
            ))  # double space after f...

            logger.info("omega = {}  {}j,  |domega| = {:.2e}\n\n".format(
                s.format(omega_adj[k + 1].real), s.format(omega_adj[k + 1].imag), abs(domega_adj)
            ))

        else:

            logger.info('')

    if not two_sided:
        return E, E_adj
    else:
        return E
