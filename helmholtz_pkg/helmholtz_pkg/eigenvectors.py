"""
    normalize_1 and normalize_2 have been tested in parallel
    normalize_eigenvector and normalize_adjoint have not
"""

import dolfin as dolf
import numpy as np
from slepc4py import SLEPc

from .petsc4py_utils import multiply, vector_matrix_vector


def normalize_1(mesh, vr, degree=1, constrained_domain=None):

    V = dolf.FunctionSpace(mesh, 'CG', degree, constrained_domain=constrained_domain)

    # print(vr.getArray())

    # index, value = abs(vr).max()
    # # print(index, value)
    #
    # vr /= value
    # # print(vr.getArray())

    p = dolf.Function(V)
    p.vector().set_local(vr.getArray())
    p.vector().apply('insert')

    dx = dolf.Measure('dx')
    meas = dolf.assemble(p * p * dx)
    meas = np.sqrt(meas)
    # print(meas)
    vr /= meas
    p.vector().set_local(vr.getArray())
    p.vector().apply('insert')
    # meas = dolf.assemble(p * p * dx)
    # meas = np.sqrt(meas)
    # print(meas)

    return p


def normalize_2(mesh, vr, vi, degree=1, constrained_domain=None):

    # local_y = np.max(np.abs(y))
    # global_y = MPI.COMM_WORLD.allreduce(sendobj=local_y, op=MPI.MAX)
    # y /= global_y

    V = dolf.FunctionSpace(mesh, "CG", degree, constrained_domain=constrained_domain)
    CG = dolf.FiniteElement("CG", mesh.ufl_cell(), degree)
    W = dolf.FunctionSpace(mesh, dolf.MixedElement([CG, CG]), constrained_domain=constrained_domain)

    vr_1 = vr.getArray()[0::2]
    vr_2 = vr.getArray()[1::2]
    vi_1 = vi.getArray()[0::2]
    vi_2 = vi.getArray()[1::2]

    x = vr_1 - vi_2 + 1j * (vr_2 + vi_1)

    x_r = x.real
    x_i = x.imag

    p_r = dolf.Function(V)
    p_i = dolf.Function(V)

    p_r.vector().set_local(x_r)
    p_r.vector().apply('insert')
    p_i.vector().set_local(x_i)
    p_i.vector().apply('insert')

    dx = dolf.Measure('dx')
    meas = dolf.assemble((p_r * p_r + p_i * p_i) * dx)
    meas = np.sqrt(meas)
    # print(meas)

    x /= meas

    x_r = x.real
    x_i = x.imag

    # p_r.vector().set_local(x_r)
    # p_r.vector().apply('insert')
    # p_i.vector().set_local(x_i)
    # p_i.vector().apply('insert')
    #
    # meas = dolf.assemble((p_r * p_r + p_i * p_i) * dx)
    # meas = np.sqrt(meas)
    # # print(meas)

    x = vr.copy()

    istart, iend = x.getOwnershipRange()

    x[istart:iend:2] = x_r
    x[istart+1:iend+1:2] = x_i

    p = dolf.Function(W)
    p.vector().set_local(x.getArray())
    p.vector().apply('insert')

    return p


def normalize_eigenvector(mesh, obj, i, degree=1, which='right'):

    omega = 0.
    A = obj.getOperators()[0]
    vr, vi = A.createVecs()

    if isinstance(obj, SLEPc.EPS):
        eig = obj.getEigenvalue(i)
        omega = np.sqrt(eig)
        if which == 'right':
            obj.getEigenvector(i, vr, vi)
        elif which == 'left':
            obj.getLeftEigenvector(i, vr, vi)

    elif isinstance(obj, SLEPc.PEP):
        eig = obj.getEigenpair(i, vr, vi)
        omega = eig

    if not np.any(vi):  # all zeros
        omega = omega.real
        p = normalize_1(mesh, vr, degree)
    else:
        p = normalize_2(mesh, vr, vi, degree)

    return omega, p


def normalize_adjoint(omega, p_dir, p_adj, matrices, D=None):
    """
    p_dir and p_adj  are both: <class 'dolfin.function.function.Function'>
    p_dir_vec and p_adj_vec are both: <class 'petsc4py.PETSc.Vec'>

    """

    B = matrices.B

    p_dir_vec = p_dir.vector().vec()
    p_adj_vec = p_adj.vector().vec()

    if not B and not D:
        print('not B and not D: return None')
        return None
    elif B and not D:
        # B + 2 \omega C
        dL_domega = (B +
                     matrices.assemble_zC(2 * omega))
    elif D and not B:
        # 2 \omega C - D'(\omega)
        dL_domega = (matrices.assemble_zC(2 * omega) -
                     D.get_derivative(omega))
    else:
        # B + 2 \omega C - D'(\omega)
        dL_domega = (B +
                     matrices.assemble_zC(2 * omega) -
                     D.get_derivative(omega))

    meas = vector_matrix_vector(p_adj_vec, dL_domega, p_dir_vec)
    # print(meas)

    p_adj_vec = multiply(p_adj_vec, 1 / meas.conjugate())

    # meas = vector_matrix_vector(p_adj_vec, dL_domega, p_dir_vec)
    # print(meas)

    p_adj1 = p_adj.copy(deepcopy=True)
    p_adj1.vector().set_local(p_adj_vec.getArray())
    p_adj1.vector().apply('insert')

    return p_adj1
