import dolfin as dolf
import numpy as np
import scipy.linalg

from ufl import Dn

grad = dolf.grad
dot = dolf.dot
div = dolf.div

# ________________________________________________________________________________


def _shape_gradient_Dirichlet_real(c, p_dir, p_adj):
    return - c**2 * dot(grad(p_adj), grad(p_dir))


def _shape_gradient_Dirichlet_complex(c, p_dir, p_adj):
    p_dir = p_dir.split(True)
    p_adj = p_adj.split(True)
    real = - c**2 * (dot(grad(p_adj[0]), grad(p_dir[0])) +
                     dot(grad(p_adj[1]), grad(p_dir[1])))
    imag = - c**2 * (dot(grad(p_adj[0]), grad(p_dir[1])) -
                     dot(grad(p_adj[1]), grad(p_dir[0])))
    return real, imag


def _shape_gradient_Neumann_real(c, omega, p_dir, p_adj):
    return c**2 * dot(grad(p_adj), grad(p_dir)) - omega**2 * p_adj * p_dir


def _shape_gradient_Neumann_complex(c, omega, p_dir, p_adj):
    k = omega**2
    k_r, k_i = k.real, k.imag
    p_dir = p_dir.split(True)
    p_adj = p_adj.split(True)
    real = c**2 * (dot(grad(p_adj[0]), grad(p_dir[0])) +
                   dot(grad(p_adj[1]), grad(p_dir[1])))
    imag = c**2 * (dot(grad(p_adj[0]), grad(p_dir[1])) -
                   dot(grad(p_adj[1]), grad(p_dir[0])))
    f_r = p_adj[0] * p_dir[0] + p_adj[1] * p_dir[1]
    f_i = p_adj[0] * p_dir[1] - p_adj[1] * p_dir[0]
    real -= k_r * f_r - k_i * f_i
    imag -= k_r * f_i + k_i * f_r
    return real, imag


def _shape_gradient_Robin_complex(mesh, c, omega, p_dir, p_adj):

    V = dolf.FunctionSpace(mesh, 'CG', 1)
    c = dolf.interpolate(c, V)

    k = omega ** 2
    k_r, k_i = k.real, k.imag

    p_dir = p_dir.split(True)
    p_adj = p_adj.split(True)

    real = c ** 2 * (dot(grad(p_adj[0]), grad(p_dir[0])) +
                     dot(grad(p_adj[1]), grad(p_dir[1]))) - \
        2 * c ** 2 * (Dn(p_adj[0]) * Dn(p_dir[0]) +
                      Dn(p_adj[1]) * Dn(p_dir[1])) - \
        Dn(c ** 2) * (p_adj[0] * Dn(p_dir[0]) +
                      p_adj[1] * Dn(p_dir[1]))

    imag = c ** 2 * (dot(grad(p_adj[0]), grad(p_dir[1])) -
                     dot(grad(p_adj[1]), grad(p_dir[0]))) - \
        2 * c ** 2 * (Dn(p_adj[0]) * Dn(p_dir[1]) -
                      Dn(p_adj[1]) * Dn(p_dir[0])) - \
        Dn(c ** 2) * (p_adj[0] * Dn(p_dir[1]) -
                      p_adj[1] * Dn(p_dir[0]))

    f_r = p_adj[0] * p_dir[0] + p_adj[1] * p_dir[1]
    f_i = p_adj[0] * p_dir[1] - p_adj[1] * p_dir[0]

    real -= k_r * f_r - k_i * f_i
    imag -= k_r * f_i + k_i * f_r

    return real, imag

# ________________________________________________________________________________


def shape_derivatives(geometry, boundary_conditions, omega, p_dir, p_adj, c, local=False):
    """
    Supports real and complex, simple and 2-fold degenerate eigenpairs

    :param geometry: <class '..._pkg.mshr.Mesh'>
    :param boundary_conditions: <class 'dict'>
    :param omega: <class 'float'> or <class 'complex'>
    :param p_dir: <class 'tuple'> (or <class 'list'>)
    :param p_adj: <class 'tuple'> (or <class 'list'>)
    :param c: <class 'dolfin.function.constant.Constant'>, <....expression.Expression'> or <....function.Function'>
    :param local: <class 'bool'>, If True, divide by the surface area of the patch to obtain an average local shape
    derivative
    :return: <class 'dict'> containing the shape derivatives

    """

    mesh = geometry.mesh
    boundaries = geometry.boundaries
    # n = dolf.FacetNormal(mesh)
    V = dolf.FunctionSpace(mesh, 'CG', 1)
    ds = dolf.Measure('ds', subdomain_data=boundaries)

    G_Dir = []
    G_Neu = []
    G_Rob = []

    results = {}

    num_sub_spaces = p_dir[0].function_space().num_sub_spaces()  # real or complex

    shape_gradient_Dirichlet = _shape_gradient_Dirichlet_real
    shape_gradient_Neumann   = _shape_gradient_Neumann_real
    if num_sub_spaces == 2:
        shape_gradient_Dirichlet = _shape_gradient_Dirichlet_complex
        shape_gradient_Neumann   = _shape_gradient_Neumann_complex
    shape_gradient_Robin = _shape_gradient_Robin_complex

    # ________________________________________________________________________________
    # shape gradient, ufl forms for Dirichlet and Neumann boundary conditions

    if 'Dirichlet' in boundary_conditions.values():
        for p_adj_i in p_adj:
            for p_dir_j in p_dir:
                G_Dir.append(shape_gradient_Dirichlet(c, p_dir_j, p_adj_i))
    if 'Neumann' in boundary_conditions.values():
        for p_adj_i in p_adj:
            for p_dir_j in p_dir:
                G_Neu.append(shape_gradient_Neumann(c, omega, p_dir_j, p_adj_i))
    for value in boundary_conditions.values():
        if 'Robin' in value:
            for p_adj_i in p_adj:
                for p_dir_j in p_dir:
                    G_Rob.append(shape_gradient_Robin(mesh, c, omega, p_dir_j, p_adj_i))
            break

    # ________________________________________________________________________________
    # shape derivatives
    #
    # if the number of subspaces equals 0, the eigenvectors are real
    # (else) if the number of subspaces equals 2, mixed function space,
    # the eigenvectors are complex

    # for i, value in boundary_conditions.items():
    #
    #     C = dolf.Constant(1)
    #     if local:
    #         C = dolf.interpolate(C, V)
    #         A = dolf.assemble(C * ds(i))
    #         C = 1 / A
    #
    #     g = []
    #     if num_sub_spaces == 0:
    #         if value == 'Dirichlet':
    #             g = [dolf.assemble(C * G_ij * ds(i)) for G_ij in G_Dir]
    #         elif value == 'Neumann':
    #             g = [dolf.assemble(C * G_ij * ds(i)) for G_ij in G_Neu]
    #     elif num_sub_spaces == 2:
    #         if value == 'Dirichlet':
    #             g = [(dolf.assemble(C * G_ij[0] * ds(i)) +
    #                   1j * dolf.assemble(C * G_ij[1] * ds(i)))
    #                  for G_ij in G_Dir]
    #         elif value == 'Neumann':
    #             g = [(dolf.assemble(C * G_ij[0] * ds(i)) +
    #                   1j * dolf.assemble(C * G_ij[1] * ds(i)))
    #                  for G_ij in G_Neu]
    #     if len(g) == 1:
    #         # the eigenvalue is simple
    #         results[i] = g
    #     elif len(g) == 4:
    #         # the eigenvalues are 2-fold degenerate
    #         A = np.array(([g[0], g[1]],
    #                       [g[2], g[3]]))
    #         eig = scipy.linalg.eigvals(A)
    #         results[i] = eig.tolist()

    for i, value in boundary_conditions.items():

        if value == 'Dirichlet':
            G = G_Dir
        elif value == 'Neumann':
            G = G_Neu
        elif list(value.keys())[0] == 'Robin':
            G = G_Rob

        C = dolf.Constant(1)
        if local:
            C = dolf.interpolate(C, V)
            A = dolf.assemble(C * ds(i))
            C = 1 / A

        g = []
        if num_sub_spaces == 0:
            g = [dolf.assemble(C * G_ij * ds(i)) for G_ij in G]
        elif num_sub_spaces == 2:
            g = [(dolf.assemble(C * G_ij[0] * ds(i)) +
                  1j * dolf.assemble(C * G_ij[1] * ds(i)))
                 for G_ij in G]
        if len(g) == 1:
            # the eigenvalue is simple
            results[i] = g
        elif len(g) == 4:
            # the eigenvalues are 2-fold degenerate
            A = np.array(([g[0], g[1]],
                          [g[2], g[3]]))
            eig = scipy.linalg.eigvals(A)
            results[i] = eig.tolist()

    return results


def shape_derivatives_2(geometry, boundary_conditions, boundary_mark, omega, p_dir, p_adj, c):
    """
    Supports real and complex, simple and 2-fold degenerate eigenpairs

    :param geometry: <class '..._pkg.mshr.Mesh'>
    :param boundary_conditions: <class 'dict'>
    :param boundary_mark: <class 'int'>
    :param omega: <class 'float'> or <class 'complex'>
    :param p_dir: <class 'tuple'> (or <class 'list'>)
    :param p_adj: <class 'tuple'> (or <class 'list'>)
    :param c: <class 'dolfin.function.constant.Constant'>, <....expression.Expression'> or <....function.Function'>
    :return: <class 'dict'> containing the shape derivatives

    """

    i = boundary_mark

    mesh = geometry.mesh
    boundaries = geometry.boundaries
    V = dolf.FunctionSpace(mesh, 'CG', 1)
    ds = dolf.Measure('ds', subdomain_data=boundaries)

    G = []
    results = {}

    shape_gradient_Dirichlet = _shape_gradient_Dirichlet_complex
    shape_gradient_Neumann   = _shape_gradient_Neumann_complex
    shape_gradient_Robin = _shape_gradient_Robin_complex

    # ________________________________________________________________________________
    # shape gradient, ufl forms for Dirichlet and Neumann boundary conditions

    if boundary_conditions[i] == 'Dirichlet':
        for p_adj_i in p_adj:
            for p_dir_j in p_dir:
                G.append(shape_gradient_Dirichlet(c, p_dir_j, p_adj_i))
    elif boundary_conditions[i] == 'Neumann':
        for p_adj_i in p_adj:
            for p_dir_j in p_dir:
                G.append(shape_gradient_Neumann(c, omega, p_dir_j, p_adj_i))
    elif boundary_conditions[i] == 'Robin':
        for p_adj_i in p_adj:
            for p_dir_j in p_dir:
                G.append(shape_gradient_Robin(mesh, c, omega, p_dir_j, p_adj_i))

    # ________________________________________________________________________________
    # shape derivatives

    key = 0

    for _G_ij in G:

        for k in range(2):

            key += 1

            C = _G_ij[k]
            const = dolf.assemble(abs(C) * ds(i))
            print(dolf.assemble(C * ds(i)))

            g = [(dolf.assemble(C / const * G_ij[0] * ds(i)) +
                  1j * dolf.assemble(C / const * G_ij[1] * ds(i)))
                 for G_ij in G]

            # the eigenvalues are 2-fold degenerate
            A = np.array(([g[0], g[1]],
                          [g[2], g[3]]))
            eig = scipy.linalg.eigvals(A)

            results[key] = eig.tolist()  # change key

    return results
