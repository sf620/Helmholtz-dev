import dolfin as dolf
from math import pi

from helmholtz_pkg.passive_flame import PassiveFlame
from helmholtz_pkg.flame_transfer_function import n_tau
from helmholtz_pkg.active_flame import ActiveFlame
from helmholtz_pkg.eigensolvers import fixed_point_iteration_pep
from helmholtz_pkg.eigenvectors import normalize_eigenvector

import params


def mshr(el):

    mesh = dolf.UnitIntervalMesh(el)

    def l_boundary_func(x, on_boundary):
        x = x[0]
        return on_boundary and dolf.near(x, 0.)

    def r_boundary_func(x, on_boundary):
        x = x[0]
        return on_boundary and dolf.near(x, 1.)

    boundaries = dolf.MeshFunction('size_t', mesh, mesh.topology().dim() - 1)

    l_boundary = dolf.AutoSubDomain(l_boundary_func)
    r_boundary = dolf.AutoSubDomain(r_boundary_func)

    l_boundary.mark(boundaries, 1)
    r_boundary.mark(boundaries, 2)

    # ________________________________________________________________________________

    def fl_subdomain_func(x):
        x = x[0]
        x_f = params.x_f[0][0]
        a_f = params.a_f
        return x_f - a_f - dolf.DOLFIN_EPS <= x <= x_f + a_f + dolf.DOLFIN_EPS

    subdomains = dolf.MeshFunction('size_t', mesh, mesh.topology().dim())

    subdomains.set_all(1)

    fl_subdomain = dolf.AutoSubDomain(fl_subdomain_func)
    fl_subdomain.mark(subdomains, 0)

    return mesh, boundaries, subdomains


def run():

    degree = 1

    mesh, boundaries, subdomains = mshr(400)

    boundary_conditions = {1: {'Robin': params.Y_in},  # inlet
                           2: {'Robin': params.Y_out}}  # outlet

    foo = PassiveFlame(mesh, boundaries, boundary_conditions,
                       c=params.c,
                       degree=degree)
    foo.assemble_A()
    foo.assemble_B()
    foo.assemble_C()

    ftf = n_tau(params.n, params.tau)

    D = ActiveFlame(mesh, subdomains,
                    params.x_f, params.x_r, params.rho_in, 1., 1., ftf,
                    degree=degree)
    D.assemble_submatrices()

    E = fixed_point_iteration_pep(foo, D, pi, nev=2, i=0)

    omega, p = normalize_eigenvector(mesh, E, i=0, degree=degree)

    # p_r, p_i = p.split()
    #
    # dolf.plot(p_r)
    # plt.show()
    # dolf.plot(p_i)
    # plt.show()


if __name__ == '__main__':

    import time

    t0 = time.time()

    run()

    print('eigenvalues computed in {:.3f} sec'.format(time.time() - t0))
