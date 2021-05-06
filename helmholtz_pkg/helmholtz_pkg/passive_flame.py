import dolfin as dolf
from petsc4py import PETSc


def passive_flame(mesh, boundaries, boundary_conditions,
                  c, degree=1, constrained_domain=None):

    dx = dolf.Measure('dx', domain=mesh)

    V = dolf.FunctionSpace(mesh, 'CG', degree, constrained_domain=constrained_domain)

    bcs = []

    for i in boundary_conditions:
        if 'Dirichlet' in boundary_conditions[i]:
            bc = dolf.DirichletBC(V, 0.0, boundaries, i)
            bcs.append(bc)

    u = dolf.TrialFunction(V)
    v = dolf.TestFunction(V)

    a_ = - c ** 2 * dolf.dot(dolf.grad(u), dolf.grad(v)) * dx
    c_ = u * v * dx

    dummy = v * dx

    A, b = dolf.assemble_system(a_, dummy, bcs)
    A = dolf.as_backend_type(A).mat()

    C, b = dolf.assemble_system(c_, dummy, bcs)
    [bc.zero(C) for bc in bcs]
    C = dolf.as_backend_type(C).mat()

    N = V.dim()  # global size
    istart, iend = V.dofmap().ownership_range()
    n = iend - istart  # local size

    B = PETSc.Mat().create()
    B.setSizes([(n, N), (n, N)])
    B.setFromOptions()
    B.setUp()
    B.assemble()

    return A, B, C


class PassiveFlame:

    def __init__(self, mesh, boundaries, boundary_conditions,
                 c, degree=1, constrained_domain=None):

        self.mesh = mesh
        # self.boundaries = boundaries
        self.boundary_conditions = boundary_conditions
        self.c = c
        self.degree = degree

        self.dx = dolf.Measure('dx', domain=mesh)
        self.ds = dolf.Measure('ds', domain=mesh, subdomain_data=boundaries)

        CG = dolf.FiniteElement('CG', mesh.ufl_cell(), degree)
        W = dolf.FunctionSpace(mesh, CG * CG, constrained_domain=constrained_domain)
        self.function_space = W  #

        self.u = dolf.TrialFunction(W)
        self.v = dolf.TestFunction(W)

        self.bcs = []

        for i in boundary_conditions:
            if 'Dirichlet' in boundary_conditions[i]:
                for j in range(2):
                    bc = dolf.DirichletBC(W.sub(j), 0.0, boundaries, i)
                    self.bcs.append(bc)

        self._A = None
        self._B = None
        self._B_astuple = None
        self._B_adj = None
        self._B_adj_astuple = None
        self._C = None
        self._C_astuple = None

    @property
    def A(self):
        return self._A

    @property
    def B(self):
        return self._B

    @property
    def B_astuple(self):
        return self._B_astuple

    @property
    def B_adj(self):
        return self._B_adj

    @property
    def B_adj_astuple(self):
        return self._B_adj_astuple

    @property
    def C(self):
        return self._C

    @property
    def C_astuple(self):
        return self._C_astuple

    def assemble_A(self):

        (u_1, u_2) = self.u
        (v_1, v_2) = self.v

        a_11 = - self.c ** 2 * dolf.dot(dolf.grad(v_1), dolf.grad(u_1)) * self.dx
        a_22 = - self.c ** 2 * dolf.dot(dolf.grad(v_2), dolf.grad(u_2)) * self.dx
        a_ = a_11 + a_22

        dummy = (v_1 + v_2) * self.dx

        A, b = dolf.assemble_system(a_, dummy, self.bcs)
        A = dolf.as_backend_type(A).mat()
        self._A = A

    def assemble_B(self):

        (u_1, u_2) = self.u
        (v_1, v_2) = self.v

        N = self.function_space.dim()  # global size
        istart, iend = self.function_space.dofmap().ownership_range()
        n = iend - istart  # local size

        integrals_R = []

        for i in self.boundary_conditions:
            if 'Robin' in self.boundary_conditions[i]:

                Y = self.boundary_conditions[i]['Robin']
                Y_r, Y_i = Y.real, Y.imag

                b_11 = - Y_i * self.c * v_1 * u_1 * self.ds(i)
                b_12 = - Y_r * self.c * v_1 * u_2 * self.ds(i)
                b_21 = + Y_r * self.c * v_2 * u_1 * self.ds(i)
                b_22 = - Y_i * self.c * v_2 * u_2 * self.ds(i)

                b_ = b_11 + b_12 + b_21 + b_22

                integrals_R.append(b_)

        if integrals_R:

            b_ = sum(integrals_R)
            B = dolf.assemble(b_)
            B = dolf.as_backend_type(B).mat()

        else:

            B = PETSc.Mat().create()
            B.setSizes([(n, N), (n, N)])
            B.setFromOptions()
            B.setUp()
            B.assemble()

        B_adj = B.copy()
        B.transpose(B_adj)

        self._B = B
        self._B_adj = B_adj

    def assemble_B_astuple(self):

        (u_1, u_2) = self.u
        (v_1, v_2) = self.v

        (b_11_r, b_12_r, b_21_r, b_22_r, b_11_i, b_12_i, b_21_i, b_22_i) = 8 * (0, )

        for i in self.boundary_conditions:
            if 'Robin' in self.boundary_conditions[i]:

                Y = self.boundary_conditions[i]['Robin']
                Y_r, Y_i = Y.real, Y.imag

                b_11_r += - Y_i * self.c * v_1 * u_1 * self.ds(i)
                b_12_r += - Y_i * self.c * v_1 * u_2 * self.ds(i)
                b_21_r += - Y_i * self.c * v_2 * u_1 * self.ds(i)
                b_22_r += - Y_i * self.c * v_2 * u_2 * self.ds(i)

                b_11_i += + Y_r * self.c * v_1 * u_1 * self.ds(i)
                b_12_i += + Y_r * self.c * v_1 * u_2 * self.ds(i)
                b_21_i += + Y_r * self.c * v_2 * u_1 * self.ds(i)
                b_22_i += + Y_r * self.c * v_2 * u_2 * self.ds(i)

        b_ = (b_11_r, b_12_r, b_21_r, b_22_r, b_11_i, b_12_i, b_21_i, b_22_i)

        B = [dolf.as_backend_type(dolf.assemble(ufl_form)).mat() for ufl_form in b_]

        [B_11_r, B_12_r, B_21_r, B_22_r, B_11_i, B_12_i, B_21_i, B_22_i] = B
        B_adj = [B_11_r, B_12_r, B_21_r, B_22_r, - B_11_i, - B_12_i, - B_21_i, - B_22_i]

        self._B_astuple = (*B, )
        self._B_adj_astuple = (*B_adj, )

    def assemble_zB(self, z):

        z_r, z_i = z.real, z.imag

        (u_1, u_2) = self.u
        (v_1, v_2) = self.v

        integrals_R = []

        for i in self.boundary_conditions:
            if 'Robin' in self.boundary_conditions[i]:

                Y = self.boundary_conditions[i]['Robin']
                Y_r, Y_i = Y.real, Y.imag

                coeff_r = - z_r * Y_i - z_i * Y_r
                coeff_i = + z_r * Y_r - z_i * Y_i

                b_11 = + coeff_r * self.c * v_1 * u_1 * self.ds(i)
                b_12 = - coeff_i * self.c * v_1 * u_2 * self.ds(i)
                b_21 = + coeff_i * self.c * v_2 * u_1 * self.ds(i)
                b_22 = + coeff_r * self.c * v_2 * u_2 * self.ds(i)

                b_ = b_11 + b_12 + b_21 + b_22

                integrals_R.append(b_)

        if integrals_R:
            b_ = sum(integrals_R)
            Matrix = dolf.assemble(b_)
            Mat = dolf.as_backend_type(Matrix).mat()

            return Mat

    def assemble_C(self):

        (u_1, u_2) = self.u
        (v_1, v_2) = self.v

        c_11 = v_1 * u_1 * self.dx
        c_22 = v_2 * u_2 * self.dx
        c_   = c_11 + c_22

        dummy = (v_1 + v_2) * self.dx

        C, b = dolf.assemble_system(c_, dummy, self.bcs)
        [bc.zero(C) for bc in self.bcs]
        C = dolf.as_backend_type(C).mat()
        self._C = C

    def assemble_C_astuple(self):

        (u_1, u_2) = self.u
        (v_1, v_2) = self.v

        c_11 = v_1 * u_1 * self.dx
        c_12 = v_1 * u_2 * self.dx
        c_21 = v_2 * u_1 * self.dx
        c_22 = v_2 * u_2 * self.dx

        c_ = (c_11, c_12, c_21, c_22)

        dummy = (v_1 + v_2) * self.dx

        C = []

        for ufl_form in c_:
            Matrix, b = dolf.assemble_system(ufl_form, dummy, self.bcs)
            [bc.zero(Matrix) for bc in self.bcs]
            Mat = dolf.as_backend_type(Matrix).mat()
            C.append(Mat)

        self._C_astuple = (*C, )

    def assemble_zC(self, z):

        z_r, z_i = z.real, z.imag

        (u_1, u_2) = self.u
        (v_1, v_2) = self.v

        c_11 = + z_r * v_1 * u_1 * self.dx
        c_12 = - z_i * v_1 * u_2 * self.dx
        c_21 = + z_i * v_2 * u_1 * self.dx
        c_22 = + z_r * v_2 * u_2 * self.dx

        c_ = c_11 + c_12 + c_21 + c_22

        dummy = (v_1 + v_2) * self.dx

        Matrix, b = dolf.assemble_system(c_, dummy, self.bcs)
        [bc.zero(Matrix) for bc in self.bcs]
        Mat = dolf.as_backend_type(Matrix).mat()

        return Mat

# if __name__ == '__main__':
