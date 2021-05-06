"""
All the functions in this module work in parallel
"""


def multiply(x0, z):
    """
    multiply a complex vector by a complex scalar
    """
    a = z.real
    b = z.imag

    x1 = x0.copy()
    x2 = x0.copy()

    x1[0::2] = - x0[1::2]
    x1[1::2] = + x0[0::2]

    x2.axpby(b, a, x1)

    return x2


def vector_vector(y0, x0):

    y1 = y0.copy()
    y1[0::2] = - y0[1::2]
    y1[1::2] = + y0[0::2]

    a = y0.dot(x0)
    b = y1.dot(x0)

    z = a + b * 1j

    return z


def vector_matrix_vector(y0, A, x0):

    x1 = x0.copy()
    A.mult(x0, x1)
    z = vector_vector(y0, x1)

    return z


def mult_complex_scalar_real_matrix(z, A_11, A_12, A_21, A_22):
    """
    multiply a real matrix by a complex scalar
    """
    a = z.real
    b = z.imag

    B = a*A_11 - b*A_12 + b*A_21 + a*A_22

    return B


def mult_complex_scalar_complex_matrix(z, A_r_11, A_r_12, A_r_21, A_r_22, A_i_11, A_i_12, A_i_21, A_i_22):
    """
    multiply a complex matrix by a complex scalar
    """
    a = z.real
    b = z.imag

    # (A_r_11, A_r_12, A_r_21, A_r_22, A_i_11, A_i_12, A_i_21, A_i_22) = A

    B =     a*A_r_11 - b*A_i_11   - a*A_i_12 - b*A_r_12   \
          + b*A_r_21 + a*A_i_21   - b*A_i_22 + a*A_r_22

    return B

# if __name__ == '__main__':
