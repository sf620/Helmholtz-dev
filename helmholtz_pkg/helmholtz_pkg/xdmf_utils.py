import dolfin as dolf


def save_xdmf(filename, function, comm=None):
    """
    Does not support mixed vector functions
    filename can be either 'f.xdmf' or 'f'
    """

    filename = filename.split('.')[0]
    filename = '{}.{}'.format(filename, 'xdmf')

    if comm:
        f_out = dolf.XDMFFile(comm, filename)
    else:
        f_out = dolf.XDMFFile(filename)

    f = function

    if f.ufl_element().family() == 'Mixed':
        f_r, f_i = f.split(deepcopy=True)
        f_out.write_checkpoint(f_r, 'f', 0, dolf.XDMFFile.Encoding.HDF5, False)
        f_out.write_checkpoint(f_i, 'f', 1, dolf.XDMFFile.Encoding.HDF5, True)  # append to file
    else:
        f_out.write_checkpoint(f, 'f', 0, dolf.XDMFFile.Encoding.HDF5, False)

    f_out.close()


def load_xdmf(filename, function_space, comm=None):
    """
    Does not support mixed vector functions
    filename can be either 'f.xdmf' or 'f'
    """

    filename = filename.split('.')[0]
    filename = '{}.{}'.format(filename, 'xdmf')

    if comm:
        f_in = dolf.XDMFFile(comm, filename)
    else:
        f_in = dolf.XDMFFile(filename)

    f = dolf.Function(function_space)

    if f.ufl_element().family() == 'Mixed':
        f_r, f_i = f.split(deepcopy=True)
        f_in.read_checkpoint(f_r, 'f', 0)
        f_in.read_checkpoint(f_i, 'f', 1)
        dolf.assign(f.sub(0), f_r)
        dolf.assign(f.sub(1), f_i)
    else:
        f_in.read_checkpoint(f, 'f', 0)

    f_in.close()

    return f
