import dolfin as dolf
from dolfin_utils.meshconvert import meshconvert
import meshio
import os


class MeshXDMF:

    def __init__(self, file, write_xdmf_file=True, rank=0):

        self.file = file
        self.write_xdmf_file = write_xdmf_file
        self.rank = rank

        self._mesh = None
        self._subdomains = None
        self._boundaries = None

    def __call__(self):

        msh_file = '{}.{}'.format(self.file, 'msh')
        xdmf_file = '{}.{}'.format(self.file, 'xdmf')

        if self.rank == 0:
            if self.write_xdmf_file:
                convert2xdmf(msh_file)
            else:
                if not os.path.isfile(xdmf_file):
                    convert2xdmf(msh_file)

        self._mesh = dolf.Mesh()
        with dolf.XDMFFile('{}.{}'.format(self.file, 'xdmf')) as infile:
            infile.read(self._mesh)

        dim = self._mesh.geometric_dimension()

        facet_region = '{}_{}.{}'.format(self.file, 'facet_region', 'xdmf')
        if os.path.isfile(facet_region):
            # mesh value collection boundaries
            mvcb = dolf.MeshValueCollection("size_t", self._mesh, dim - 1)
            with dolf.XDMFFile(facet_region) as infile:
                infile.read(mvcb, "name_to_read")
                self._boundaries = dolf.cpp.mesh.MeshFunctionSizet(self._mesh, mvcb)
        else:
            print('No {}'.format(facet_region))
            print('boundaries = None')

        physical_region = '{}_{}.{}'.format(self.file, 'physical_region', 'xdmf')
        if os.path.isfile(physical_region):
            # mesh value collection subdomains
            mvcs = dolf.MeshValueCollection("size_t", self._mesh, dim)
            with dolf.XDMFFile(physical_region) as infile:
                infile.read(mvcs, "name_to_read")
            self._subdomains = dolf.cpp.mesh.MeshFunctionSizet(self._mesh, mvcs)
        else:
            print('No {}'.format(physical_region))
            print('subdomains = None')

    @property
    def mesh(self):
        return self._mesh

    @property
    def subdomains(self):
        return self._subdomains

    @property
    def boundaries(self):
        return self._boundaries


def convert2xdmf(msh_file):
    """
    call _2d or _3d depending on the cell type (triangle or tetra)
    _2d and _3d convert the msh file into a xdmf file
    if there are physical groups, xdmf files are created for the boundaries and the subdomains
    """

    file = msh_file.split('.')[0]
    msh = meshio.read(msh_file)

    dim = 2  # 2d mesh (default)
    for cell in msh.cells:
        if cell.type == "tetra":
            dim = 3  # 3d mesh

    if dim == 2:
        _2d(msh, file)
    elif dim == 3:
        _3d(msh, file)


def _2d(msh, file):

    for cell in msh.cells:
        if cell.type == "line":
            line_cells = cell.data
        elif cell.type == "triangle":
            triangle_cells = cell.data

    for key, value in msh.cell_data_dict["gmsh:physical"].items():
        if key == "line":
            line_data = value
        elif key == "triangle":
            triangle_data = value

    mesh = meshio.Mesh(points=msh.points, cells={"triangle": triangle_cells})

    facet_region = meshio.Mesh(points=msh.points,
                               cells=[("line", line_cells)],
                               cell_data={"name_to_read": [line_data]})

    physical_region = meshio.Mesh(points=msh.points,
                                  cells=[("triangle", triangle_cells)],
                                  cell_data={"name_to_read": [triangle_data]})

    # Remove third (z) component of points if it is 0 everywhere
    mesh.prune_z_0()
    facet_region.prune_z_0()
    physical_region.prune_z_0()

    meshio.write('{}.{}'.format(file, 'xdmf'), mesh)
    meshio.write('{}_{}.{}'.format(file, 'facet_region', 'xdmf'), facet_region)
    meshio.write('{}_{}.{}'.format(file, 'physical_region', 'xdmf'), physical_region)


def _3d(msh, file):

    for cell in msh.cells:
        if cell.type == "triangle":
            triangle_cells = cell.data
        elif cell.type == "tetra":
            tetra_cells = cell.data

    for key, value in msh.cell_data_dict["gmsh:physical"].items():
        if key == "triangle":
            triangle_data = value
        elif key == "tetra":
            tetra_data = value

    mesh = meshio.Mesh(points=msh.points, cells={"tetra": tetra_cells})

    facet_region = meshio.Mesh(points=msh.points,
                               cells=[("triangle", triangle_cells)],
                               cell_data={"name_to_read": [triangle_data]})

    physical_region = meshio.Mesh(points=msh.points,
                                  cells=[("tetra", tetra_cells)],
                                  cell_data={"name_to_read": [tetra_data]})

    meshio.write('{}.{}'.format(file, 'xdmf'), mesh)
    meshio.write('{}_{}.{}'.format(file, 'facet_region', 'xdmf'), facet_region)
    meshio.write('{}_{}.{}'.format(file, 'physical_region', 'xdmf'), physical_region)

# ______________________________________________________________________________________________________________________


class MeshXML:

    def __init__(self, file, write_xml_file=True):

        self.file = file
        self.write_xml_file = write_xml_file

        self._mesh = None
        self._subdomains = None
        self._boundaries = None

    def __call__(self):

        msh_file = '{}.{}'.format(self.file, 'msh')
        xml_file = '{}.{}'.format(self.file, 'xml')

        if self.write_xml_file:
            meshconvert.convert2xml(msh_file, xml_file)
        else:
            if not os.path.isfile(xml_file):
                meshconvert.convert2xml(msh_file, xml_file)

        self._mesh = dolf.Mesh(xml_file)

        facet_region = '{}_{}.{}'.format(self.file, 'facet_region', 'xml')
        if os.path.isfile(facet_region):
            self._boundaries = dolf.MeshFunction('size_t', self._mesh, facet_region)
        else:
            print('No {}'.format(facet_region))
            print('boundaries = None')

        physical_region = '{}_{}.{}'.format(self.file, 'physical_region', 'xml')
        if os.path.isfile(physical_region):
            self._subdomains = dolf.MeshFunction('size_t', self._mesh, physical_region)
        else:
            print('No {}'.format(physical_region))
            print('subdomains = None')

    @property
    def mesh(self):
        return self._mesh

    @property
    def subdomains(self):
        return self._subdomains

    @property
    def boundaries(self):
        return self._boundaries


def mesh_xml(file, write_xml_file=True):

    mesh = None
    boundaries = None
    subdomains = None

    def call():

        msh_file = '{}.{}'.format(file, 'msh')
        xml_file = '{}.{}'.format(file, 'xml')

        if write_xml_file:
            meshconvert.convert2xml(msh_file, xml_file)
        else:
            if not os.path.isfile(xml_file):
                meshconvert.convert2xml(msh_file, xml_file)

        nonlocal mesh
        mesh = dolf.Mesh(xml_file)

        facet_region = '{}_{}.{}'.format(file, 'facet_region', 'xml')
        if os.path.isfile(facet_region):
            nonlocal boundaries
            boundaries = dolf.MeshFunction('size_t', mesh, facet_region)
        else:
            print('No {}'.format(facet_region))
            print('boundaries = None')

        physical_region = '{}_{}.{}'.format(file, 'physical_region', 'xml')
        if os.path.isfile(physical_region):
            nonlocal subdomains
            subdomains = dolf.MeshFunction('size_t', mesh, physical_region)
        else:
            print('No {}'.format(physical_region))
            print('subdomains = None')

    def get_mesh():
        return mesh

    def get_boundaries():
        return boundaries

    def get_subdomains():
        return subdomains

    call.get_mesh = get_mesh
    call.get_boundaries = get_boundaries
    call.get_subdomains = get_subdomains

    return call
