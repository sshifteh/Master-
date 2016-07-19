# -*- coding: utf-8 -*-
# Copyright (C) 2008-2012 Joachim B. Haga and Fredrik Valdmanis
#
# This file is part of DOLFIN.
#
# DOLFIN is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# DOLFIN is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with DOLFIN. If not, see <http://www.gnu.org/licenses/>.
#
# Modified by Martin Sandve Alnæs 2008-2015
# Modified by Anders Logg 2008-2015
#
# First added:  2008-03-05
# Last changed: 2015-11-26

from __future__ import print_function
import os
import dolfin
import dolfin.cpp as cpp
import ufl
import numpy as np

__all__ = ['plot']

_meshfunction_types = (cpp.MeshFunction, cpp.MeshFunctionBool, cpp.MeshFunctionInt,
                      cpp.MeshFunctionDouble, cpp.MeshFunctionSizet)
_plottable_types = (cpp.Function,
                    cpp.Expression,
                    cpp.Mesh,
                    cpp.MultiMesh,
                    cpp.DirichletBC) + _meshfunction_types

# Add to parameter system (don't need this in C++)
cpp.parameters.add("plotting_backend", "vtk")

def compute_cell_values(f):
    mesh = f.function_space().mesh()
    V = dolfin.FunctionSpace(mesh, "DG", 0)
    v = dolfin.TestFunction(V)
    h = dolfin.CellVolume(mesh)
    return dolfin.assemble(f*v/h*dolfin.dx).array()

def mesh2triang(mesh):
    import matplotlib.tri as tri
    xy = mesh.coordinates()
    return tri.Triangulation(xy[:, 0], xy[:, 1], mesh.cells())

def mplot_mesh(ax, mesh, **kwargs):
    tdim = mesh.topology().dim()
    gdim = mesh.geometry().dim()
    if gdim == 2 and tdim == 2:
        return ax.triplot(mesh2triang(mesh), color='#808080')
    elif gdim == 3 and tdim == 3:
        bmesh = dolfin.BoundaryMesh(mesh, "exterior", order=False)
        mplot_mesh(ax, bmesh, **kwargs)
    elif gdim == 3 and tdim == 2:
        xy = mesh.coordinates()
        return ax.plot_trisurf(*[xy[:,i] for i in range(gdim)], triangles=mesh.cells())
    elif tdim == 1 and gdim == 1:
        x = mesh.coordinates()[:,0]
        return ax.plot(x, 0*x, 'o')
    elif tdim == 1:
        x = [mesh.coordinates()[:,i] for i in range(gdim)]
        return ax.plot(*x)
    else:
        raise AttributeError('Matplotlib plotting backend only supports 2D mesh.')

# TODO: This is duplicated somewhere else
def create_cg1_function_space(mesh, sh):
    r = len(sh)
    if r == 0:
        V = dolfin.FunctionSpace(mesh, "CG", 1)
    elif r == 1:
        V = dolfin.VectorFunctionSpace(mesh, "CG", 1, dim=sh[0])
    else:
        V = dolfin.TensorFunctionSpace(mesh, "CG", 1, shape=sh)
    return V

def mplot_expression(ax, f, mesh, **kwargs):
    # TODO: Can probably avoid creating the function space here by restructuring
    #       mplot_function a bit so it can handle Expression natively
    V = create_cg1_function_space(mesh, f.ufl_shape)
    g = dolfin.interpolate(f, V)
    return mplot_function(ax, g, **kwargs)

def mplot_function(ax, f, **kwargs):
    mesh = f.function_space().mesh()
    gdim = mesh.geometry().dim()
    tdim = mesh.topology().dim()

    if f.vector().size() == mesh.num_cells():
        # DG0 cellwise function
        C = f.vector().array() # NB! Assuming here dof ordering matching cell numbering
        if gdim == 2 and tdim == 2:
            return ax.tripcolor(mesh2triang(mesh), C)
        elif gdim == 3 and tdim == 2: # surface in 3d
            # FIXME: Not tested, probably broken
            xy = mesh.coordinates()
            return ax.plot_trisurf(mesh2triang(mesh), xy[:,2], C, shade=True)
        elif gdim == 1 and tdim == 1:
            x = mesh.coordinates()[:,0]
            nv = len(x)
            # Insert duplicate points to get piecewise constant plot
            xp = np.zeros(2*nv-2)
            xp[0] = x[0]
            xp[-1] = x[-1]
            xp[1:2*nv-3:2] = x[1:-1]
            xp[2:2*nv-2:2] = x[1:-1]
            Cp = np.zeros(len(xp))
            Cp[0:len(Cp)-1:2] = C
            Cp[1:len(Cp):2] = C
            return ax.plot(xp, Cp)
        #elif tdim == 1: # FIXME: Plot embedded line
        else:
            raise AttributeError('Matplotlib plotting backend only supports 2D mesh for scalar functions.')

    elif f.value_rank() == 0:
        # Scalar function, interpolated to vertices
        # TODO: Handle DG1?
        C = f.compute_vertex_values(mesh)
        if gdim == 2 and tdim == 2:
            mode = kwargs.get("mode", "surface")
            if mode == "surface":
                return ax.tripcolor(mesh2triang(mesh), C, shading='gouraud')
            elif mode == "warp":
                # FIXME: Make this work
                cellvalues = compute_cell_values(f)
                ax.plot_trisurf(mesh2triang(mesh), cellvalues, shade=True)
            elif mode == "wireframe":
                return ax.triplot(mesh2triang(mesh))
            elif mode == "contour":
                return ax.tricontour(mesh2triang(mesh), C)
        elif gdim == 3 and tdim == 2: # surface in 3d
            # FIXME: Not tested
            return ax.plot_trisurf(mesh2triang(mesh), C, shade=True)
        elif gdim == 3 and tdim == 3:
            # Volume
            # TODO: Isosurfaces?
            # Vertex point cloud
            X = [mesh.coordinates()[:, i] for i in range(gdim)]
            return ax.scatter(*X, c=C)
        elif gdim == 1 and tdim == 1:
            x = mesh.coordinates()[:,0]
            return ax.plot(x, C)
        #elif tdim == 1: # FIXME: Plot embedded line
        else:
            raise AttributeError('Matplotlib plotting backend only supports 2D mesh for scalar functions.')

    elif f.value_rank() == 1:
        # Vector function, interpolated to vertices
        w0 = f.compute_vertex_values(mesh)
        nv = mesh.num_vertices()
        if len(w0) != gdim*nv:
            raise AttributeError('Vector length must match geometric dimension.')
        X = mesh.coordinates()
        X = [X[:, i] for i in range(gdim)]
        U = [w0[i*nv: (i+1)*nv] for i in range(gdim)]

        # Compute magnitude
        C = U[0]**2
        for i in range(1,gdim):
            C += U[i]**2
        C = np.sqrt(C)

        mode = kwargs.get("mode", "glyphs")
        if mode == "glyphs":
            args = X + U + [C]
            if gdim == 3:
                return ax.quiver(*args, length=0.1) # TODO: Configurable length?
            else:
                return ax.quiver(*args)
        elif mode == "deformation":
            Xdef = [X[i] + U[i] for i in range(gdim)]
            import matplotlib.tri as tri
            triang = tri.Triangulation(Xdef[0], Xdef[1], mesh.cells())
            if gdim == 2 and tdim == 2:
                # FIXME: Not tested
                return ax.tripcolor(triang, C, shading='gouraud')
            else:
                raise AttributeError('Matplotlib plotting backend does not support deformation for %d in %d.' % (tdim, gdim))

def mplot_meshfunction(ax, obj, **kwargs):
    mesh = obj.mesh()
    tdim = mesh.topology().dim()
    d = obj.dim()
    if tdim == 2 and d == 2:
        C = obj.array()
        triang = mesh2triang(mesh)
        return ax.tripcolor(triang, facecolors=C)
    #elif tdim == 1 and d == 1:
    #elif d == 0: # vertices
    else:
        # FIXME:
        raise AttributeError("Matplotlib plotting backend doesn't handle meshfunction of dim %d" % d)

def mplot_dirichletbc(ax, obj, **kwargs):
    raise AttributeError("Matplotlib plotting backend doesn't handle DirichletBC.")

def _plot_matplotlib(obj, mesh, kwargs):
    # Avoid importing until used
    import matplotlib.pyplot as plt

    gdim = mesh.geometry().dim()
    if gdim == 3 or kwargs.get("mode") in ("warp",):
        # Importing this toolkit has side effects enabling 3d support
        from mpl_toolkits.mplot3d import axes3d
        # Enabling the 3d toolbox requires some additional arguments
        ax = plt.gca(projection='3d')
    else:
        ax = plt.gca()
    ax.set_aspect('equal')

    title = kwargs.pop("title", None)
    if title is not None:
        ax.set_title(title)

    if isinstance(obj, cpp.Function):
        return mplot_function(ax, obj, **kwargs)
    elif isinstance(obj, cpp.Expression):
        return mplot_expression(ax, obj, mesh, **kwargs)
    elif isinstance(obj, cpp.Mesh):
        return mplot_mesh(ax, obj, **kwargs)
    elif isinstance(obj, cpp.DirichletBC):
        return mplot_dirichletbc(ax, obj, **kwargs)
    elif isinstance(obj, _meshfunction_types):
        return mplot_meshfunction(ax, obj, **kwargs)
    else:
        raise AttributeError('Failed to plot %s' % type(obj))


# Compatibility with book
def _VTKPlotter_write_ps(self, *args, **kwargs) :
    print("*** Warning: VTKPlotter::write_ps() is not implemented -- use write_pdf instead")

_objects_referenced_from_plot_windows = {}
def _plot_cpp(object, mesh, kwargs):
    # Convert kwargs to cpp format
    p = cpp.Parameters()
    for key in kwargs:
        try:
            p.add(key, kwargs[key])
        except TypeError:
            cpp.warning("Incompatible type for keyword argument \"%s\". Ignoring." % key)

    if isinstance(object, cpp.Expression):
        plot_object = cpp.plot(object, mesh, p)
    elif isinstance(object, cpp.MultiMesh):
        return cpp.plot_multimesh(object)
    else:
        plot_object = cpp.plot(object, p)

    # Compatibility with book
    plot_object.write_ps = _VTKPlotter_write_ps

    # Avoid premature deletion of plotted objects if they go out of scope
    # before the plot window is closed. The plotter itself is safe, since it's
    # created in the plot() C++ function, not directly from Python. But the
    # Python plotter proxy may disappear, so we can't store the references
    # there.
    global _objects_referenced_from_plot_windows
    _objects_referenced_from_plot_windows[plot_object.key()] = (object, mesh, p)

    return plot_object

def plot(object, *args, **kwargs):
    """
    Plot given object.

    *Arguments*
        object
            a :py:class:`Mesh <dolfin.cpp.Mesh>`, a :py:class:`MeshFunction
            <dolfin.cpp.MeshFunction>`, a :py:class:`Function
            <dolfin.functions.function.Function>`, a :py:class:`Expression`
            <dolfin.cpp.Expression>, a :py:class:`DirichletBC`
            <dolfin.cpp.DirichletBC> or a :py:class:`FiniteElement
            <ufl.FiniteElement>`.

    *Examples of usage*
        In the simplest case, to plot only e.g. a mesh, simply use

        .. code-block:: python

            mesh = UnitSquare(4,4)
            plot(mesh)

        Use the ``title`` argument to specify title of the plot

        .. code-block:: python

            plot(mesh, tite="Finite element mesh")

        It is also possible to plot an element

        .. code-block:: python

            element = FiniteElement("BDM", tetrahedron, 3)
            plot(element)

        Vector valued functions can be visualized with an alternative mode

        .. code-block:: python

            plot(u, mode = "glyphs")

        A more advanced example

        .. code-block:: python

            plot(u,
                 wireframe = True,              # use wireframe rendering
                 interactive = False,           # do not hold plot on screen
                 scalarbar = False,             # hide the color mapping bar
                 hardcopy_prefix = "myplot",    # default plotfile name
                 scale = 2.0                    # scale the warping/glyphs
                 title = "Fancy plot"           # Set your own title
                 )

    """

    # Plot element
    if isinstance(object, ufl.FiniteElementBase):
        if os.environ.get("DOLFIN_NOPLOT", "0") != "0":
            return
        import ffc
        return ffc.plot(object, *args, **kwargs)

    # Get mesh from explicit mesh kwarg, only positional arg, or via object
    mesh = kwargs.pop('mesh', None)
    if isinstance(object, cpp.Mesh):
        if mesh is not None and mesh.id() != object.id():
            cpp.dolfin_error("plotting.py",
                             "plot mesh",
                             "Got different mesh in plot object and keyword argument")
        mesh = object
    if mesh is None:
        if isinstance(object, cpp.Function):
            mesh = object.function_space().mesh()
        elif hasattr(object, "mesh"):
            mesh = object.mesh()

    # Expressions do not carry their own mesh
    if isinstance(object, cpp.Expression) and mesh is None:
        cpp.dolfin_error("plotting.py",
                         "plot expression",
                         "Expecting a mesh as keyword argument")

    # Try to project if object is not a standard plottable type
    if not isinstance(object, _plottable_types):
        from dolfin.fem.projection import project
        try:
            cpp.info("Object cannot be plotted directly, projecting to"\
                     " piecewise linears.")
            object = project(object, mesh=mesh)
        except Exception as e:
            msg = "Don't know how to plot given object:\n  %s\n" \
                  "and projection failed:\n  %s" % (str(object), str(e))
            cpp.dolfin_error("plotting.py", "plot object", msg)

    # Select backend
    backend = cpp.parameters["plotting_backend"]
    if backend == "vtk":
        return _plot_cpp(object, mesh, kwargs)
    elif backend == "matplotlib":
        return _plot_matplotlib(object, mesh, kwargs)
