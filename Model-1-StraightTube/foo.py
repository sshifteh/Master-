from dolfin import *
from WSS_bdry import bdry
import numpy as np

mesh = UnitSquareMesh(20, 20)
DG1 =FunctionSpace(mesh, 'DG', 1) 	# ShearStressSpace
DG0 = FunctionSpace(mesh, 'DG', 0)	# DPspace CHANGE this name to DG0 much more intuitive than discont. pressure space

w = 0.45
K_1  = Expression('(x[0] > w - DOLFIN_EPS && x[0] < 1 - w + DOLFIN_EPS) ? 0.0 : 1.0', w=w)   
K_1 = interpolate(K_1, DG0)

# Which edges are such that the connected cells have different values
facet_f = FacetFunction('size_t', mesh, 0)
K1_values = K_1.vector().array()
mesh.init(1, 2)
for facet in facets(mesh):
    connected_cells = facet.entities(2)
    dofs = [DG0.dofmap().cell_dofs(cell)[0] for cell in connected_cells]
    if len(dofs) > 1:
        one, two = K1_values[dofs]
        if abs(two-one) > 0:
            facet_f[facet] = 1

# Now that I found them, which cells are connected to them
cell_f = CellFunction('size_t', mesh, 0)
facet_f_values = facet_f.array()
mesh.init(2, 1)
for cell in cells(mesh):
    if np.sum(facet_f_values[cell.entities(1)]) > 0:
        cell_f[cell] = 1

plot(facet_f)
plot(K_1)
plot(cell_f)
interactive()
