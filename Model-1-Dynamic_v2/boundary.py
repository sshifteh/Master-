
from dolfin import *



def boundary(K_Func, DG0):

	import numpy as np

	mesh = K_Func.function_space().mesh()

	# Which edges are such that the connected cells have different values
	facet_f = FacetFunction('size_t', mesh, 0)
	K1_values = K_Func.vector().array()
	mesh.init(1, 2)   # in the mesh loop up cells that are connected to en edge mesh.init(2,3) means which side is connected to which pyramid
	# this bulids a table with index for
	# fenics knows for each cell which vertices make it up. 
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

        f = Function(DG0)
        f.vector()[:] = cell_f.array()

	return f






