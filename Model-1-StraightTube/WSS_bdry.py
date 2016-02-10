
from dolfin import *

def WSS_bdry(K_Func, DPspace, WSS):


	bdry= grad(K_Func)**2
	WSS_bdry = WSS*bdry # only if ur at the bdry AND the wss is nonzero
	# keept the bdry and wss separate
	
	WSS_bdry_Function =  project(WSS_bdry, DPspace)
	plot(WSS_bdry_Function, interactive = True, title = 'WSS_bdry_Function')

	return WSS_bdry_Function



def bdry_old(K_Func, DG1):
        n = FacetNormal(DG1.mesh())

        u = TrialFunction(DG1)
        v = TestFunction(DG1)

	a = avg(u)*avg(v)*dS() + v*u*ds()
	L = inner(jump(K_Func,n), jump(K_Func, n))*avg(v)*dS()
        #L = jump(K_Func)**2 * avg(v) * dS()

	bdry_func = Function(DG1)
        solve(a == L, bdry_func)
	bdry_cg1 = interpolate(bdry_func, FunctionSpace(refine(bdry_func.function_space().mesh()), "CG", 1))
	from IPython import embed; embed()
	plot(bdry_cg1, interactive = True, title = 'bdry')
	return bdry_func




def bdry(mesh, K_Func, DG0):

	import numpy as np

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






