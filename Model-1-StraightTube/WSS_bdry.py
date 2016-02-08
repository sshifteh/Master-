
from dolfin import *

def WSS_bdry(K_Func, DPspace, WSS):

	"""

	The implicit function, K, creates separates the R**2. The zerocontour gives and interface.
	The outer region is assigned value 1, while the inner is assigned value 0.
	This can be expressed as a function of one variable by the Heaviside function
	By definition the directional derivative of the Heaviside function is the dirac delta function.
	The volume or surface integral of a funcion over a domain is defined as the product of the function and the drac delta function.
	The dirac delta function picks out the boundary.

	Something similar I am trying to do here.

	"""


	bdry= grad(K_Func)**2
	WSS_bdry = WSS*bdry # only if ur at the bdry AND the wss is nonzero
	# keept the bdry and wss separate
	# if at the bdry AND wss  larger grow the domain
	# if at bdry AND wss is less then shrink the domain


	WSS_bdry_Function =  project(WSS_bdry, DPspace)
	plot(WSS_bdry_Function, interactive = True, title = 'WSS_bdry_Function')

	return WSS_bdry_Function


"""
def bdry(K_Func, DG1):
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

"""


def bdry(mesh, K_Func, DG0):

	import numpy as np

	# Which edges are such that the connected cells have different values
	facet_f = FacetFunction('size_t', mesh, 0)
	K1_values = K_Func.vector().array()
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

        f = Function(DG0)
        f.vector()[:] = cell_f.array()

	return f






