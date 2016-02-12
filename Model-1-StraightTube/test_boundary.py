def test_boundary():
	n = 10
	w= 0.45
	mesh = UnitSquareMesh(n,n) #'crossed')
	
	DG1 =FunctionSpace(mesh, 'DG', 1)
	DG0 = FunctionSpace(mesh, 'DG', 0)
	
	Pspace = FunctionSpace(mesh, 'CG',1)
	Vspace = VectorFunctionSpace(mesh, 'Lagrange', 2)
	Pspace = FunctionSpace(mesh, 'Lagrange', 1)
	Wspace = MixedFunctionSpace([Vspace, Pspace])
	RT = FunctionSpace(mesh, 'CR',1)
	
	Kspace = DG0
	WSSspace = DG1
	
	K = Expression('(x[0] > w - DOLFIN_EPS && x[0] < 1 - w + DOLFIN_EPS) ? 0.0 : 1.0', w=w)
	K_func = interpolate(K,Pspace)	
	plot(K_func)
	
	n = FacetNormal(mesh)
	

	
		
	
	interactive()
if __name__ == '__main__':
	from dolfin import * 	
	test_boundary()
