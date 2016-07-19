def tube_domain(w, n = 20):
	mesh = UnitSquareMesh(n,n)       # Mesh over a unit square
	V = FunctionSpace(mesh, 'DG', 1) # Funcspace w DG as basisfunctions, i.e V = span{Ui*phi_i}, Phi_i DG basisfuncs for i = 0,..,n 
	f  = Expression('(x[0] > w && x[0] < 1 - w ) ? 3.0 : 10.0', w=w)   # This is a CompiledExpression Object
	c = 0.00005	
	f1 = Expression('x[0]< 1 - tanh(c*x[0])', c=c)
	k = interpolate(f1,V)                                               # This is a Function object k = sum(Ui*phi_i) and can be plotted 
	plot(k ,interactive = True, title = 'Tube sobdomain')

if __name__ == '__main__':
	
	from dolfin import * 
	w = 0.3
	tube_domain(w = w)


