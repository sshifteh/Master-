from dolfin import * 

n = 2
mesh = UnitSquareMesh(n,n)
#plot(mesh, interactive = True)
V = FunctionSpace(mesh, 'CG',1)

u = TrialFunction(V)
v = TestFunction(V)
f = Constant(-6) 
a = inner(grad(u), grad(v))*dx  
L = inner(f, v)*dx
u_solved = Function(V)

def Domain_boundary(x,on_bnd):
	return on_bnd

bcs  = DirichletBC(V, Constant((0.0)), Domain_boundary)

solve(a == L, u_solved, bcs)

plot(u_solved, interactive = True)
for u in u_solved.vector():
	print u 
print 'len of array', len(u_solved.vector()) # 1681





