from dolfin import *

n = 5 
mesh = UnitIntervalMesh(n)
V = FunctionSpace(mesh, 'Lagrange', 1)
u = TrialFunction(V)
v = TestFunction(V)
f = Constant(1.0) # this means that the heat is distributed evenly
g = Constant(1.0) # what does this mean 
u_ = Function(V)  # a function in a finite element function space is a function that is a linear combination of the basisfunctions which make up or in math terms make up that space. 

a = inner(grad(u), grad(v))*dx
L = f*v*dx + g*v*ds
solve(a == L, u_)
plot(u_, interactive = True, title = 'my solution')

