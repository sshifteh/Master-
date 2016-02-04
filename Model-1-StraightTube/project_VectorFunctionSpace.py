
from dolfin import *

mesh = UnitSquareMesh(10,10)

V = FunctionSpace(mesh, 'CG', 1)
u = TrialFunction(V)
v = TestFunction(V)
f = Constant(-6.0)
a = inner(grad(u), grad(v))*dx
L = f*v*dx

def u0_boundary(x, on_bdry):
	return on_bdry 

u0 = Constant((1.,0.))
bs = DirichletBC(V, u0, u0_boundary)
u = Function(V)
solve(a==L, u, bc)



# Gradient 
V_g = VectorFunctionSpace(mesh, 'Lagrange', 1)
w = TrialFunction(V_g)
v = TestFunction(V_g)

a = inner(w,v)
L = inner(grad(u), v)*dx 
grad_u = function(V_g)
solve(a == L, grad_u)


