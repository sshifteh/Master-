# -*- coding: utf-8 -*-
"""
Created on Wed Aug 19 15:23:46 2015

@author: shiftehs
"""


from dolfin import * 

# Making a mesh-------------------------------------------------------------- 
n = 10 
mesh = UnitSquareMesh(n,n)
def boundary(x):
    return x[0] == 1E-8 or x[1] == 1E-8 or x[0] == 1-1E-8 or x[1] == 1-1E-8 

# Making an expression object for f---------------------------------------- 
f = Constant(0.0)  

# Defining the variational problem--------------------------------------------

V = FunctionSpace(mesh, "Lagrange", 1)
bcs = DirichletBC(V,Constant(0.0), boundary)
u = TrialFunction(V)
v = TestFunction(V)
a = inner(grad(u),grad(v))*dx
L = f*v*dx


# Solving the variational problem by a linear algebra package-----------------
# Compute solution
u_ = Function(V)
problem = LinearVariationalProblem(a, L, u_, bcs)
solver = LinearVariationalSolver(problem)
solver.parameters["linear_solver"] = "cg"    #conjugate gradient
solver.parameters["preconditioner"] = "ilu"  # incomplete LU factorization
solver.solve()


 
mesh = UnitSquareMesh(n,n)
V_g = VectorFunctionSpace(mesh, "Lagrange", 1)
w = TrialFunction(V_g)
v = TestFunction(V_g)
a = inner(w, v)*dx
L = inner(grad(u_), v)*dx
grad_u = Function(V_g)
solve(a == L, grad_u)
plot(grad_u, title="grad(u)")