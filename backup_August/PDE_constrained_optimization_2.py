# -*- coding: utf-8 -*-

""" 
Solver optimial control model for varme lignignen 

""" 

# Importer pakker.................................................
from dolfin import * 
from dolfin_adjoint import * 

# Lag mesh over geometrien........................................

n = 200
mesh = RectangleMesh(-1, -1, 1, 1, n, n)
v = FunctionSpace(mesh, 'CG', 1)
u = Function(V, name = 'state')
m = Function(V, name = 'control')
v = TestFunction(V)

# Run the forward model once to create the simulation record......
F = (inner(grad(u), grad(v)) - m*v)*dx
bc = DirichletBC(V, 0.0, "on_boundary")
solve(F == 0, u, bc)

# The functional of interest is the normed difference between desired
# and simulated temperature profile (misfit) 
x = triangle.x
u_desired = exp(-1/(1-x[0]*x[0])-1/(1-x[1]*x[1]))  # exact solution
J = Functional((0.5*inner(u-u_desired, u-u_desired))*dx*dt[FINISH_TIME])

# Run the optimisation
reduced_functional = ReducedFunctional(J, Control(m, value=m))
# Make sure you have scipy >= 0.11 installed
m_opt = minimize(reduced_functional, method = "L-BFGS-B",
                 tol=2e-08, bounds = (-1, 1), options = {"disp": True})