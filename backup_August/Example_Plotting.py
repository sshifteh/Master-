# -*- coding: utf-8 -*-
"""
Created on Wed Aug 19 14:43:22 2015

@author: shiftehs

For plotting in Fenics we can use something called Viper.
Viper uses a package called VTK.
Viper is a thin layer on top of the Python interface from VTK.

I am making an example program to test the plotting functionalities.
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


# Plotting ------------------------------------------------------------------
# Ploting automatically  
plot(mesh, title = "Finite Element Mesh ")
plot(u_, wirefram = True, title = "Solution")

# Plotting by grabbing the Viper object and making adjustments 

viz_u = plot(u_,wirefram = False, title = "Poisson numerical solution", rescale = False, axes = True, basename = "Poisson")
viz_u.update(u_)
# bring settings above into action
viz_u.write_png("Poisson_plotting.png")
interactive(True)