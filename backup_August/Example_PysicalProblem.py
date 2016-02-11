# -*- coding: utf-8 -*-
"""
Created on Wed Aug 19 13:59:34 2015

@author: shiftehs
"""

# import statements -----------------------------------------
from dolfin import *
import numpy


# Step 1) Initializing the variables------------------------ 
T = float(10) #Tension
A = float(1)  #pressure amplitude 
R = 0.3       #radius of domain 
theta = 0.2   # angle 
x0 = 0.6*R*cos(theta)
y0 = 0.6*R*sin(theta)
s = 0.025               # sigma
s = 50                  # larger values for verification 



# Step 2) creating a mesh with bcs---------------------------- 
n = 40
mesh = UnitCircle(n)
V = FunctionSpace(mesh, "Lagrange", 1)

def boundary(x):
    return x[0] == 1E-8 or x[1] == 1E-8 or x[0] == 1-1E-8 or x[1] == 1-1E-8
     
bcs = DirichletBC(V, Constant(0.0), boundary)


# Step 3) Make an expression object for f ----------------------
f = Expression("4*exp(-0.5* (R*x[0]-x0/s)*(R*x[0]-x0/s) - 0.5*(R*x[1]-y0/s)*(R*x[1]-y0/s) )",R = R, x0 = x0, y0 = y0, s = s  )
f = interpolate(f,V)

# Step 4) Define variational problem --------------------------
w = TrialFunction(V)
v = TestFunction(V)
a = inner(grad(w), grad(v))*dx
L = f*v*dx

# Compute solution
w = Function(V)
problem = LinearVariationalProblem(a,L,w,bcs)
solver = LinearVariationalSolver(problem)
solver.parameters["linear_solver"] = "cg"   # conjugate gradient
solver.parameters["preconditioner"] = "ilu" # incomplete lu factorization
solver.solve()


# Plot scaled solution, mesh and pressure------------------------------
plot(mesh, title = "Mesh over the scaled domain")
plot(w, title = "Scaled deflection")
plot(f, title = "scaled pressure")

# Find maximum real deflection THIS PART I DIDNT DO!!!------------------ 
max_D = A*max_w/(8*pi*sigma*T)
print ’Maximum real deflection is’, max_D
# Verification for "flat" pressure (large sigma)
if sigma >= 50:
w_e = Expression("1 - x[0]*x[0] - x[1]*x[1]")
w_e = interpolate(w_e, V)
dev = numpy.abs(w_e.vector().array() - w.vector().array()).max()
print ’sigma=%g: max deviation=%e’ % (sigma, dev)
# Should be at the end
interactive()



























