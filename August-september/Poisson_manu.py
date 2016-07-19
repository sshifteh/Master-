# -*- coding: utf-8 -*-
"""
Created on Tue Aug 18 14:24:50 2015

@author: shiftehs
"""
from dolfin import * 


mesh = UnitSquareMesh(3,3)
V = FunctionSpace(mesh, "Lagrange", 1)
#Define boundary conditions
u0 = Expression(" 1 + x[0]*x[0] +2*x[1]*x[1]")
def u0_boundary(x, on_boundary):
    return on_boundary
bcs = DirichletBC(V,u0,u0_boundary)

# an alternative 
#def u0_boundary(x):
#    return x[0] == 1e-8 or x[1]== 1e-8 or x[0] == 1-1e-8 or x[1] == 1-1e-8

"------------------------------------"


#Define variational problem 
u = TrialFunction(V)
v = TestFunction(V)
f = Constant(-6.0)
a = inner(grad(u),grad(v))*dx
L = f*v*dx 

"-----------------------------------"



# Compute the solution control the parameters and plot 
u_ = Function(V)
#p = parameters["krylov_solver"]
#p["absolute_tolerance"] = 1E-12
#p["relative_tolerance"] = 1E-6
#p["maximum_iterations"] = 2
#set_log_level(PROGRESS)
#set_log_level(DEBUG)

solve(a == L, u_, bcs,
    solver_parameters = dict(linear_solver= "cg", preconditioner= "ilu"))   
parameters["linear_algebra_backend"]= "STL"
#info(parameters,True)

"----------------------------------"
# grab all nodal values 
u_nodal_values = u_.vector()
u_array = u_nodal_values.array()

# test for if the value of exact solution at coord is same as u_array
#coord = mesh.coordinates()
#if mesh.num_vertices() == len(u_array):
#    for i in range(mesh.num_vertices()):
#        print "u(%8g,%8g) = %g, u0 = %g " %(coord[i][0], coord[i][1], u_array[i], u0(coord[i]))
               
#print mesh.num_vertices()

#alternativ way of calc the max error 
u_e = interpolate(u0, V)
u_e_vec = u_e.vector()
u_e_array = u_e_vec.array()
#alternativt 
u_e_array = u_e.vector().array()
error = abs(u_e_array - u_array).max()
print "Error in inf norm is :", error 
# Plot the solution 
#plot(u_)
#plot(mesh)

center = (0.5, 0.5)
print "numerical u at the center point", u_(center)
print "exact u at the center point ", u0(center )

"---------------------------------"
# normalize the solution 

print "u_array", u_array
max_u = u_array.max() #get max u_ from the array
print "max_u", max_u
u_array /= max_u
print "u_array",u_array





"---------------------------------"


#Dump the solution to file in VTK format
file = File("Poisson.pvd")
file << u_

#Hold the plot 
#interactive(True)
