# -*- coding: utf-8 -*-
"""
Created on Wed Aug 19 17:36:29 2015

@author: shiftehs
"""
from dolfin import * 

def boundary(x, on_boundary):
    tol = 1e-8
    return on_boundary and abs(x[0] < tol)
    
# Alternatively make an instans of object of a subclass of SubDomain
# The inside method is an alternative to the stand alone function above

class Boundary(SubDomain):
    def inside(self, x, on_boundary):
        tol = 1e-8
        return on_boundary and abs(x[0] < tol)

boundary_ = Boundary() # Make and instans of the object 

bcs = DirichletBC(V, Constant(0.0), boundary_) 

class Omega0(SubDomain):
    def inside(self, x, on_boundary):
        return True if x[1] <= x[0]*x[0] else False 
class Omega1(SubDomain):
    def inside(self, x, on_boundary):
        return True if x[1] >= x[0]*x[0] else False 
    
n = 5    
mesh = UnitSquareMesh(n,n)    
subdomains = MeshFunction("uint", mesh,2)

subdomain0 = Omega0()
subdomain0.mark(subdomains, 0)

subdomain1 = Omega1()
subdomain1.mark(subdomains, 1)

V0 = FunctionSpace(mesh, "DG",1)
k = Function(V0) # en funksjon som er constant over hvert subdomene 
# naa skal vi fylle k med den rette verdien for hvert domene 
k_values = [1.5, 10]
for cell_no in range(len(subdomains.array()))



#subdomain.array() returnerer subdomain valuene som et array
# subdomain.array()[i] gir verdien i celle i i subdomenet  




