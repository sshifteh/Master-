# -*- coding: utf-8 -*-
"""
Created on Mon Sep  7 14:21:53 2015

@author: shiftehs
"""

from dolfin import *


N = 40
w = 0.25

mesh = UnitSquareMesh(N, N)


class Stem(SubDomain):
    def __init__(self, only_boundary=False):
        self.only_boundary = only_boundary
        SubDomain.__init__(self)
        
    def inside(self, x, on_boundary):
        inside = x[1] < 0.5 + DOLFIN_EPS and \
               0.5-w-DOLFIN_EPS < x[0] < 0.5+w+DOLFIN_EPS
               
        if self.only_boundary:
            return inside and on_boundary
        else:
            return inside
             

class LeftBranch(SubDomain):
    def inside(self, x, on_boundary):
        return x[1] > 0.5 - DOLFIN_EPS and x[0] < 0.5 + DOLFIN_EPS\
               and not (x[1] < -x[0] + 1 - w)\
               and not (x[1] > -x[0] + 1 + w)

class RightBranch(SubDomain):
    def inside(self, x, on_boundary):
        return x[1] > 0.5 - DOLFIN_EPS and x[0] > 0.5 - DOLFIN_EPS\
               and not (x[1] < x[0] - w)\
               and not (x[1] > x[0] + w)

subdomains = CellFunction('size_t', mesh, 1)
subdomains.set_all(1)
# The CellFunction gives each cell the value 0 
# When we make subdomains we loop through the mesh and assign values  
Stem().mark(subdomains, 0)
LeftBranch().mark(subdomains, 0)
RightBranch().mark(subdomains, 0)

plot(subdomains, interactive = True, title = 'Y-shape')




class LeftWall(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[0], 0.0)

class RightWall(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[0], 1.0)

class Bottom(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[1], 0.0) and not Stem().inside(x, on_boundary)

class Top(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[1], 1.0) and not RightBranch().inside(x, on_boundary) and not LeftBranch().inside(x, on_boundary) 


V = VectorFunctionSpace(mesh, 'CG', 2)
Q = FunctionSpace(mesh, 'CG', 1)
W = V*Q

# Convert the subdomains into a dolfin.Function such we can use it 
# in the variational form
S = FunctionSpace(mesh, 'DG', 0)
s = Function(S)
s.vector()[:] = subdomains.array()
#plot(s, interactive=True)
# TODO: Use the s function in your Stokes variational form

# Dirichlet Pressure Bcs
bc_p_rightB      = DirichletBC(W.sub(1), Constant(0.0), LeftBranch())
bc_p_leftB       = DirichletBC(W.sub(1), Constant(0.0), RightBranch())

# Dirichlet Velocity Bcs 
bc_stem        = DirichletBC(W.sub(0), Constant((0.0, 1.0)), Stem(only_boundary=True))

bc_rightWall   = DirichletBC(W.sub(0), Constant((0.0, 0.0)), RightWall())
bc_leftWall    = DirichletBC(W.sub(0), Constant((0.0, 0.0)), LeftWall())
bc_top         = DirichletBC(W.sub(0), Constant((0.0, 0.0)), Top())
bc_bottom      = DirichletBC(W.sub(0), Constant((0.0, 0.0)), Bottom())

bcs = [bc_p_rightB, bc_p_leftB,bc_stem, bc_rightWall, bc_leftWall, bc_top, bc_bottom]

f= project(Constant((0, 0, 1)), W, bcs) # (x,y,z)
plot(f, interactive=True, title = 'Dirichlet boundary conditions')

