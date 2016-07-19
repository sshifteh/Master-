# -*- coding: utf-8 -*-
"""
Created on Mon Sep  7 12:53:48 2015
@author: shiftehs
"""

from dolfin import *

N = 40
w = 0.25

mesh = UnitSquareMesh(N, N)
subdomains = CellFunction('size_t', mesh, 0)
subdomains.set_all(0)

# The CellFunction gives each cell the value 0 
# When we make subdomains we loop over all cells and assign them a different value 


class Stem(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and x[1] < 0.5 + DOLFIN_EPS and \
               0.5-w-DOLFIN_EPS < x[0] < 0.5+w+DOLFIN_EPS 

class LeftBranch(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and x[1] > 0.5 - DOLFIN_EPS and x[0] < 0.5 + DOLFIN_EPS\
               and not (x[1] < -x[0] + 1 - w)\
               and not (x[1] > -x[0] + 1 + w)

class RightBranch(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and x[1] > 0.5 - DOLFIN_EPS and x[0] > 0.5 - DOLFIN_EPS\
               and not (x[1] < x[0] - w)\
               and not (x[1] > x[0] + w)

Stem().mark(subdomains, 1)
LeftBranch().mark(subdomains, 1)
RightBranch().mark(subdomains, 1)


plot(subdomains, interactive = True, title = 'Y-shape')
import sys; sys.exit()


"""
# The Square and bc on the square

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

# Dirichlet Pressure Bcs
bc_p_rightB      = DirichletBC(W.sub(1), Constant(0.0), subdomains, 1)
bc_p_leftB       = DirichletBC(W.sub(1), Constant(0.0), subdomains, 1)

# Dirichlet Velocity Bcs 
bc_stem        = DirichletBC(W.sub(0), Constant((0.0, 1.0)), subdomains, 0)

bc_rightWall   = DirichletBC(W.sub(0), Constant((0.0, 0.0)), RightWall())
bc_leftWall    = DirichletBC(W.sub(0), Constant((0.0, 0.0)), LeftWall())
bc_top         = DirichletBC(W.sub(0), Constant((0.0, 0.0)), Top())
bc_bottom      = DirichletBC(W.sub(0), Constant((0.0, 0.0)), Bottom())

bcs = [bc_p_rightB, bc_p_leftB,bc_stem, bc_rightWall, bc_leftWall, bc_top, bc_bottom]

f= project(Constant((0, 0, 1)), W, bcs) # (x,y,z)
plot(f, interactive=True, title = 'Dirichlet boundary conditions')
""" 