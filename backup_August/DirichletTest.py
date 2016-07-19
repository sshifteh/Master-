# -*- coding: utf-8 -*-
"""
Created on Fri Sep  4 14:19:28 2015
@author: shiftehs

Boundary betingelser, Dirichlet paa enehets firekanten,
for at fluidet skal stromme nedover 

p = 0 paa grensene for at fluidet skal trekkes fra omraade hoyt trykk og ned 
dit (omraade med lavt trykk)

"""

from dolfin import * 
#from subdomain_Y import LeftBranch,RightBranch, Stem 


n = 50
w = 0.2
mesh = UnitSquareMesh(n, n)

V = VectorFunctionSpace(mesh, 'CG', 2)
Q = FunctionSpace(mesh, 'CG', 1)
W = V*Q

class RightBranch(SubDomain):
    """
    
    The class RightBranch inherits from SubDomain which is of object type CellFuntion.
    SubDomain assigns a value to the midpoint of the cell. 
    RightBranch defines the subdomain in the inside method. 

    on_boundary is included in the return statement to avoid including interior points.
    
    """
    def inside(self, x, on_boundary):
        return  on_boundary and x[1] > 0.5 - DOLFIN_EPS and x[0] > 0.5 - DOLFIN_EPS\
               and not (x[1] < x[0] - w)\
               and not (x[1] > x[0] + w)


class LeftBranch(SubDomain):
    """
    LeftBranch inherits functionality of the class Subdomain which is a CellFunction
    i.e. SubDomain assigns value to the midpoint of cells. 
    Here the LeftBranch point is defined as one subdomain. 
    
    on_boundary is included in the return statement so the velocity is only on the 
    boundary and not in interior cells points. 
    
    """
    def inside(self, x, on_boundary):
        return on_boundary and x[1] > 0.5 - DOLFIN_EPS and x[0] < 0.5 + DOLFIN_EPS\
               and not (x[1] < -x[0] + 1 - w)\
               and not (x[1] > -x[0] + 1 + w)
               
class Stem(SubDomain):
    """

    The Stem 
    
    """
    def inside(self, x, on_boundary):
        return on_boundary and x[1] < 0.5 + DOLFIN_EPS and \
                   0.5-w-DOLFIN_EPS < x[0] < 0.5+w+DOLFIN_EPS 

class LeftWall(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[0], 0.0)

class RightWall(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[0], 1.0)

class Bottom(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[1], 0.0)

class Top(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[1], 1.0)



# velocity bcs 
bc_rightB      = DirichletBC(W.sub(0), Constant((0.0, 1.0)), RightBranch())
bc_leftB       = DirichletBC(W.sub(0), Constant((0.0, 1.0)), LeftBranch())
bc_stem        = DirichletBC(W.sub(0), Constant((0.0, 1.0)), Stem())

bc_rightWall   = DirichletBC(W.sub(0), Constant((0.0, 0.0)), RightWall())
bc_leftWall    = DirichletBC(W.sub(0), Constant((0.0, 0.0)), LeftWall())
bc_top         = DirichletBC(W.sub(0), Constant((0.0, 0.0)), Top())
bc_bottom      = DirichletBC(W.sub(0), Constant((0.0, 0.0)), Bottom())

bcs = [bc_rightB, bc_leftB, bc_stem, bc_rightWall, bc_leftWall, bc_top, bc_bottom]

f= project(Constant((0, 0, 0)), W, bcs)
plot(f, interactive=True, title = 'Dirichlet boundary conditions' )