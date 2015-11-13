# -*- coding: utf-8 -*-
"""
Created on Fri Sep  4 15:27:17 2015
@author: shiftehs
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
    Subclass of SubDomain.
    Object type CellFunction    
    """
    def inside(self, x, on_boundary):
        return  on_boundary and x[1] > 0.5 - DOLFIN_EPS and x[0] > 0.5 - DOLFIN_EPS\
               and not (x[1] < x[0] - w)\
               and not (x[1] > x[0] + w)


class LeftBranch(SubDomain):
    """   
    Subclass of SubDomain.
    Object type CellFunction    
    """
    def inside(self, x, on_boundary):
        return on_boundary and x[1] > 0.5 - DOLFIN_EPS and x[0] < 0.5 + DOLFIN_EPS\
               and not (x[1] < -x[0] + 1 - w)\
               and not (x[1] > -x[0] + 1 + w)
               
class Stem(SubDomain):
    """    
    Subclass of SubDomain.
    Object type CellFunction    
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
        return near(x[1], 0.0) and not Stem().inside(x, on_boundary)

class Top(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[1], 1.0) and not RightBranch().inside(x, on_boundary) and not LeftBranch().inside(x, on_boundary) 




# Pressure Bcs
 
bc_p_rightB      = DirichletBC(W.sub(1), Constant(0.0), RightBranch())
bc_p_leftB       = DirichletBC(W.sub(1), Constant(0.0), LeftBranch())

# Velocity bcs 

bc_stem        = DirichletBC(W.sub(0), Constant((0.0, 1.0)), Stem())

bc_rightWall   = DirichletBC(W.sub(0), Constant((0.0, 0.0)), RightWall())
bc_leftWall    = DirichletBC(W.sub(0), Constant((0.0, 0.0)), LeftWall())
bc_top         = DirichletBC(W.sub(0), Constant((0.0, 0.0)), Top())
bc_bottom      = DirichletBC(W.sub(0), Constant((0.0, 0.0)), Bottom())

bcs = [bc_p_rightB, bc_p_leftB,bc_stem, bc_rightWall, bc_leftWall, bc_top, bc_bottom]

#f= project(Constant((0, 0, 1)), W, bcs) # (x,y,z)
#plot(f, interactive=True, title = 'Dirichlet boundary conditions')


# Test and Trial functions 

u,p = TrialFunctions(W)
v,q = TestFunctions(W)

# Variational formulation 

a = inner(grad(u), grad(v))*dx + div(v)*p*dx + q*div(u)*dx
L = inner(Constant((0.0,0.0)), v)*dx
w_ = Function(W)
solve(a == L, w_, bcs )
#plot(w_, interactive = True, title = 'Solution')

u_, p_ = w_.split()
File("velocity.xdmf") << u_
File("pressure.xdmf") << p_