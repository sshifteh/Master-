# -*- coding: utf-8 -*-
"""
Created on Mon Sep  7 10:53:20 2015

@author: shiftehs
"""
# -*- coding: utf-8 -*-


from dolfin import * 
#from Subdomain_Y_boundary import RightBranch(), LeftBranch(), Stem()
#from Subdomain_Square_boundary import LeftWall(), RightWall(),Bottom(), Top() 


# Mesh and FunctionSpace 
n = 100
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




# Marking the subdomains
# Initialize mesh function for interior domains 
Ydomain = CellFunction('size_t', mesh)
Ydomain.set_all(0)
Stem().mark(Ydomain, 1)
LeftBranch().mark(Ydomain, 1)
RightBranch().mark(Ydomain, 1)

#plot(Ydomain)
#interactive()
#import sys; sys.exit()

# Initialize mesh function for boundary domains
square = FacetFunction("size_t", mesh)
square.set_all(0)
LeftWall().mark(square, 1)
RightWall().mark(square, 2)
Top().mark(square, 3)
Bottom().mark(square, 4)



plot(square)
interactive()
import sys; sys.exit()


# Dirichlet Pressure Bcs
bc_p_rightB      = DirichletBC(W.sub(1), Constant(0.0), Ydomain, 0)
bc_p_leftB       = DirichletBC(W.sub(1), Constant(0.0), Ydomain, 0)

# Dirichlet Velocity Bcs 
bc_stem        = DirichletBC(W.sub(0), Constant((0.0, 1.0)), Ydomain, 0)

bc_rightWall   = DirichletBC(W.sub(0), Constant((0.0, 0.0)), square, 2)
bc_leftWall    = DirichletBC(W.sub(0), Constant((0.0, 0.0)), square, 1)
bc_top         = DirichletBC(W.sub(0), Constant((0.0, 0.0)), square, 3)
bc_bottom      = DirichletBC(W.sub(0), Constant((0.0, 0.0)), square, 4)

bcs = [bc_p_rightB, bc_p_leftB,bc_stem, bc_rightWall, bc_leftWall, bc_top, bc_bottom]

f= project(Constant((0, 0, 1)), W, bcs) # (x,y,z)
plot(f, interactive=True, title = 'Dirichlet boundary conditions')



"""
# Test and Trial functions 

u,p = TrialFunctions(W)
v,q = TestFunctions(W)


# Variational formulation 

a0 = Constant((0.0, 0.0))     # inside Y 
a1 = Constant((1000, 1000))   # outside

a = (inner(a0*grad(u), grad(v))*dx(0) + inner(a1*grad(u), grad(v))*dx(1)
L = inner(Constant((0.0))*v)*ds(0) - inner(Constant((0.0, 0.0))*v)*ds(1)


#w_ = Function(W)
#solve(a == L, w_, bcs )
#plot(w_, interactive = True, title = 'Solution')

#u_, p_ = w_.split()
#File("velocity.xdmf") << u_
#File("pressure.xdmf") << p_
"""