from dolfin import * 

n = 50 
mesh = UnitSquareMesh(n,n)

V = VectorFunctionSpace(mesh, 'Lagrange', 2)
P = FunctionSpace(mesh, 'Lagrange', 1)
#plot(mesh, interactive = True, title = 'Mymesh')

print'No of nodes in V = span{phi_i}', V.dim() #no of basisfunc for V
print'No of nodes in P = span{phi_j}', P.dim() #no of basisfunc for P
# j is less than i, but its the same type of basisfunctions if both are p1, p1,
# piecevise linear polynomials 

W = MixedFunctionSpace([V,P])
print 'dimension of the mixed functionspace W is:', W.dim() # 27 

u, p = TrialFunctions(W)
v, q = TestFunctions(W)

f = Constant([0.0,0.0]) # This is a source for the velocity

# Lets not have any analytical terms for now 

a = inner(grad(u), grad(v))*dx + div(u)*q*dx + div(v)*p*dx
L = inner(f, v)*dx

# remember that the the Function object is the same as out function approximation u_h = sum(U_i*phi_i)
# so this is us creating a linear combination of the basis functions from W, that is W = span{phi_0, phi_1, .... , phi_27}
# What is remaining is to find the coefficients
UP = Function(W)

#print UP.dim(), eventhough this doesnt work in fenics, it is of dim 27
print' UP count:', UP.count()
U, P = UP.split(True)
#U, P = UP.split()
print 'U count: ', U.count(), ',\tP count: ', P.count()
print "P: ", P.vector()
print "U: ", U.vector()


# Making boundary conditions 
# So I want a square with zero velocity on R and L walls 
# low pressure on the top 
# and u = (0,1) at bottom, so the so the flow is going forwards


# Defining what is the L and R walls 
# Saying when x is close to 0, then we have y = all values to 1
#  Also x is clsoe to 1, then y is all values up to 1
def u_boundaryLeftRightWall(x, on_bnd):
	return  x[0] < DOLFIN_EPS or x[0] > 1- DOLFIN_EPS and on_bnd
 	

# Defining what is the inlet part
# Saying, when y is close to o(computer precision)
# then we have all x values
def u_boundaryBottom(x, on_bnd):
	return x[1] < DOLFIN_EPS and on_bnd
 
# We include on_bnd variable, which is for the internal memory of the nodes. The nodes on the boundary will then recognize that they are on bnd, and those not on bnd will not recognize it an be automatically excluded. 

# Defining what is the outlet part 
def p_boundaryTop(x, on_bnd):
	return x[1] < 1 - DOLFIN_EPS and on_bnd

#print 'testing if the subspace has the right dim: S.sub(0):',W.sub(0).dim(0) <-- not working should divide it if want it later


bcs_u_LR = DirichletBC(W.sub(0), Constant((0.0,0.0)), u_boundaryLeftRightWall)
bc_u_B  = DirichletBC(W.sub(0), Constant((0.0,1.0)),u_boundaryBottom)
bc_p_T  = DirichletBC(W.sub(1), Constant(0.0), p_boundaryTop)

bcs = [bcs_u_LR, bc_u_B, bc_p_T]


"""
# this was a bad idea because it doesnt include the bnd of pressure
# Lets try something else

bcs_up_LR = DirichletBC(W, Constant((0,0,1)), u_boundaryLeftRightWall)
bc_up_B   = DirichletBC(W, Constant((0,1,1)), u_boundaryBottom)
bc_p_T    = DirichletBC(W, Constant((0,0,0)), p_boundaryTop)

bcs = [bcs_up_LR, bc_up_B, bc_p_T]

#This didnt work ! 
"""






A, b = assemble_system(a, L, bcs)
solve(A, UP.vector(), b, "lu")

# Thiiiiiiiiiiiiiis woooooorkesssssssssssssss! Yeeeah!!! :):) 

U, P = UP.split()
plot(U, interactive = True, title="velocity , vector valued function, approximated by FE")
plot(P, interactive = True,  title="pressure, scalar valued function, appriximated by FE")

#U_analytical = project(u_analytical, V)
#P_analytical = project(p_analytical, Q)



