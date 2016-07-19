from dolfin import *
from boundary import boundary 
from stokes_solver import stokes_solver 
from numpy import * 





# Analytical wss 
def tau_wall(R, dpdz):
	return  abs(R*0.5*dpdz)	


def twoDshear(U):
	shear = 0.5*(U[0].dx(1) + U[1].dx(0)) 
	return shear


def wss(bdry, shear):
	
	shear_ = shear.vector().array()
        bdry = interpolate(bdry, shear.function_space())
	bdry_ = bdry.vector().array()

        assert len(shear_) == len(bdry_)

	# Calculating the min and max values of the shear over the boundary cells
        bdry_idx = where(bdry_ > DOLFIN_EPS)[0]
        wallshearstress = shear_[bdry_idx]
	
	print ''
        print "Max / min WSS at boundary: {} / {}".format(max(abs(wallshearstress)),
                min(abs(wallshearstress)))
	print ''

	return max(abs(wallshearstress)), min(abs(wallshearstress))


def main(x0, R, N):
	mesh = UnitSquareMesh(N,N) 
	
	# Function spaces for the approximations to the functions  
	CG1 =FunctionSpace(mesh, 'CG', 1)
	DG1 =FunctionSpace(mesh, 'DG', 1)
	DG0 = FunctionSpace(mesh, 'DG', 0)
	
	Pspace = CG1
	Vspace = VectorFunctionSpace(mesh, 'Lagrange', 2)
	Pspace = FunctionSpace(mesh, 'Lagrange', 1)
	Wspace = MixedFunctionSpace([Vspace, Pspace])
	
	Kspace = DG0
	WSSspace = DG1
	
	# The vessel geometry and interpolation into the Kspace 
	K_expr = Expression('(x[0] > (x0-R) - DOLFIN_EPS && x[0] < (x0+R) + DOLFIN_EPS) ? 0.0 : 1.0', x0 = x0, R = R)
	K = interpolate(K_expr, Kspace)

	U,P = stokes_solver(x0,R, mesh, Vspace, Pspace, Wspace, K)
	
	shear_ = twoDshear(U)


	shear = project(shear_, CG1)
	shear = interpolate(shear, Kspace)
	plot(shear, interactive = True)

	boundary_ = boundary(K, DG0)
	wss(boundary_, shear)
	
	return wss(boundary_, shear)



if __name__ == '__main__':

	mu = 0.00345 # Pascal x second 
	dpdz = -2.1
	x0 = 0.5
	R = 0.1 
	N = 40 
	
	main(x0,R,N)
	print tau_wall(R,dpdz)

	

	""" 
	Max / min WSS at boundary: 0.000122569370392 / 4.34618616546e-05
	Max / min WSS at boundary: 0.000122569370392 / 4.34618616546e-05


	tau_wall() gir : 0.105

	"""




