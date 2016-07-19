from alpha import alpha
from dolfin import *

side_width = 0.1

# Analytical velocity - slow Python version  
class MyExpression(Expression):
	def __init__(self, x0 = 0.5, R = 0.05, nu = 0.00345, dpdz = -2.1):
		self.x0 = x0
		self.R = R
		self.nu = nu
		self.dpdz = dpdz

	def eval(self, values, x): 
		if  self.x0 - self.R <= x[0] <= self.x0 + self.R:
			# inside pipe
			values[0] = 0			
			values[1] = 0.25*self.nu*self.dpdz*((x[0]-self.x0)*(x[0]-self.x0) - self.R*self.R)
		else:
			# outside pipe 
			values[0] = 0
			values[1] = 0

	# Spesify the value shape of the expression 
	def value_shape(self):
		return (2,)

#V = VectorFunctionSpace(mesh2, 'Lagrange', 2)
#g = MyExpression(R = 0.25) 
mesh2 = UnitSquareMesh(10,10)
#plot(g, mesh = mesh2)
#interactive()


# Faster C++ version    
u_analytical = '''  
class AnalyticalVelocity : public Expression
{



	public:
		double x0, R, nu, dpdz;

        void eval(Array<double>& values, const Array<double>& x) const
	{
		if(((x0 - R) <= x[0])   &&  (x[0] <= (x0 + R)))
		{  
			double velocity = 0.25*nu*dpdz*((x[0]-x0)*(x[0]-x0) - R*R); 		
			values[0] = 0.0;
                        values[1] = velocity; 
		}
		else 
		{
			values[0] = 0.0; 
                        values[1] = 0.0;
		}
	}

        size_t value_rank() const{
            return 1;
        }

        size_t value_dimension(size_t i) const{
            return 2;
        }

};''' 


# Create a fine mesh for comparing analytical and numerical solution
#V = VectorFunctionSpace(mesh2, 'Lagrange', 2)
#g = MyExpression(R = 0.25) 
fine_mesh = UnitSquareMesh(128,128)
fine_V = VectorFunctionSpace(fine_mesh, 'Lagrange', 2)


nu = 0.00345 # Pascal x second 
dpdz = -2.1
x0 = 0.50
#R = 0.25	
R = side_width



u_c = Expression(u_analytical) 
u_c.x0 = float(x0);
u_c.R  = float(R);
u_c.nu = float(nu);
u_c.dpdz = float(dpdz);

print u_c.value_rank()

# Experimental solver and error 

def main(w, N):  
	"""
	Solves the Stokes equation with an additional inverse permeability term. 
	Calculates the error between analytical and experimental velocity.	
	"""
	
	mesh = 	UnitSquareMesh(N,N)
	Vspace = VectorFunctionSpace(mesh, 'Lagrange', 2)
	Pspace = FunctionSpace(mesh, 'Lagrange', 1)
	Wspace = MixedFunctionSpace([Vspace, Pspace]) 
	Kspace = FunctionSpace(mesh, 'DG', 0)

	# The vessel geometry and interpolation into the Kspace 
	K_expr = Expression('(x[0] > w - DOLFIN_EPS && x[0] < 1 - w + DOLFIN_EPS) ? 0.0 : 1.0', w=w)
	K = interpolate(K_expr, Kspace)	

	# Test and trial spaces over the Taylor-Hood elements  
	u, p = TrialFunctions(Wspace)
	v, q = TestFunctions(Wspace)

	# Source term
	f = Constant([0.0,0.0])
	
	# Parameters 
	nu = Constant(0.00345) # Pascal x second 
	dpdz = Constant(-2.1)
	R = w #Constant(0.05)

	# Variational form 
 	a = inner(alpha(u, K), v)*dx + nu*inner(grad(u), grad(v))*dx + div(u)*q*dx + div(v)*p*dx
	L = inner(f, v)*dx
	UP = Function(Wspace)

	# Parabolic inlet velocty Dirichlet condition 
        speed = 10
	#velocityFunc = Expression(["0","-speed*(x[0]-0.45)*(x[0]-0.55)"], speed=speed)
	analyticalVelocity = Expression(["0"," 0.25*nu*dpdz*((x[0]-0.50)*(x[0]-0.50) - R*R) "], speed=speed, 		nu=nu, R=R, dpdz=dpdz)


	# Defining the boundaries of the vessel 
	def u_boundaryBottom(x, on_bnd):
		return x[1] < DOLFIN_EPS and x[0] > w - DOLFIN_EPS and x[0] < 1-w + DOLFIN_EPS and on_bnd

	def restBottom(x, on_bnd):

		return x[1] < DOLFIN_EPS and  not (x[0]> w -DOLFIN_EPS) and not (x[0] < 1-w + DOLFIN_EPS) and 			on_bnd

    	def leftWall(x, on_boundary):
        	return near(x[0], 0.0)

	def rightWall(x, on_boundary):
        	return near(x[0], 1.0)


	def top_rest(x, on_bnd):
		return x[1] < 1 - DOLFIN_EPS and not (x[0] > w) and not (x[0] < 1-w) and on_bnd

	inlet  = DirichletBC(Wspace.sub(0), analyticalVelocity, u_boundaryBottom)
	restBottom = DirichletBC(Wspace.sub(0), Constant((0.0,0.0)), restBottom)
	
	right = DirichletBC(Wspace.sub(0), Constant((0.0,0.0)), rightWall)
	left = DirichletBC(Wspace.sub(0), Constant((0.0,0.0)), leftWall)
	top = DirichletBC(Wspace.sub(0), Constant((0.0, 0.0)), top_rest )
	bcs = [inlet, restBottom, right, left, top] 

	# Solving the system of linear equations 
	A, b = assemble_system(a, L, bcs)
        solver = LUSolver(A)
	solver.solve(UP.vector(), b) #instead of saving the solution in U immediatly to plot it successively
	U, P = UP.split()

	# Hack to save time perhaps, interpolate and use degree rise 0 
	#V = FunctionSpace(mesh, 'Lagrange', 2)
	#V1 = VectorFunctionSpace(mesh, 'Lagrange', 2)
	
	#u_e = interpolate(u_c, V1)#Wspace.sub(0))

	#u_e_proj = project(u_e, Wspace.sub(0))  
	#bcs.apply(u_e.vector())		
	#plot(u_c, mesh, title = 'Exact Velocity Field')
	#plot(U,mesh, title = 'Experimental Velocity Field')
	print 'success inside main '
        #interactive()
	#u_e = ()
        
        error  = interpolate(U, fine_V)
	error.vector().axpy(-1, interpolate(u_c, fine_V).vector())
        plot(error, title=str(N))
	#from IPython import embed; embed()

	plot(interpolate(U, fine_V),mesh, title = 'Experimental Velocity Field')
	plot(interpolate(u_c, fine_V),mesh, title = 'Analytical Velocity Field')
	interactive()
	u_error = errornorm(u_c, U, norm_type='L2', degree_rise=4, mesh=mesh)
	#l2 is the vector norm, while L2 is the norm for the function space

	#u_error = errornorm(u_c, U, norm_type = )#, degree_rise = 4)
        
	return u_error, mesh.hmin()



E = []; h = []; 
for n in [64]:
	ei, hi = main( w = side_width, N = n )
	E.append(ei)
	h.append(hi)
	

from math import log as ln 
for i in range(1, len(E)):
	r = ln( E[i] / E[i - 1]) / ln(h[i] / h[i - 1])
	print " h = % 2.2E  E = % 2.4E   r = %.2f"  % ( h[i] , E[i] , r )
	


"""
comments:
	If the error is already zero there is nothing to converge to. 
	Also make a mistake in the boundary velocity and see what happens.
	Also they are both approximated by p2 functions so it will be exact. 


Questions left to ansver:
	does it become stabile after a while. doesnt change even if we change the mesh? 
	

"""

