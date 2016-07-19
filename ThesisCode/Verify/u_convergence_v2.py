from alpha import alpha
from dolfin import *


	
# Analytical velocity - slow Python version  
class MyExpression(Expression):
	def __init__(self, x0 = 0.5, R = 0.05, nu = 0.00345, dpdy = -2.1):
		"""
		Defining the parameters in the constructor.
		x0 : the center of the pipe 
		R  : the radius our of the center of the pipe 
		nu : dynamic viscosity 
		dpdz: the pressure gradient along the pipe 
		""" 			
		self.x0 = x0
		self.R = R
		self.nu = nu
		self.dpdy = dpdy

	def eval(self, values, x): 
		"""
		Method for the analyrical expression of the velocity in a cylindrical pipe.
		The velocity profile doesnt change in the x-direction, so x = 0.
		The velocity profile is always parabolic in the y-direction and is 
		
		y = 0.25*nu*dpdz ((x-x0)**2 - R**2)		
	

		"""
		if  self.x0 - self.R <= x[0] <= self.x0 + self.R:
			# inside pipe
			values[0] = 0			
			values[1] = 0.25*self.nu*self.dpdy*((x[0]-self.x0)*(x[0]-self.x0) - self.R*self.R)
		else:
			# outside pipe 
			values[0] = 0
			values[1] = 0

	def value_shape(self):
		"""	
		Method for specifying that we are in two dimensional space. 
		"""
		return (2,)






# Faster C++ version    
u_analytical = '''  
class AnalyticalVelocity : public Expression
{
	// Parameters 
	public:
		double x0, R, mu, dpdy;

	// The eval function for the function expression 
        void eval(Array<double>& values, const Array<double>& x) const
	{


		if(  ((x0 - R) <= x[0]) && (x[0] <= (x0+R)))
		{  
			double velocity = 0.25*mu*dpdy*(-(x[0])*(x[0]) + R*R); 		
			values[0] = 0.0;
                        values[1] = velocity; 
		}
		else 
		{
			values[0] = 0.0; 
                        values[1] = 0.0;
		}
	}
	// Making the array 2 dimensional 
        size_t value_rank() const{
            return 1;
        }

        size_t value_dimension(size_t i) const{
            return 2;
        }

};''' 





 
def stokes_solver(mu, dpdy, N, x0, R):  


	print 'solving for %s' %N

	mesh = RectangleMesh(Point(-1.0,-1.0) ,Point( 1.0, 1.0), N, N, "right")
	h = 1./N # The element size
		
	Vspace = VectorFunctionSpace(mesh, 'Lagrange', 2) 
	Pspace = FunctionSpace(mesh, 'Lagrange', 1)       
	Wspace = MixedFunctionSpace([Vspace, Pspace]) 
	Kspace = FunctionSpace(mesh, 'DG', 0)

	u, p = TrialFunctions(Wspace)
	v, q = TestFunctions(Wspace)
	f = Constant([0.0,0.0])

	
	phi_I = Expression('( x[0] > (x0 - R - DOLFIN_EPS) && x[0] < (x0 + R + DOLFIN_EPS)) ? 0.0 : 1.0', x0 = x0, R = R )
	K = interpolate(phi_I, Pspace)	
	#plot(K, interactive = True, title = 'K')
	#File('K.xdmf') << K


 	a = inner(alpha(u, K), v)*dx + mu*inner(grad(u), grad(v))*dx + div(u)*q*dx + div(v)*p*dx
	L = inner(f, v)*dx
	UP = Function(Wspace)


	# BOUNDARY CONDITION
	# Only need bcs on the inlet; surrounding homegenous Dirichlet bcs will result in fenics not finding facets mathing domain for bcs   
	# DEFINING PHYSICAL BOUNDARIES 
	def inlet_boundary(x, on_bnd):
		return (x[1] - (-1)) < DOLFIN_EPS  and x[0] > - R  and x[0] <  R  
  

	# ASSIGNING THE PHYSICAL BOUNDARY CONDITIONS  
	# ANALYTICAL SOLUTION IS INLET BCS (grenseflatebetingelse - dynamisk)  
	analyticalVelocity = Expression(["0","0.25*mu*dpdy*(-(x[0])*(x[0]) + R*R) "], mu=mu, R=R, dpdy=dpdy)
	inlet  = DirichletBC(Wspace.sub(0), analyticalVelocity, inlet_boundary)
	bcs = [inlet] 


	A, b = assemble_system(a, L, bcs)
        solver = LUSolver(A)
	solver.solve(UP.vector(), b) 
	U, P = UP.split()  
	
       
	# **************  CALCULATING THE ERROR/DISCREPANCY FIELD  *********************
	

	# COMPARISON ON THE FINE MESH; 128X1278 

	u_experimental_onFineMesh = interpolate(U, fine_V)
	u_exact_onFineMesh = interpolate(u_c, fine_V)
	#plot(u_exact_onFineMesh, interactive = True, title = 'exact')

	# Make an available vector in the same space to hold the error feld 
	error_onFineMesh = Function(fine_V)	
	#The error field is the discrepancy between u_exact - u 
	error_onFineMesh.vector()[:] = u_exact_onFineMesh.vector()[:] - u_experimental_onFineMesh.vector()[:] 
	#plot(error_onFineMesh, interactive = True, title = 'discrepancy field')
	


	
	#File('u_analytisk.xdmf') << u_exact_onFineMesh
	#File('u_discrepancy_N{}.xdmf'.format( N)) << error_onFineMesh



	# ***************** CALCULATING THE ERRORNORM *****************
	# 2) USING ERRORNORM METHOD FROM THE FENICS LIBRARY 
	u_error = errornorm(u_c, U, norm_type='H1', degree_rise=4, mesh=mesh)
	#l2 is the vector norm, while L2 is the norm for the function space
      



	return u_error, mesh.hmin(), h








if __name__ == '__main__':


	# DEFINING THE MESH PARAMETERS AND FLUID FLOW PARAMETERS 
	center = DOLFIN_EPS       # Center 
	radius = 0.2              # Radius
	mu = 1.0                  # Dynamic viscosity, MATERIAL PROPERTY  
	dpdy = 2.1                # Pressure gradient of the blood flow 
	#radius = 1.0         # then it is the entire domain  



	# SENDING PARAMETERS TO C++ FUNCTION 
	u_c = Expression(u_analytical)        # Analytical solution written in C++
	u_c.x0 = float(center);
	u_c.R  = float(radius);
	u_c.mu = float(mu);
	u_c.dpdy = float(dpdy);



	# FINE MESH FOR ANALYTICAL SOLUTION AND COMPARISON BETWEEN ANALYTICAL AND EXPERIMENTAL SOLUTION  	
	fine_mesh = RectangleMesh(Point(-1.0,-1.0) ,Point( 1.0, 1.0), 128, 128, "right")	
	fine_V = VectorFunctionSpace(fine_mesh, 'Lagrange', 2)
	

	elist = []
	hminlist =[]
	hlist = []
	for n in [16,32,64,128]:
		e, hmin, h = stokes_solver(mu = mu, dpdy = dpdy, N = n, x0 = center, R = radius) 
		print e, hmin, h












	"""
	C = 1E5 
	error_values = []	
	h_values = []   # mesh str er ikke avhengig av C. det er den ene og samme i alle iterasjonene

	for n in [16,32,64,128]:  	
		error_i, h_i = stokes_solver(C=C, mu = mu, dpdy = dpdy, N = n, x0 = center, R = radius) 
		error_values.append(error_i)
		h_values.append(h_i)
	
	
	
	print '  h=1/N          H1 errornorm  '
	for i, j in zip(h_values, error_values):
		print '%5.3f        %10.5e' %(i,j) 

	""" 	
	

	

	

	




