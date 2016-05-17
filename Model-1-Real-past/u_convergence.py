from alpha import alpha #, C
from dolfin import *
import numpy as np

	
# Analytical velocity - slow Python version  
class MyExpression(Expression):
	def __init__(self, x0 = 0.5, R = 0.05, nu = 0.00345, dpdz = -2.1):
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
		self.dpdz = dpdz

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
			values[1] = 0.25*self.nu*self.dpdz*(-(x[0]-self.x0)*(x[0]-self.x0) - self.R*self.R)
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
		double x0, R, nu, dpdz;

	// The eval function for the function expression 
        void eval(Array<double>& values, const Array<double>& x) const
	{
		if(((x0 - R - DOLFIN_EPS) <= x[0])   &&  (x[0] <= (x0 + R + DOLFIN_EPS)))
		{  
			double velocity = 0.25*nu*dpdz*(-(x[0]-x0)*(x[0]-x0) + R*R); 		
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












def stokes_solver(nu, dpdz, N, x0, R, viz = False ):  


	print 'solving for %s' %N
	
	# THE MESH AND FUNCTION SPACES OVER THE MESH 
	mesh = 	UnitSquareMesh(N,N)
	h = 1./N 
		

	V = VectorFunctionSpace(mesh, 'Lagrange', 2) 
	P = FunctionSpace(mesh, 'Lagrange', 1)       
	W = MixedFunctionSpace([V, P]) 
	Kspace = FunctionSpace(mesh, 'DG', 0)



	u, p = TrialFunctions(W)
	v, q = TestFunctions(W)
	f = Constant([0.0,0.0])




	
	K_expr = Expression('( x[0] > (x0 - R - DOLFIN_EPS) && x[0] < (x0 + R + DOLFIN_EPS)) ? 0.0 : 1.0', x0 = x0, R = R )
	K = interpolate(K_expr, P)	 # annen mulighet er Kspace
	
	

 	a = inner(alpha(u, K), v)*dx + nu*inner(grad(u), grad(v))*dx + div(u)*q*dx + div(v)*p*dx
	L = inner(f, v)*dx
	UP = Function(W)




	def u_boundaryBottom(x, on_bnd):
		return x[1] < DOLFIN_EPS and x[0] > (x0 - R - DOLFIN_EPS) and x[0] < (x0 + R + DOLFIN_EPS) #and on_bnd

	

	analyticalVelocity = Expression(["0","0.25*nu*dpdz*(-(x[0]-x0)*(x[0]-x0) + R*R) "], nu=nu, x0=x0, R=R, dpdz=dpdz)
	inlet  = DirichletBC(W.sub(0), analyticalVelocity, u_boundaryBottom)
	bcs = [inlet] 



	
	A, b = assemble_system(a, L, bcs)
        solver = LUSolver(A)
	solver.solve(UP.vector(), b) 
	U, P = UP.split()  

	u = interpolate(U, V)

	#plot(K, interactive=viz, title = 'K')
	#plot(P, interactive=viz, title = 'pressure')
        #plot(U, interactive=viz, title = 'velocity')
	


	# CALCULATING THE ERROR/DISCREPANCY FIELD
	u_fineMesh = interpolate(U, fine_V)
	u_exact_fineMesh = interpolate(u_c, fine_V)
	#plot(u_exact_fineMesh, interactive = viz, title = 'eksakt paa fine mesh')
	
	


	# Make an available vector in the same space to hold the error feld 
	error_fineMesh = Function(fine_V)	
	# The error field is the discrepancy between u_exact - u 
	error_fineMesh.vector()[:] = u_exact_fineMesh.vector()[:] - u_fineMesh.vector()[:] 
	#plot(error_fineMesh, interactive = viz, title = 'e = u_e - u')
	
	#File("U_C.xdmf") <<  U
	#File("K_C.xdmf") <<  K


	# ***************** CALCULATING THE ERRORNORM *****************
	
	# 1) EXPLICITLY  

	# might have to do something with Martins analytical formula
	# can use 6 and 8 and se if it makes a difference? 
	# then in the error I think  i might see a difference in the decimals 	
	q = 2 
	degree = 2*(q + 2)	
	error_norm = sqrt(assemble((u_c - U)**2*dx(degree =degree)))
	

	
	# 2 ) fra biblioteket.  funker ikke 
	# bruker errornormen fra FEniCS biblioteket 
	errornorm_u = errornorm(u_c, U, norm_type="H1", degree_rise = 2, mesh = None)
	#l2 is the vector norm, while L2 is the norm for the function space
        


	# 3) selvlagd H1 errornorm. gir feil meld 
	#error = sqrt(assemble( (error_fineMesh)**2 + (grad(error_fineMesh))**2)*dx)
	


	return u, errornorm_u, h 
	


if __name__ == '__main__':
	import numpy as np 

	# DEFINING THE MESH PARAMETERS AND FLUID FLOW PARAMETERS 
	center = 0.50      # CENTER OF THE TUBE 
	side_width = 0.1   # RADIUS OUT FROM THE CENTER 
	nu = 1.0 #0.00345       # DYNAMIC FLUID VISCOSITY, MATERIAL PROPERTY  
	dpdz = 2.1        # PRESSURE GRADIENT OF THE FLUID FLOW 
	#side_width = 0.5   # RADIUS OUT FROM THE CENTER, then it is the entire domain  



	# SENDING PARAMETERS TO C++ FUNCTION 
	u_c = Expression(u_analytical) 
	u_c.x0 = float(center);
	u_c.R  = float(side_width);
	u_c.nu = float(nu);
	u_c.dpdz = float(dpdz);


	# FINE MESH FOR ANALYTICAL SOLUTION AND COMPARISON BETWEEN ANALYTICAL AND EXPERIMENTAL SOLUTION  	
	fine_mesh = UnitSquareMesh(128,128)
	fine_V = VectorFunctionSpace(fine_mesh, 'Lagrange', 2)
	


	e = [] ; h_ = [];
	for i in [16,32,64, 128]: # 128
		u, error, h = stokes_solver(nu = nu, dpdz = dpdz, N = i, x0 = center, R = side_width )
		print 'u max', max(u.vector())
		e.append(error)
		h_.append(h)	

	# formel 1 
	#k_verdier = []
	#for i in range(1, len(e)):
	#	k = np.log(e[i-1]/e[i]) /np.log(2.)
	#	k_verdier.append(k)
	# gav ikke 3 

	# formel 2

	print 'Mesh size, h=1/N          Error  '
	for i, j in zip(h_, e):
		print '%5.3f        %20.5e' %(i,j) 
	
	#print k_verdier 

	
	for i in range(1, len(e)):
		r = np.log(e[i]/e[i-1]) / np.log(h_[i]/h_[i-1]) 
		print 'r=', r  







	
	# PLOTTING A LOGARITMIC PLOT OF THE ERROR AS A FUNCTION OF THE MESH SIZE 
	import pylab as pylab 
	enable_plot = True
	if enable_plot:
		pylab.loglog(h_, e, 'g-o' ) 
		#pylab.plot(h_values, error_values, 'b-o')	
		pylab.xlabel('h')
		pylab.ylabel('error')
		pylab.grid(True)
		pylab.savefig('loglog_velocity')		
		pylab.show()		
		
		
	

	

	




