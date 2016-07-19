from alpha import alpha
from dolfin import *


	
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
			values[1] = 0.25*self.nu*self.dpdz*((x[0]-self.x0)*(x[0]-self.x0) - self.R*self.R)
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
	// Making the array 2 dimensional 
        size_t value_rank() const{
            return 1;
        }

        size_t value_dimension(size_t i) const{
            return 2;
        }

};''' 





 
def experimentalSolver(nu,dpdz, N,x0 ,R ):  

	"""
	************************** THE STONG FORMULATION OF THE EQUATIONS  ***************************** 
 


	Solves the Stokes equations with inverse permeability,
 
		nu*alpha(u,k) - laplace(u) + del(p)  = f in Omega, and
		
		div(u)  = 0 ,


	the first equation is termed the momentum equation and the second the continuity equation,
	where alpha is defined as 
	
		alpha = CKu

	where C is a material constant defined as 1E15. 
	K is the interface (phi in the notes) defined as,

		K = 0       if x0-R < (x, y) < x0+R  then defined as point in the fluid domain 
		    1 	    else point in the solid domain  

		
	where x0 and  R are parameters defining the geomtry, where x0 is the center point about which the solution is symmetric,
 	chosen as    
		
		x0 = 0.5, 

	and R is the radius extending from the center point, chosen as, 

		R = 0.1,  

	f is a body force (gravity or magnetic force if charged particles(iones) in the fluid) here it is defined as,
		f = 0, 

	nu is the dynamic viscosity (not kinematic I think), defined as, 
		
		nu = 0.00345.


	The boundary conditions(grenseflatebetingelser) are ,
		
		u = g_D  on the Dirichlet part of the boundary, and  
		du/dn - pn = g_N  on the Neuman part of the boundary,  

	gD is a kinematic boundary condition enforced strongly, defined as the exact solution of the Hagen-Poiseuille flow in a 	cylindrical tube, 
	
		gD =  (0, 1/4*nu*dpdz ((x- x0)^2 - R^2)  )
	
	where dpdz is the pressure gradient along the tube, chosen to be 
		dpdz = -1,
		
	gN is a dynamic boundary condition also termed a 'lazy', 'do-nothing' pressure condition which is weakly enforced, here
		gN = 0,  

	where u and p are velocity and pressure respectively - the unknowns to be determined. 



	******************  THE WEAK FORMULATION OF THE EQUAITONS   ****************************************
	

	To arrive at the weak form of the equations, we multiply by a testfunction and integrate over the domain.
	
	First the momentum equations, 

		alpha(u,K)*v*dx - nu*inner(laplace(u), v)*dx + inner(del(u), v)*dx = f*v*dx

	The second term needs to be integrated by parts, 

		- inner(laplace(u), v)*dx = - (- inner(grad(u), grad(v))*dx + inner(du/dn, v)*dS  )  

		- inner(laplace(u), v)*dx = inner(grad(u), grad(v))*dx - inner(grad(u),v)*dS
		 
	where dS is the lineintegral over facets in the interior. The lineintegral over the facets can be 
	divided into two parts, one for the Dirichlet part of the boundary and one for the Neumann part of the boundary 
	
		-inner(grad(u),v)*dS = -(inner(grad(u), v)*dS_D + inner(du/dn,v)*dS_N)
		
	now we can use the boundary conditions and both terms vanish. 
	
	 	


	u = g_D  on the Dirichlet part of the boundary, and  
	du/dn - pn = g_N  on the Neuman part of the boundary,  


	**************** THE DISCRETIZATION OF THE WEAK FORMULATION OF THE EQUATIONS ***********************
	
	
	V_h is a subset of V = span{Lagrange basis funcitons of order 2} 
	Q_h is a subset of Q = span{Lagrange basis functions of order 1}
	W_h is a mixed 	
	
	
	



	****************** THE IMPLEMENTAION PART  **************************



	The input parameters are nu,dpdz, N,x0 ,R  

		nu - a constant dynamic viscosity, type float 
		dpdz - a constan pressure gradient, type float 
		N  - number of cells in the x- and y- direction
		x0 - center point of the geometry around which the solution is symmetric
		R - the radius of from the center point 


	The output is 


	"""

	print 'solving for %s' %N




		
	# *********************** SOLVING THE WEAK FORMULATION ****************************
	


	# PRECONDITION AND SOLVE PARAMETERS 
	# Test for PETSc or Tpetra
	if not has_linear_algebra_backend("PETSc") and not has_linear_algebra_backend("Tpetra"):
	    info("DOLFIN has not been configured with Trilinos or PETSc. Exiting.")
	    exit()

	if not has_krylov_solver_preconditioner("amg"):
	    info("Sorry, this demo is only available when DOLFIN is compiled with AMG "
		 "preconditioner, Hypre or ML.")
	    exit()

	if has_krylov_solver_method("minres"):
	    krylov_method = "minres"
	elif has_krylov_solver_method("tfqmr"):
	    krylov_method = "tfqmr"
	else:
	    info("Default linear algebra backend was not compiled with MINRES or TFQMR "
	         "Krylov subspace method. Terminating.")
	    exit()

	# THE MESH AND FUNCTION SPACES OVER THE MESH 
	mesh = 	UnitSquareMesh(N,N)
	h = 1./N # The element size 
		
	Vspace = VectorFunctionSpace(mesh, 'Lagrange', 2) 
	Pspace = FunctionSpace(mesh, 'Lagrange', 1)       
	Wspace = MixedFunctionSpace([Vspace, Pspace]) 
	Kspace = FunctionSpace(mesh, 'DG', 0)

	# Test and trial spaces over the elements  
	u, p = TrialFunctions(Wspace)
	v, q = TestFunctions(Wspace)
	# Source term
	f = Constant([0.0,0.0])




	# DEFINITON OF THE INTERFACE FUNCTION 
	K_expr = Expression('( x[0] > (x0 - R - DOLFIN_EPS) && x[0] < (x0 + R + DOLFIN_EPS)) ? 0.0 : 1.0', x0 = x0, R = R )
	# INTERPOLATED FOR VISIALIZATION	
	K = interpolate(K_expr, Kspace)	
		

	
	# VARIATIONAL FORMULATION 
 	a = inner(alpha(u, K), v)*dx + nu*inner(grad(u), grad(v))*dx + div(u)*q*dx + div(v)*p*dx
	L = inner(f, v)*dx
	UP = Function(Wspace)




	# BOUNDARY CONDITION
	# Only need bcs on the inlet; surrounding homegenous Dirichlet bcs will result in fenics not finding facets mathing domain for bcs   
	# DEFINING PHYSICAL BOUNDARIES 
	def u_boundaryBottom(x, on_bnd):
		return x[1] < DOLFIN_EPS and x[0] > (x0 - R - DOLFIN_EPS) and x[0] < (x0 + R + DOLFIN_EPS) #and on_bnd

	# ASSIGNING THE PHYSICAL BOUNDARY CONDITIONS  
	# ANALYTICAL SOLUTION IS INLET BCS (grenseflatebetingelse - dynamisk)  
	analyticalVelocity = Expression(["0","0.25*nu*dpdz*((x[0]-0.50)*(x[0]-0.50) - R*R) "], nu=nu, R=R, dpdz=dpdz)
	inlet  = DirichletBC(Wspace.sub(0), analyticalVelocity, u_boundaryBottom)
	bcs = [inlet] 





	# SOLVING THE LINEAR SYSTEM OF EQUATIONS  

	A, b = assemble_system(a, L, bcs)
        solver = LUSolver(A)
	solver.solve(UP.vector(), b) 
	U, P = UP.split()  
        








	# **************  CALCULATING THE ERROR/DISCREPANCY FIELD  *********************
	

	# COMPARISON ON THE FINE MESH; 128X1278 

	u_experimental_onFineMesh = interpolate(U, fine_V)
	u_exact_onFineMesh = interpolate(u_c, fine_V)
	# Make an available vector in the same space to hold the error feld 
	error_onFineMesh = Function(fine_V)	
	# The error field is the discrepancy between u_exact - u 
	error_onFineMesh.vector()[:] = u_exact_onFineMesh.vector()[:] - u_experimental_onFineMesh.vector()[:] 



	# COMPARISON ON MESH; NXN FOR N IN [8,16,32,64,128] 
	
	u_experimental_NxNmesh = interpolate(U, Vspace )	
	#u_exact_NxNmesh = interpolate(u_c, Vspace )
	# Again make a vector available in the same space 
	#error_NxNmesh = Function(Vspace)	
	#error_NxNmesh.vector()[:] = u_exact_NxNmesh.vector()[:] - u_experimental_NxNmesh.vector()[:]  



	# FASTER VECTORIZED VERSION; ALSO AXPY WORK IN PARALLEL; NOT THE BEST FOR PLOTTING  	
	#error_3  = interpolate(U, fine_V)
	#error_3.vector().axpy(-1, interpolate(u_c, fine_V).vector())
        #plot(error, title=str(N))
	








	# ******************* PLOTTING ****************************** 

	# Plotting for fine mesh  
	File('K{}.xdmf'.format(N)) << K
	#File("u_experimental_onFineMesh{}.xdmf".format(N)) <<  u_experimental_onFineMesh 	
	File('u_exact_onFineMesh{}.xdmf'.format(N)) << u_exact_onFineMesh
	File('u_error_onFineMesh{}.xdmf'.format(N)) << error_onFineMesh
		
	
	#Plotting for nxn mesh
	File("u_experimental_NxNmesh{}.xdmf".format(N)) <<  u_experimental_NxNmesh 	
	#File('u_exact_NxNmesh{}.xdmf'.format(N)) << u_exact_NxNmesh
	#File('u_discrepancy_NxNmesh{}.xdmf'.format(N)) << error_NxNmesh
	
 





	# ***************** CALCULATING THE ERRORNORM *****************
	
	# 1) EXPLICITLY  

	# might have to do something with Martins analytical formula
	# can use 6 and 8 and se if it makes a difference? 
	# then in the error I think  i might see a difference in the decimals 	
	q = 2 
	degree = 2*(q + 2)	
	norm_of_error = sqrt(assemble((u_c - U)**2*dx(degree =degree)))
	# numerisk tilnaelse til integralet 
	# polynom med order degree 
	# specify the accuracy of the integration , the higher the more accurate 




	# 2) USING ERRORNORM METHOD FROM THE FENICS LIBRARY 
	u_error = errornorm(u_c, U, norm_type='L2', degree_rise=4, mesh=mesh)
	#l2 is the vector norm, while L2 is the norm for the function space
        


	return norm_of_error, h, u_experimental_NxNmesh, u_exact_onFineMesh, error_onFineMesh 

	








if __name__ == '__main__':


	# DEFINING THE MESH PARAMETERS AND FLUID FLOW PARAMETERS 
	center = 0.50      # CENTER OF THE TUBE 
	side_width = 0.1   # RADIUS OUT FROM THE CENTER 
	nu = 0.00345       # DYNAMIC FLUID VISCOSITY, MATERIAL PROPERTY  
	dpdz = -2.1        # PRESSURE GRADIENT OF THE FLUID FLOW 
	

	# SENDING PARAMETERS TO C++ FUNCTION 
	u_c = Expression(u_analytical) 
	u_c.x0 = float(center);
	u_c.R  = float(side_width);
	u_c.nu = float(nu);
	u_c.dpdz = float(dpdz);



	# FINE MESH FOR ANALYTICAL SOLUTION AND COMPARISON BETWEEN ANALYTICAL AND EXPERIMENTAL SOLUTION  	
	fine_mesh = UnitSquareMesh(64,64) #128,128)
	fine_V = VectorFunctionSpace(fine_mesh, 'Lagrange', 2)
	fine_P = FunctionSpace(fine_mesh, 'Lagrange',1)		




	
	# CALCULATING ERROR VALUES AS A FUNCTION OF MESH SIZE 
	error_values = []
	h_values = []
	for n in [8,16,32,64, 128]: #256 ,512]: 
		error_i, h_i, U, U_exact, U_discrepancy =   experimentalSolver(nu=nu, dpdz=dpdz, N=n,x0=center ,R=side_width) 		
		error_values.append(error_i)
 		h_values.append(h_i)

	print 'Mesh size, h=1/N          Error  '
	for i, j in zip(h_values, error_values):
		print '%5.3f        %10.5e' %(i,j) 





	# PLOTTING A LOGARITMIC PLOT OF THE ERROR AS A FUNCTION OF THE MESH SIZE 
	import pylab as pylab 
	enable_plot = True
	if enable_plot:
		pylab.loglog(h_values, error_values, 'g-o' ) 
		#pylab.plot(h_values, error_values, 'b-o')	
		pylab.xlabel('h')
		pylab.ylabel('error')
		pylab.grid(True)
		pylab.show()		
		
 		
	
	

	

	




