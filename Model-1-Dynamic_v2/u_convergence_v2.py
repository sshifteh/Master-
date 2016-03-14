from alpha import alpha
from dolfin import *

	
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
	// Making the array 2 dimensional 
        size_t value_rank() const{
            return 1;
        }

        size_t value_dimension(size_t i) const{
            return 2;
        }

};''' 


# Experimental solver and error 
def main(nu,dpdz, N,x0 ,R ):  
	"""
	Solves the Stokes equation with an additional inverse permeability term. 
	Calculates the error between analytical and experimental velocity.	
	"""
	
	mesh = 	UnitSquareMesh(N,N)
	Vspace = VectorFunctionSpace(mesh, 'Lagrange', 1)
	Pspace = FunctionSpace(mesh, 'Lagrange', 1)
	Wspace = MixedFunctionSpace([Vspace, Pspace]) 
	Kspace = FunctionSpace(mesh, 'DG', 0)

	# The vessel geometry and interpolation into the Kspace 
	K_expr = Expression('( x[0] > (x0-R) - DOLFIN_EPS && x[0] < (x0+R) + DOLFIN_EPS) ? 0.0 : 1.0', x0 = x0, R = R )
	
	
	K = interpolate(K_expr, Kspace)	

	# Test and trial spaces over the Taylor-Hood elements  
	u, p = TrialFunctions(Wspace)
	v, q = TestFunctions(Wspace)

	# Source term
	f = Constant([0.0,0.0])

	# Variational form 
 	a = inner(alpha(u, K), v)*dx + nu*inner(grad(u), grad(v))*dx + div(u)*q*dx + div(v)*p*dx
	L = inner(f, v)*dx
	UP = Function(Wspace)

	# Parabolic inlet velocty Dirichlet condition 
        speed = 10
	analyticalVelocity = Expression(["0"," 0.25*nu*dpdz*((x[0]-0.50)*(x[0]-0.50) - R*R) "], nu=nu, R=R, dpdz=dpdz)

	# Defining the boundaries of the vessel 
	def u_boundaryBottom(x, on_bnd):
		return x[1] < DOLFIN_EPS and x[0] > (x0-R) - DOLFIN_EPS and x[0] < (x0+R) + DOLFIN_EPS and on_bnd

	def restBottom(x, on_bnd):

		return x[1] < DOLFIN_EPS and  not (x[0]> (x0-R) -DOLFIN_EPS) and not (x[0] < (x0+R) + DOLFIN_EPS) and on_bnd

    	def leftWall(x, on_boundary):
        	return near(x[0], 0.0)

	def rightWall(x, on_boundary):
        	return near(x[0], 1.0)


	def top_rest(x, on_bnd):
		return x[1] < 1 - DOLFIN_EPS and not (x[0] > (x0-R)) and not (x[0] < (x0+R)) and on_bnd

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
        
	# Calculating the difference between the analytical and experimental solution on a fine mesh  
        error  = interpolate(U, fine_V)
	error.vector().axpy(-1, interpolate(u_c, fine_V).vector())
        plot(error, title=str(N))
	#from IPython import embed; embed()

	# Plotting the experimental solution on a superfine mesh with 128 cells 
	plot(interpolate(U, fine_V),mesh, title = 'Experimental Velocity Field, N')
	# Plotting the analytical solution on a fine mesh 
	plot(interpolate(u_c, fine_V),mesh, title = 'Analytical Velocity Field, N')
	interactive()
	# Calculating the errornorm
	u_error = errornorm(u_c, U, norm_type='L2', degree_rise=4, mesh=mesh)
	#l2 is the vector norm, while L2 is the norm for the function space
        
	return u_error, mesh.hmin()




if __name__ == '__main__':

	center = 0.50
	side_width = 0.1
	nu = 0.00345
	dpdz = -2.1
	
	u_c = Expression(u_analytical) 
	u_c.x0 = float(center);
	u_c.R  = float(side_width);
	u_c.nu = float(nu);
	u_c.dpdz = float(dpdz);


	# Coarse or varying mesh 	
	#V = VectorFunctionSpace(mesh, 'Lagrange', 2)
	#g = MyExpression(R = 0.25) 

	# Create a fine mesh for comparing analytical and numerical solution	
	fine_mesh = UnitSquareMesh(128,128)
	fine_V = VectorFunctionSpace(fine_mesh, 'Lagrange', 1)
	


	E = []; h = []; 
	for n in [8,16,32,64]:
		ei, hi = main(nu = nu,dpdz = dpdz, N = n,x0 = center ,R =side_width)  
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



Simon: 
	poeng1:tilnaermet en kvadratisk funksjon med kvadratiske funksjoner, fikk dermed noeyaktig representasjon av funksjonen
	poeng2: geometrien paa avstand 0.25 passet akkurat til at vi oekte meshet med det dobbelte fra et partall 16,32 osv
	da vi endret paa dette saa vi med en gang feil, og svakheter i koden, som naa er rettet 
	
	P1-P1 elementer gir 2 grads congervens rate 
	mens Taylor-Hood elementer ga sub convergens, foerste gangs konvergens rate men det var paa grunn at eksakt loesning var samme som eksperimentell pga kvadratisk tilnaerming og geometri. 

		
	 	

"""

