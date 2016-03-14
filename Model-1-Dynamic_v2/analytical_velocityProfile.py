from dolfin import * 



# Del 1 ------------------------------------------------------------------------------ 
""" 
mesh = UnitSquareMesh(2,2)
class Omega0(SubDomain):
	def inside(self, x, on_bnd):
		return True if x[1] <= 0.5 else False 
	
class Omega1(SubDomain):
	def inside(self, x, on_bnd):
		return True if x[1] >= 0.5 else False 

subdomains = MeshFunction('size_t', mesh, 2) # 'size_t' = 'uint'

#marker hvert subdomenene med 0 og 1 
subdomain0 = Omega0()
subdomain0.mark(subdomains, 0)

subdomain1 = Omega1()
subdomain1.mark(subdomains, 1)

#print type(subdomains) # MeshFunction
#plot(subdomains, interactive = True)


V0 = FunctionSpace(mesh, 'DG', 0)
k = Function(V0)


# verdien av k i de to subdomenene
k_values = [1.5, 50]
#print len(subdomains.array()) # 8 celler vil det si 
# itererer over cellene nummerene 0,1,2,3,4,5,6,7 
for cell_no in range(len(subdomains.array())):
	#print cell_no	

	#definerer en variable med navn subdomain_no
	#print subdomains.array()[cell_no]
	cell_value = subdomains.array()[cell_no]  # itererer over alle cellene 
	#print k.vector()[cell_no]
	#print k_values[cell_value]

	k.vector()[cell_no] = k_values[cell_value] # k tar verdien til k_values i den cellen 


#plot(k, interactive = True)

""" 


# Del II---------------------------------------------------------------------------------------------------------- 


mesh2 = UnitSquareMesh(32*4, 32*4) 



class MyExpression(Expression):
	def __init__(self, x0 = 0.5, R = 0.05, nu = 0.00345, dpdz = 2.1):
		self.x0 = x0
		self.R = R
		self.nu = nu
		self.dpdz = dpdz

	def eval(self, values, x):
		#print type(values), values.flags, 'hei paa deg'		 
		# TODO hvorfor er den satt til False i biblioteket? en god grunn? en feil? 		
		# HACK: 		
		#values.flags.writeable = True
		# x er immutable. fordi funksjonen i biblioteket tar en vektor med verdier og posisjonsvektoren
		# vektoren med verdier skal bli oppdatert mhp posisjon 
		# v(r) og r er alltid (x,y,z)
 
		if  self.x0 - self.R <= x[0] <= self.x0 + self.R:
			# inside pipe
			values[0] = 0.0			
			values[1] = 0.25*self.nu*self.dpdz*((x[0]-self.x0)*(x[0]-self.x0) - self.R*self.R)

		else:
			# outside pipe 
			values[0] = 0.0
			values[1] = 0.0

	# Spesify the value shape of the expression 
	def value_shape(self):
		return (2,)



# Naar den skrives i C ++ . Siden dolfin er i C ++. Saa er veien mellom de to kodene kortere muligens. 
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
			values[1] = velocity ; 
		}
		else 
		{
			values[0] = 0.0;
			values[1] = 0.0; 
		}
	}
};''' 


nu = Constant(0.00345) # Pascal x second 
dpdz = Constant(-2.1)
x0 = Constant(0.50)
R = Constant(0.25)	

g_C = Expression(u_analytical) 
g_C.x0 = float(x0(0));
g_C.R  = float(R(0));
g_C.nu = float(nu(0));
g_C.dpdz = float(dpdz(0));


"""
# Python code 
# Velocity parameters 
nu = Constant(0.00345) # Pascal x second 
dpdz = Constant(-2.1)
V = VectorFunctionSpace(mesh2, 'Lagrange', 2)
#P = FunctionSpace(mesh2, 'Lagrange', 1)
#W = V*P
g = MyExpression(R = 0.25, dpdz = dpdz) # Myexpression er en funksjon som gitt et punkt gir et verdi 
# g har gitt en graf, den kan vaere diskont, eller kvadratisk eller hva som helst 

# men naar vi interpolerer det far vi en funksjon i det rommet vi har bestemt den skal vaere , her V. 
f = interpolate(g, V)  # linear kombinasjoner av basisen for V rommet. 
# V inneholder kvadratiske basisfunksjoner saa f blir nlyaktig lik g. 
# da blir den reprodusert eksakt 
plot(g, mesh = mesh2)


#v
#diff = errornorm(g,v)
# mulig errornorm interpolerer til et hoyere dim rom for den regner ut erroren i en eller annen norm 


plot(f)
interactive()
""" 



