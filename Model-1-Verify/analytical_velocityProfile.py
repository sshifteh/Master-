from dolfin import * 



# Del 1 ------------------------------------------------------------------------------ 

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




# Del II---------------------------------------------------------------------------------------------------------- 

mesh2 = UnitSquareMesh(10,10) 
class MyExpression(Expression):
	def eval(self, x, values):
		#print type(values), values.flags, 'hei paa deg'		 
		# TODO hvorfor er den satt til False i biblioteket? en god grunn? en feil? 		
		# HACK: 		
		values.flags.writeable = True 
		if  0.45 <= x[0] <= 0.55:
			# inside pipe
			values[0] = 0.25*nu*dpdz*((x[0]-0.50)*(x[0]-0.50) - R*R)			
			values[1] = 0

		else:
			# outside pipe 
			values[0] = 0
			values[1] = 0

	# Spesify the value shape of the expression 
	def value_shape(self):
		return (2,)


# Velocity parameters 
nu = Constant(0.00345) # Pascal x second 
dpdz = Constant(2.1)
R = Constant(0.05)

V = VectorFunctionSpace(mesh2, 'Lagrange', 2)
#P = FunctionSpace(mesh2, 'Lagrange', 1)
#W = V*P
	
g = MyExpression()
f = interpolate(g, V)
plot(f)
interactive()

