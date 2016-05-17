


from dolfin import *

N = 32
x0 = 0.5
R = 0.1 # ops den er assymetrisk 11 celer paa lhs og 12 paa rhs 
viz = True
mesh = UnitSquareMesh(N,N)
dg0 = FunctionSpace(mesh, 'DG',0)
cg1 = FunctionSpace(mesh, 'CG',1)
dg1 = FunctionSpace(mesh, 'DG',1)


K = Expression('( x[0] > (x0 - R ) && x[0] < (x0 + R )) ? 0.0 : 1.0', x0 = x0, R = R )
K = interpolate(K, cg1) # alternativt dg0 
plot(K, interactive = viz, title = 'K' )			






# Alternativ 1 	
interface = grad(K)**2
interface_ = project(interface, dg1)
#Det er 3 alternativer for rommet til interface, eller vent kanskje konverteres til et annet rom for plotting 
# dg0 funker 
# cg1 funker 
# dg1 funker 

print type(interface_) 
plot(interface_, interactive = viz, title = 'interface ')	






	
# Alternativ 2 
#n = FacetNormal(mesh) 
#boundary = grad(K)*n  
# da matte jeg foerst ha markert cellene paa interfacet med Cell Funcition og mark tror jeg. Og da maatte jeg gaa gjennom meshet uansett. 

	

