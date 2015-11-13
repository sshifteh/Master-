

from dolfin import * 

n = 100 
mesh = UnitSquareMesh(n,n)                        # Antallet elementer er 8

V0 = FunctionSpace(mesh, "DG", 1)                 # Funksjonsrom over meshet med diskontinuerlig galerkin av grad 1 
f1 = Expression('(x[1] > 0.5 && x[1] < 1-x[0]  ) ? 3.0 : 10.0')
#f2 = 
f3  = Expression('(x[0] > 0.3 && x[0] < 0.7 && x[1]< 0.5 ) ? 3.0 : 10.0')

k  = interpolate(f1, V0)
k2 = interpolate(f3, V0)

plot(k,k2, interactive = True, title = "Subdomains")     # Viser homogent omraade, k_values er ikke fordelt ennaa                                               

