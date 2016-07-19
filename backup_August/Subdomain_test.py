from dolfin import * 

n = 100 # number of cells  

mesh = UnitSquareMesh(n,n)                        # Antallet elementer er 8

V0 = FunctionSpace(mesh, "DG", 1)                 # Funksjonsrom over meshet med diskontinuerlig galerkin av grad 1 
f  = Expression('(x[0] > 0.3 && x[0] < 0.7 ) ? 3.0 : 10.0')


k  = interpolate(f, V0)

plot(k, interactive = True, title = "Subdomains")     # Viser homogent omraade, k_values er ikke fordelt ennaa                                               

