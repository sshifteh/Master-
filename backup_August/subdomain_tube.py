from dolfin import * 

n = 100 # number of cells  
w = 0.2 # widt 

mesh = UnitSquareMesh(n,n)                      # unit mesh of n lengh and n widt                         
subdomain = CellFunction('size_t', mesh, 0)    # Hva skjer her, computer structure

class Tube(SubDomain):                          # Tube arver Subdomain, som har en inside metode som definerer en subklasset av unitsquare 
    def inside(self, x, on_boundary):
        return  0.5-w-DOLFIN_EPS  < x[0] < 0.5 + w + DOLFIN_EPS 
        
Tube().mark(subdomain, 0)                      # mark funksjonen gaar over alle cellene og sjekker om inside metoden for den cellen returnerer true
plot(subdomain)
interactive(True)



V0 = FunctionSpace(mesh, "DG", 1)                 # Funksjonsrom over meshet med diskontinuerlig galerkin av grad 1 
f  = Expression('(x[0] > 0.3 && x[0] < 0.7 && x[1]< 0.5 ) ? 3.0 : 10.0')


k  = interpolate(f, V0)

plot(k, interactive = True, title = "Subdomains")     # Viser homogent omraade, k_values er ikke fordelt ennaa                                               

