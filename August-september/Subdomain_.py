from dolfin import * 

n = 2 
mesh = UnitSquareMesh(n,n)                      # Antallet elementer er 8

class Omega0(SubDomain):
    def inside(self, x, on_boundary):            # Lager grenser for subdomene1
        return True if x[1] <= 0.5 else False 
class Omega1(SubDomain):
    def inside(self, x, on_boundary):            # Lager grenser for subdomene2
        return True if x[1] >= 0.5 else False 
        
subdomains = CellFunction("size_t", mesh)        # Lager et mesh over subdomenene
print 'Length of the subdomains.array is = ',len(subdomains.array())     # Denne har 8 celler 
print 'The subdomains array ois : ',subdomains.array()

subdomain0 = Omega0()                            # Subdomains er en liste av to subdomener indeksert med 0 og 1, i.e. [subdomain0, subdomain1] 
subdomain0.mark(subdomains, 0)                              

subdomain1 = Omega1()
subdomain1.mark(subdomains, 1)

print 'The subdomains array ois : ',subdomains.array()

V0 = FunctionSpace(mesh, "DG", 0)                 # Funksjonsrom over meshet med diskontinuerlig galerkin av grad 1 
k = Function(V0)                                 # Konstantfunksjon over hvert subdomene 
k_values = [1.5, 100]                            # Liste over material verdiene for hvert subdomene 

plot(k, interactive = True, title = "Subdomains")     # Viser homogent omraade, k_values er ikke fordelt ennaa                                               
print 'k.vector().array() : ',k.vector().array()      # k.vector.array is a zero array to be filled with 1.5 or 100,depending on which subdomain it belongs to
print 'len(k.vector().array())', len(k.vector().array()) # 24 med Raviart elementer 


# Cells have to be distributed to one of the two subdomains
for i in range(len(subdomains.array())):
    subdomain_no = subdomains.array()[i]
    k.vector()[i] = k_values[subdomain_no]

print 'k.vector().array() : ',k.vector().array()
plot(k, interactive = True, title = "Subdomains")