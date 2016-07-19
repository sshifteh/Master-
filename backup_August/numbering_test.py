from dolfin import * 
import numpy 

n = 2
mesh = UnitSquareMesh(n,n)
f = FacetFunction("size_t", mesh)  # a vector of size n
f.set_values(numpy.arange(f.size(), dtype = "uintp")) # numbering all the edges in 2d 
print f.array()                    # 
plot(f, domain = mesh, interactive = True)
