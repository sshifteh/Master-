from dolfin import * 
import sys 

element = FiniteElement('BDM', tetrahedron, 3)   # Brezzi-Douglas-Marini
#element = FiniteElement('TH', tetrahedron, 3)   #
plot(element)
sys.exit()


