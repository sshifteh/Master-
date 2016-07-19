# -*- coding: utf-8 -*-
"""
Created on Wed Jun 24 22:30:02 2015

@author: shiftehs

Plan:
First, create the geometry 
Second, create the matrix and create a discrete function over the geometry
Third, code the variational formulation, assemble and solve
obs the integrals ! 

"""

from dolfin import * 
N = 10 
mesh = UnitSquareMesh(N,N)
plot(mesh, interactive = True)






# Creating the matrix K 
# Create a discrete function over the mesh for the matrix 
# it symmetric use A = U+D+L , save only upper triangular part 
a00 = MeshFunction('double', mesh, 2)  # double is the object type
a01 = MeshFunction('double', mesh, 2)  # mesh the MeshFunction is defined on
a11 = MeshFunction('double', mesh, 2)  # 2 stands for 2D

# Test for the domain that we are in
a = 0.5
for cell in cells(mesh):
    if cell.midpoint().x() < a:        
        a00[cell] = 1.0
        a01[cell] = 0.0
        a11[cell] = 1.0
    else: 
        a00[cell] = 0.0
        a01[cell] = 0.0
        a11[cell] = 0.0

# Store the data. Dont understand this part really. Why is it necessary? 
# Store to file
mesh_file = File("mesh.xml.gz")
a00_file = File("a00.xml.gz")
a01_file = File("a01.xml.gz")
a11_file = File("a11.xml.gz")

mesh_file << mesh
a00_file << a00
a01_file << a01
a11_file << a11


# Plot mesh functions
plot(a00, title="a00", interactive = True)
plot(a01, title="a01", interactive = True)
plot(a11, title="a11", interactive = True)

#interactive

