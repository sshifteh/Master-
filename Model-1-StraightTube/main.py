from alpha import alpha 
from stokes_solver import stokes_solver 
from domain1 import domain
#from domain2 import domain2 
from dolfin import * 

def main():
	""" The old main  """"


	# Parameters
	n = 40
	w   = 0.3
	show = True 
	
	# Mesh and Functionspaces
	mesh = UnitSquareMesh(n,n)
	Vspace = VectorFunctionSpace(mesh, 'Lagrange', 2) 
	Pspace = FunctionSpace(mesh, 'Lagrange', 1)
	Wspace = MixedFunctionSpace([Vspace, Pspace]) 
		
	ShearStressSpace=FunctionSpace(mesh, 'DG', 1)
	TensorSpace = TensorFunctionSpace(mesh, "DG", 1) # Space for the shearstess,could be named ShearStressSpace			
	DPspace = FunctionSpace(mesh, 'DG', 0)	         # Space for the gradient of the wss 
	
	# Expression for the initial T shape 
	K_1 = Expression('(x[0] > w && x[0] < 1- w && x[1] < 0.5 || x[1] >= 0.5 && x[1] <= 0.7 )? 0.0:1.0', w=w) #K in math.formulation	
	
	# Function calls
	U, P, K_Func = stokes_solver(w=w, mesh=mesh, Vspace=Vspace,Pspace=Pspace, Wspace=Wspace, TensorSpace=TensorSpace, 			DPspace=DPspace, ShearStressSpace=ShearStressSpace , K_array=K_1, n=n)

	#K_Func = domain(mesh = mesh, TensorSpace=TensorSpace, DPspace=DPspace, ShearStressSpace=ShearStressSpace, Pspace=Pspace, U=U, K_Func = K_Func, n=n)

	stokes_solver(w,mesh,Vspace,Pspace,Wspace,TensorSpace,DPspace,ShearStressSpace, K_array=K_Func, n=n)
 	
	#domain2(Vspace=Vspace, mesh=mesh, TensorSpace=TensorSpace, DPspace=DPspace,ShearStressSpace=ShearStressSpace, Pspace=Pspace, U=U,K_Func=K_Func)
	
	

if __name__ == '__main__':
	# Import statements	
	import time 
	import numpy as np
	import matplotlib.pyplot as plt 
	from dolfin import * 	
	main()

	
