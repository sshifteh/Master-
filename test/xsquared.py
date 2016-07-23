

def f2(x):
	return x*x


def f3(x):
	x*x*x

def f4(x):
	x*x*x*x



if name == "__main__":
	from numpy import *
	x = linspace(1, 10,1000)
	xsquared = f2(x)
	xcubed = f3(x)
	xfourth = f4(x)
	
