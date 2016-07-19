
#mesh size
h = [0.062, 0.031, 0.016, 0.008]


#error values for different c 
c01   = [1.0851E-1, 1.0738E-1, 1.0712E-1, 1.0773E-1 ]
c1    = [1.0423E-1, 1.04897E-1, 1.00458E-1, 1.0517E-1]
c10   = [1.01387E-1, 1.02789E-1, 1.02375E-1, 1.02934E-1]
c100  = [9.6479E-2, 9.7527E-2, 9.6676E-2, 9.7030E-2]
c1000 = [8.388E-2, 7.910E-2, 7.560E-2, 7.442E-2]
c1E4  = [0.07149, 0.04664, 0.0379928, 0.034147]
c1E5  = [1.226E-1, 3.352E-2, 2.3182E-2, 2.468E-2]



# PLOTTING A LOGARITMIC PLOT OF THE ERROR AS A FUNCTION OF THE MESH SIZE 
import pylab as pylab 
enable_plot = True #False
if enable_plot:
		pylab.loglog(h,c01,'b',h,c1,'g',h,c10,'r',h, c100,'c',h,c1000,'m',h,c1E4,'y',h,c1E5,'k') 
		#pylab.plot(h_values, error_values, 'b-o')
		#pylab.axis([-5, 5, -5, 5 ])	
		pylab.xlabel('mesh size, h')
		pylab.ylabel('H1 Error norm ')
		pylab.title('Velocity Converge for C - values')
		pylab.legend(['C=0.1','C=1', 'C=10','C=100', 'C=1000', 'C=1E4', 'C=1E5'])
		#pylab.grid(True)
		pylab.savefig('loglog_velocity')			
		pylab.show()	




"""
import matplotlib.pyplot as plt
import numpy as np

x = np.linspace(0, 2, 100)

plt.plot(x, x, label='linear')
plt.plot(x, x**2, label='quadratic')
plt.plot(x, x**3, label='cubic')

plt.xlabel('x label')
plt.ylabel('y label')

plt.title("Simple Plot")

plt.legend()

plt.show()
""" 
	
