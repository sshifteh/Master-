
from numpy import * 
from matplotlib.pyplot import * 

def f(x):
    return  -x*x

start = 0
stop = 1
steps = 10 
x = linspace(start,stop,steps)
y = f(x)

print ' x       y    '
for i,j in zip(x,y):
    print '%2.2f %5.2f ' %(i,j) 

plot(x,y)
show()
