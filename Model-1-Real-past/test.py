from numpy import where, arange 



a = arange(10)

a = [1,2,3,4,5,6,7]

b = where(a < 5)[0]
print b

#growth = where((abs(wss_) > thresh_H) & (bdry_ > DOLFIN_EPS))[0]
