
def list():
	alist = [1.2, 1.6, 1.5, -0.1]
	blist = []

	for a in range(len(alist)):
		print alist[a]
		if alist[a] > 1.2:
			print alist[a]
			alist[a] = 'hei'
		elif alist[a] <0:
			alist[a] = 0			

		else:
			alist[a] = alist[a]	
	
	print alist 
	

list()

