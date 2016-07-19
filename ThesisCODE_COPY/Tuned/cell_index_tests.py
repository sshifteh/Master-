from dolfin import * 
import numpy as numpy 


def cell_index(dof, V):
   # Returns the index of a cell associated with the given dof index.
   # One dof can be associated with multiple cells, this functions simply returns the first found.
   fg_til_funksjonsrom = V.dofmap()	
   #print(fg_til_funksjonsrom)
   mesh_til_funksjonsrom = V.mesh()
 

   for cell in cells(mesh_til_funksjonsrom):                             # for hver celle i alle cellene i meshet 
 	if dof in fg_til_funksjonsrom.cell_dofs(cell.index()):  # hvis gitt dof er i dofene generert av funksjonsrommet (globale), med celler indeksert fra meshet 
		#print cell.index()                                      # returner indeksen til cellen med den dofen ( noden)
		return cell.index()








def neighbouring_cells(cell_index, mesh):
    # Returns a list of cells neighbouring cell with given index.
    # The given index is also included
    mesh.init()                      # husker ikke hva gjoer, se gamle notater 
    nabocelle_indekser = set()       # lager et tomt set som skal romme indeksene til alle nabocellene til input cellen 
    				     	
    # loop over facets to find neighbours
    #print ' \n Cell(mesh, cell_index) \n ',Cell(mesh, cell_index)	
    for facet_id in Cell(mesh, cell_index).entities(mesh.topology().dim()-1):    # for en sidekant i utvalgt celles sidekant array. arrayen inneholder nummber paa sidekantene fx [0,1,2]   
        facet = Facet(mesh, facet_id)                                            # facet settes til celle indeks for celler med sidekant som i arrayet    
        # facet is connected to two cells                                        # facet gjoer om tallene som er i denne facet arrayen til facet objekter       
        nabocelle_indekser.update(facet.entities(mesh.topology().dim()))         # for aa finne facet objektenes entities 2 som er celler. dvs koblin til celler 

    return nabocelle_indekser # disse cellene indeksene puttes i nabocelle_indekser. og det gjenstaar aa finne dof for disse 




def cell_dofs_by_entity_dim(cell_index, i, space):
    """
    cell_index : (nabo)celle indeks	
    i = topological dimension 
    V = Funksjonsrom 	
    """


    # Returns a dict of dofs on given cell and dim
    mesh = space.mesh()     # meshet 
    dofmap = space.dofmap() # dofmapet til meshet. ulike funksjonsrom kan vaere over samme mesh

    cell_dofs = dofmap.cell_dofs(cell_index)   # dofene til en celle er funnet ved aa ta dofmappet til meshet of mappe lokale til globale dofs for gitt celle = cell_index  
    dofs_i = set()                             # igjen lag et tomt set som skal romme dofene for den gitte cellen 
    for j in xrange(len(Cell(mesh, cell_index).entities(i))):  # cell().entities(0) returnerer en liste av dofene til den cellen fx [0,1,4]. saa vi loper gjennom disse  # i = 0 kanter
        dofs_i.update(cell_dofs[dofmap.tabulate_entity_dofs(i,j)])  # tror denne finner indeksen til dofen ? rart. er ikke dof nr. indeksen paa noden allerede? spoer om denne linjen.   

    return dofs_i



if __name__ =='__main__' :

    	
	mesh = UnitSquareMesh(2,2)
	plot(mesh, interactive = True)
	CG1 = FunctionSpace(mesh, "CG", 1)
	DG0 = FunctionSpace(mesh, "DG", 0)
	DG1 = FunctionSpace(mesh, 'DG', 1)
	CG3 = FunctionSpace(mesh, 'CG', 3)
	
	
	
	# 1: Finne indeksen pa en celle vi har dofen paa: 
	
	random_dof = numpy.random.randint(DG1.dim())  
	#print ' \n random valgt dof no: \n ', random_dof
	
	# get cell index and neighbouring cells
	dof = 1
	print ' \n dof:', dof
	cell_id = cell_index(dof, DG1)   # saa vi gir den det tallet og funksjonsrommet 
	print 'celle id: ' ,cell_id 
	
	
	
	
	# 2 : finne indeksen paa nabocellene via sidekantene til cellen vi har indeks paa :
	
	print 'cell(mesh, cell_id).entities(0):', Cell(mesh, cell_id).entities(0)
	print 'cell().entities(1):', Cell(mesh, cell_id).entities(1)
	print 'cell().entities(2):', Cell(mesh, cell_id).entities(2)
	
	facet_no = 6
	print 'Facet(mesh, facet_no = 6):',Facet(mesh, facet_no)
	# denne gir sidekantene til celle no 5 i meshet 
	print 'Facet(mesh, facet_no = 6).entities(2) dvs celler:',Facet(mesh, facet_no).entities(2)
	
	
	nabocelle_indeks_liste = neighbouring_cells(cell_id, mesh)
	print 'nabocelle_id, dvs index til alle naboceller: ', nabocelle_indeks_liste 	
	# nabocelle indeks er en liste med indeksen paa naboceller til input celle
	# for hvert celle indeks i den listen, printer vi ut dofene
	
	
	
	# STEP 3: Finne dofene til nabo cellene vi har indeksen paa :
	
	#print 'cell_dofs_by_entitity_dim(nabocelle_indeks, 0, DG1):', cell_dofs_by_entitity_dim(nabocelle_indeks, 0, DG1)
	
	for celle_indeks in nabocelle_indeks_liste:
		print celle_indeks
		print 'dof til nabocelle : ', celle_indeks, 'er', cell_dofs_by_entity_dim(celle_indeks, 2, DG1)
		
	
	
	""" Example for cg3 
	for celle_indeks in nabocelle_indeks_liste:
		print celle_indeks
		print 'dof til nabocelle : ', celle_indeks, 'er', cell_dofs_by_entitity_dim(celle_indeks, 2, CG3)
	""" 
	
	
	
	




print ''

