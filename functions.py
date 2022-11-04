import math
import numpy as np
import networkx as nx


dots = "-----------------------------------\n"

def algebraicConnectivity(i):
	""" The function that calculates the spectral radius 
	of the bottleneck matrix of the path with i vertices"""
	return 2 * (1 - math.cos(math.pi / i))	

def setInitialValue(ntree):
	"""Order parity is sufficient to determine the i of the first vertex q + i to be analyzed"""
	if(ntree == "even"):
		return 1
	else:
		return 2		

def setLimit(n, ktree):
	"""The i of the last vertex q+i to be analyzed is defined here"""
	if(n == "odd"):
		return ktree // 2
	else:
		if(ktree % 2 == 1):
			return ktree // 2
		else:
			return ktree // 2 - 1
		
	
def calcPolSeq(l, father = 0, grandfather = 0):
	"""Here the l-th term b(l) as defined in the article is calculated."""
	
	if (l == 1):
		return np.array([1])
	elif(l == 2):
		return np.array([1, 0])
	else:
		p0 = np.array([0])
		p1 = np.array([0, 0])

		# Simulates the multiplication of father by s in the calculation of b(l)
		p = np.concatenate((father, p0), axis = None) 
		
		p += np.concatenate((p0, father), axis = None)
		p -= np.concatenate((p1, grandfather), axis = None)
	
		return p

def increaseDimension(p):
	return np.concatenate((np.array([0]), p), axis = None)
		  
def calcPol(order, posicion, ktree):
	"""The function calculates the polynomial whose roots in (0,1) are 1-a(T)"""
	
	# n = 2 * q + b, b depends on the order parity of the broom tree
	b = 0
	if(order != "even"):
		b = 1
	
	# Variable that will be assigned the first term of the sequence b to the term ktree + 1 + b - 2 * posicion 
	aux_poly = 0
	# Variable that will receive the sum of all values assigned to aux_poly
	sum_polynomials = 0
	# In the calculation of the j-th term of the sequence b, father is b(j-1) and grandfather is b(j-2)
	father = 0
	grandfather = 0

	
	
	for l in range(1, ktree + 2 + b - 2 * posicion):
		aux_poly = calcPolSeq(l, father, grandfather)
		
		
		# This first part updates sum_polynomials
		if(l == 1):
			sum_polynomials = aux_poly
		else:
			sum_polynomials = increaseDimension(sum_polynomials) + aux_poly

		# This second part updates the father and grandfather
		#  variables for the aux_poly calculation in the next iteration	
		if(l == 2):
			grandfather = np.array([1])
			father = np.array([1, 0])
		elif(l >= 3):
			grandfather = father
			father = aux_poly	
	
	# Here the first argument is irrelevant, since the calculation 
	# will depend only on the variables father and grandfather
	aux_poly = calcPolSeq(3, father, grandfather)
	
	
	return - ktree * aux_poly +  np.concatenate((sum_polynomials, np.array([0])), axis = None)	
			
def calcAlgebraicConnectivity(p):
	"""This function converts the valid roots of p from 1-a(T) to a(T)"""
	l = list()
	
	for elem in np.roots(p):
		if(elem > 0 and elem < 1):
			l.append(1 - elem)
	return l	

def calcDif(x, i):
	return abs(x - algebraicConnectivity(i))

def createBroom(ntree, ktree):
	""" This function creates a broom tree, whose number of leaves hanging 
	from the star is ktree and the number of vertices is ntree"""
	   
	G = nx.Graph()
	star = nx.star_graph(ktree)
	path = nx.path_graph(ntree - ktree) 
	G = nx.disjoint_union(star, path)  
	G = nx.algorithms.minors.contracted_nodes(G, ktree + 1, 0) 
 
	return G

def compareAlgebraicConnectivity(tree_info, aux, elem, tol, list_erro):
	"""   Having elem satisfied the tolerance, we calculate the order of T by aux 
		  and then by a straight calculation we obtain the algebraic connectivity of T. 
		  Next, we compare the two values and see if it satisfies the new tolerance. """
	b = 0
	if (tree_info[0] == "odd"):
		b = 1
	# As n - (q + i) = aux, since aux is supposed to be path order.
	#  q + b - i = aux => q = aux + i - b => n = 2 * q + b = 2 * (aux + i) - b
	# Therefore, as we have the order n and k of broom, 
	# we can directly calculate a(T) and compare with elem.
	
	ntree = 2 * (aux + tree_info[2]) - b
	T = createBroom(ntree, tree_info[1])

	if (abs(nx.algebraic_connectivity(T) - elem) < tol * 10):
		list_erro.append((f"T_{ntree},{tree_info[1]}", elem))

def verifyAlgebraicConnectivity(name, tree_info, list_aT,  tol = 1e-08):
	aux = 2
	list_possibles_alg_connectivity = list()

	for elem in list_aT:
		if(not np.iscomplex(elem)):
			while(elem < algebraicConnectivity(aux)):
				aux +=1
			a = calcDif(elem, aux)
			b = calcDif(elem, aux - 1)

			if(a < tol):
				compareAlgebraicConnectivity(tree_info, aux, elem, tol, list_possibles_alg_connectivity)			 
			elif(b < tol):
				compareAlgebraicConnectivity(tree_info, aux-1, elem, tol, list_possibles_alg_connectivity)	

	if(len(list_possibles_alg_connectivity) != 0):
		name.write(f"\t\tRemaining algebraic connectivities: {list_possibles_alg_connectivity}\n")
	