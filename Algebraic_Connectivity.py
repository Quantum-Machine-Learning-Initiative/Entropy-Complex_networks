import re
import os, sys
import numpy as np
import time
import numpy.random as rnd
from scipy.sparse import csr_matrix, eye
import scipy.sparse as sparse
from scipy.sparse import linalg
import numpy as np
import matplotlib.pyplot as plt

####################################################################################

def GetAdjacencyMatrix(N,links):			
	row=[x[0] for x in links]
	col=[x[1] for x in links]
	data=[x[2] for x in links]
	adjacency=csr_matrix((data,(row, col)), shape=(N,N)) #.toarray()
	adjacency=adjacency+adjacency.transpose()
	return adjacency

def GetDegreeMatrix(A):
	N=A.shape[0]
	K=A.sum(axis=0).tolist()[0]
	D=csr_matrix( (K, (range(N),range(N)) ),shape=(N,N))
	return D

def GetLaplacianMatrix(A,D):
	L=csr_matrix(D-A)
	return L

def GetSupraLaplacian(N,links1,links2,p):
	links=[]
	for x in links1:
		links.append(x)
	for x in links2:
		links.append((x[0]+N,x[1]+N,x[2]))
	for i in range(N):
		links.append((i,i+N,p))
	supraA=GetAdjacencyMatrix(2*N,links)
	supraD=GetDegreeMatrix(supraA)
	supraL=GetLaplacianMatrix(supraA,supraD)
	return supraL

def ReplicateArenasPRL():
	links1=[(0,2),(0,3),(1,2),(1,5),(2,3),(2,4),(2,5),(3,4),(3,5)]
	links2=[(0,1),(0,4),(1,2),(1,4),(1,5),(2,5),(3,5)]

	links1=[(x[0],x[1],1.0) for x in links1]
	links2=[(x[0],x[1],1.0) for x in links2]

	print links1
	print links2

	weight=np.arange(0.05,100.0,0.1)

	results=[]
	for p in weight:
		L=GetSupraLaplacian(6,links1,links2,p)
		val,vec=sparse.linalg.eigsh(L,k=11)
		results.append(val.real)
		


	fig=plt.figure()
	ax1=fig.add_subplot(1,1,1)

	for i in range(11):
		x=[eig[i] for eig in results]
		ax1.plot(weight,x)
		
	ax1.set_xscale("log")
	ax1.set_yscale("log")
	ax1.set_ylim(0.1,100)
	plt.show()


ReplicateArenasPRL()

