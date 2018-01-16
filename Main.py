import re
import os, sys
import numpy as np
import time
import numpy.random as rnd
import subprocess

import Hyperbolic as hyp
import Network as Net

import scipy
import scipy.sparse as sparse
from scipy.sparse import csr_matrix,linalg
import matplotlib.pyplot as plt
from math import log, exp
import copy

###########################
def CalculateEntropy (network, L, beta):
	log2=log(2)
	Z=0.0
	S=0.0
	eigVal=L[1]
	for eig in eigVal:
		if eig<0.0:
			eig=0.0
		e=exp(-beta*eig)
		Z+=e
		S+=eig*e
	S=log(Z)/log2+beta/(log2*Z)*S
	S=S/log(network.numNodes)*log2
	print Z,S
	return S
###########################
def RemoveNodeAndCalculateEntropy (network, A, nodes, beta):
	log2=log(2)
	for i in nodes:
		for j in network.nodes[i].link:
			A[i][j.neighbor]=0.0
			A[j.neighbor][i]=0.0
	K = A.sum(0).reshape((1, network.numNodes)).tolist()[0]
	D = sparse.coo_matrix((K,((range(network.numNodes), range(network.numNodes)))), shape=(network.numNodes,network.numNodes)).toarray()
	L = csr_matrix(D-A)
	eigVal,eigVec=sparse.linalg.eigsh(L,k=network.numNodes-1)
	print eigVal
	S=CalculateEntropy(network,(L,eigVal),beta)
	for i in nodes:
		for j in network.nodes[i].link:
			A[i][j.neighbor]=1.0
			A[j.neighbor][i]=1.0
	return S
###########################
def CalculateEntropicDistance (network, A, x, beta):
	a,b=x[0],x[1]
	log2=log(2)
	Sa=RemoveNodeAndCalculateEntropy(network,A,[x[0]],beta)
	Sb=RemoveNodeAndCalculateEntropy(network,A,[x[1]],beta)
	Sab=RemoveNodeAndCalculateEntropy(network,A,[x[0],x[1]],beta)
	if Sa<0.0 or Sb<0.0 or Sab<0:
		print "error"
		exit(1)
	distance=Sa+Sb-Sab
	print distance
###########################
def GetAdjacencyMatrix(network):
	N=network.numNodes
	elements = []
	for i in range (N):
		for j in range (network.degree[i]):
			neighbor = network.nodes[i].link[j].neighbor
			weight = network.nodes[i].link[j].weight
			elements.append((i, neighbor, weight))
	row = [x[0] for x in elements]
	col = [x[1] for x in elements]
	data = [x[2] for x in elements]
	A=sparse.coo_matrix((data,(row, col)), shape=(N,N)).toarray()
	return A
###########################
def GetLaplacian (A):
	N=A.shape[0]
	K = A.sum(0).reshape((1,N)).tolist()[0]
	D = sparse.coo_matrix((K,(range(N),range(N))), shape=(N,N)).toarray()
	L = csr_matrix(D-A)
	eigVal,eigVec=scipy.linalg.eigh(L.A)
	#sparse.linalg.eigsh(L,k=N)
	print eigVal
	return (L, eigVal, eigVec)
###########################

netName = sys.argv[1]
cooName = sys.argv[2]
idr = int(sys.argv[3])
idt = int(sys.argv[4])

t1=time.time()

network = Net.Network(netName)
network.ReadHyperbolicNetwork(netName, cooName, idr, idt)


A=GetAdjacencyMatrix(network)
L=GetLaplacian(A)

#RemoveNodeAndCalculateEntropy (network, A, range(69), 1.0)


'''
distances=[]

for x in elements:
	if x[0]>x[1]:
		distances.append(CalculateEntropicDistance(network, copy.copy(A), x, 1.0))

for a in range(network.numNodes):
	for b in range(network.numNodes):
		if a>b:
			distances.append(CalculateEntropicDistance(network, A[:], [a,b], 10.0))
'''


t2=time.time()
print "running time",(t2-t1)

###########################
