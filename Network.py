import re
import os, sys
import numpy as np
import matplotlib as mpl
#mpl.use('Agg')
import matplotlib.pyplot as plt
import collections
import operator
from operator import itemgetter
from scipy import stats
import scipy.sparse as sparse
import numpy.linalg as linalg
import numpy.random as rnd

###########################

class Link:

	def __init__ (self, neighbor, weight=1.0):
		self.neighbor = neighbor
		self.weight = weight
		self.state = True

	def SetID (self, id):
		self.id = id

	def Deactivate (self):
		self.state = False

	def Activate (self):
		self.state = True

########################

class Node:

	def __init__ (self, coordinates=None):
		self.link = []
		self.degree = 0
		self.state = True
		self.coordinates = None

	def AddLink (self, link):
		self.link.append(link)
		self.degree += 1

	def Deactivate (self):
		self.state = False

	def Activate (self):
		self.state = True

###########################

class Network:

	#############################################
	def __init__ (self, name):
		self.name=name
		self.nodes=[]
		self.coordinates=[]
		self.links=[]
		self.degree=[]
		self.numNodes=0
		self.numLinks=0
		self.adjacency = None
		self.laplacian = None
		self.RWlaplacian = None
	#############################################
	def ActivateNodesAndLinks(self):
		for x in self.nodes:
			x.Activate()
			for y in x.link:
				y.Activate()
	#############################################
	def AddNode(self, node):
		self.numNodes += 1
		self.nodes.append(node)
	#############################################
	def findNodeID(self, list, id):
		idsplit = 0
		top = len(list)
		if top<1:
			return -1
		bottom = 0
		split = top/2

		while top>bottom:
			idsplit = list[split]
			if id>idsplit :
				bottom = split + 1
				split = (top+bottom)/2
			elif id<idsplit :
				top = split
				split = (top+bottom)/2
			else:
				return split
		return -1;
	#############################################
	def CreateNodes(self):
		for i in range(len(self.coordinates)):
			node = Node(self.coordinates[i][0])
			self.AddNode(node)
		print 'Number of nodes: ' , self.numNodes
	#############################################
	def ReadHyperbolicLinks(self, name):
		self.links = []
		file = open(name,"r")
		numLinks=0
		for row in file:
			if row and not "#" in row:
				i = int(row.split()[0])
				j = int(row.split()[1])
				self.links.append((i,j))
				numLinks += 1
		file.close()
		self.numLinks = numLinks
		print 'Number of Links:' , numLinks 
	#############################################
	def ReadLinks(self, name):
		id = set()
		self.links = []

		# read edgelist file.
		file = open(name,"r")
		for row in file:
			if row and not "#" in row:
				i = int(row.split()[0])
				j = int(row.split()[1])
				self.edgeList.append((i,j))
				id.add(i)
				id.add(j)
		file.close()

		self.nodeIDs = list(id)
		self.nodeIDs.sort()
	#############################################
	def CreateLinks(self):
		self.numLinks=0
		nodeIDs=[x[0] for x in self.coordinates]
		for x in self.links:
			i,j = x[0], x[1]
			i = self.findNodeID(nodeIDs, i)
			j = self.findNodeID(nodeIDs, j)
			if i>=0 and j>=0:
				link = Link(j)
				self.node[i].AddLink(link)
				link = Link(i)
				self.node[j].AddLink(link)
				self.numLinks += 1
			else:
				print 'Wrong node id detection.'

		self.degree = [x.degree for x in self.node]
		print 'Number of links: ' , self.numLinks
	#############################################
	def ReadNetwork(self, name1):
		

		CreateNodes()
		CreateLinks()
	#############################################
	def CreateNetworkFromList (self, edgelist):
		id=set()
		self.edgeList=edgelist[:]

		#Detect node IDs.
		for edge in self.edgeList:
			id.add(edge[0])
			id.add(edge[1])

		self.nodeIDs = list(id)
		self.nodeIDs.sort()

		# create nodes.
		for i in range(len(self.nodeIDs)):
			node = Node(self.nodeIDs[i])
			self.AddNode(node)
		
		print 'Number of nodes: ' , self.numNodes

		# create links.
		for edge in self.edgeList:
			i = edge[0]
			j = edge[1]

			i = self.findNodeID(self.nodeIDs, i)
			j = self.findNodeID(self.nodeIDs, j)

			if i>=0 and j>=0:
				link = Link(j)
				self.nodes[i].AddLink(link)
				link = Link(i)
				self.nodes[j].AddLink(link)
				self.numLinks += 1
			else:
				print 'Wrong node id detection.'

		self.degree = [x.degree for x in self.nodes]
	#############################################
	def ReadHyperbolicNetwork(self, linkName, cooName, idr, idt):
		#1.Read coordinates.
		self.nodes=[]
		self.coordinates=[]
		self.links=[]
		self.degree=[] 
		self.numNodes = 0
		self.numLinks = 0
		file = open(cooName,"r")
		for row in file:
			if row and not "#" in row:
				self.coordinates.append( (int(row.split()[0]) , float(row.split()[idr]) , float(row.split()[idt])) )
		file.close()
		self.coordinates.sort(key=lambda x:x[0])

		#2.Create nodes and assign their IDs and coordinates.
		for i in range(len(self.coordinates)):
			node = Node(self.coordinates[i])
			self.numNodes += 1
			self.nodes.append(node)
		print 'Number of nodes: ' , self.numNodes

		#3.Read links.
		file = open(linkName,"r")
		for row in file:
			if row and not "#" in row:
				i = int(row.split()[0])
				j = int(row.split()[1])
				self.links.append((i,j))
		file.close()	
		nodeIDs=[x[0] for x in self.coordinates]
		for x in self.links:
			i,j = x[0], x[1]
			i = self.findNodeID(nodeIDs, i)
			j = self.findNodeID(nodeIDs, j)
			if i>=0 and j>=0:
				link = Link(j)
				self.nodes[i].AddLink(link)
				link = Link(i)
				self.nodes[j].AddLink(link)
				self.numLinks += 1
			else:
				print 'Wrong node id detection.'
		self.degree = [x.degree for x in self.nodes]
		print 'Number of links: ' , self.numLinks
	#############################################
	def GetAdjacencyMatrix (self):
		elements = []
		for i in range (self.numNodes):
			for j in range (self.degree[i]):
				neighbor = self.nodes[i].link[j].neighbor
				weight = self.nodes[i].link[j].weight
				elements.append( (i, neighbor, weight) )
		row = [x[0] for x in elements]
		col = [x[1] for x in elements]
		data = [x[2] for x in elements]
		self.adjacency = sparse.coo_matrix((data,(row, col)), shape=(self.numNodes,self.numNodes)).toarray()
		return self.adjacency
	#############################################
	def GetLaplacianMatrix (self):
		if self.adjacency is None:
			self.GetAdjacencyMatrix()
		D = np.diag(self.degree)
		self.laplacian = D-self.adjacency
		return self.laplacian		
	#############################################
	def GetRandomWalkNormalizedLaplacianMatrix (self):

		if self.RWlaplacian is not None:
			return self.RWlaplacian

		else:

			if self.adjacency is None:
				self.GetAdjacencyMatrix()

			self.RWlaplacian = np.zeros((self.numNodes,self.numNodes))

			for i in range(self.numNodes):
				for j in range(self.numNodes):
					d_inv = np.power(float(self.degree[i]),-1)
					self.RWlaplacian[i,j] = float(d_inv*self.adjacency[i,j])
			
			self.RWlaplacian = np.identity(self.numNodes)-self.RWlaplacian
			return self.RWlaplacian	
	#############################################
	def PlotNetwork(self):
		
		fig = plt.figure()
		fig.set_size_inches(5,5)
		plt.rc('text', usetex=True)
		plt.rc('font', size=22, **{'family':'sans-serif','sans-serif':['Helvetica']})
		plt.rcParams['xtick.major.pad'] = 8
		plt.rcParams['ytick.major.pad'] = 8

		ax1 = fig.add_subplot(1,1,1,projection='polar')

		#########################

		# seperate Id, Theta and R lists.
		ID = [x[0] for x in self.coordinates]
		R = [x[1] for x in self.coordinates]
		Theta = [x[2] for x in self.coordinates]


		for x in self.edgeList:
			i = x[0]
			j = x[1]

			i = self.findNodeID(self.nodeIDs,i)
			j = self.findNodeID(self.nodeIDs,j)
			
			ax1.plot([Theta[i],Theta[j]], [R[i],R[j]],'-',color = 'maroon', linewidth = 0.01, alpha=0.2)

		ax1.plot(Theta, R,'o',color = 'orange', markeredgecolor='orange', markersize = 1.5, alpha=0.83)


		ax1.grid(False)
		ax1.spines['polar'].set_visible(False)

		ax1.set_yticklabels([])
		ax1.set_xticklabels([])

		fig.tight_layout()
		fig.savefig(self.name+'fig_graph.pdf', bbox_inches = 'tight')
		plt.close(fig)
	#############################################
	def FindComponents(self):

		visited = [False]*self.numNodes
		queue = [0]*self.numNodes

		size = 0
		gcc  = 0

		for i in range(self.numNodes):

			if visited[i] is False:

				read  = 0
				queue[read] = i
				visited[i] = True
				write = 1
				size  = 1

				while read != write:

					n = queue[read]
					read +=1

					for x in self.nodes[n].link:
						m = x.neighbor
						if visited[m] is False:
							queue[write] = m
							visited[m] = True
							write += 1
							size  += 1

			if size>gcc:
				gcc=size

		return gcc
	#############################################
	def FindCore (self):

		condition=True
		
		# Initially set all nodes to 'true' state.
		for node in self.nodes:
			node.state=True

		# remove leaves.
		while condition:
			condition = False
			for node in self.nodes:	
				if node.state:
					d=0
					for x in node.link:
						neighbor = x.neighbor
						if self.nodes[neighbor].state:
							d+=1
					if d<2:
						condition = True
						node.state = False
						for x in node.link:
							self.nodes[x.neighbor].state = False

		core=[i for i in range(self.numNodes) if self.nodes[i].state]

		return core
	#############################################
	def CreatePoissonNetwork (self, N, K):
		p = K/(N-1.0)

		for i in range(N):
			node = Node(i)
			self.AddNode(node)

		numLinks=0


		for i in range(N-1):
			for j in range (i+1, N):
				if(rnd.random_sample()<p):
					link = Link(j)
					self.nodes[i].AddLink(link)
					link = Link(i)
					self.nodes[j].AddLink(link)
					numLinks += 1

		self.degree = [0]*N

		for i in range (N):
			self.degree[i] = self.nodes[i].degree

		self.numLinks = numLinks
	#############################################

###########################


