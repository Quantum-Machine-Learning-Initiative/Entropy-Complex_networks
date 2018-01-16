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
from collections import namedtuple
import subprocess
import time
import numpy.random as rnd
from scipy.special import lambertw, erf, erfinv


#############################################
def CalculateKmin (kbar, gamma):
	return ( kbar*(gamma-2.0)/(gamma-1) )
#############################################
def CalculateC (kbar, T, gamma):
	return ( kbar*np.sin(T*np.pi)/(2.0*T)*np.power( ((gamma-2.0)/(gamma-1.0)),2.0) )
#############################################
def CalculateR (N,C):
	return (2.0*np.log(N/C))
#############################################
def SampleKappa(N, kmin1, gamma1):
	kappa=[None]*N
	for i in range(N):
		kappa[i]=kmin1*np.power(1.0-rnd.random_sample(),1.0/(1.0-gamma1))
	return kappa
#############################################
def SampleTheta(N):
	theta=[None]*N
	for i in range(N):
		theta[i]=2.0*np.pi*rnd.random_sample()
	return theta
#############################################
def SampleConditionalKappa(N, nu, kappa1, kmin1, gamma1, kmin2, gamma2):
	kappa2=[None]*N
	if nu==1:
		for i in range(N):
			kappa2[i]=kmin2*np.power( (kappa1[i]/kmin1), ((1.0-gamma1)/(1.0-gamma2)) )
	elif nu==0:
		for i in range(N):
			kappa2[i]=kmin2*np.power(1.0-rnd.random_sample(),1.0/(1.0-gamma2))
	else:
		for i in range(N):
			phi=-np.log( 1.0- np.power((kmin1/kappa1[i]),(gamma1-1.0)) )
			z=1.0/kmin1* np.power(phi,(nu/(nu-1.0)))*np.power(kappa1[i],-gamma1)
			z= z * ( kmin1*np.power(kappa1[i],gamma1) - np.power(kmin1,gamma1)*kappa1[i] )
			zr= z*rnd.random_sample()
			zr = (nu/(1.0-nu))*lambertw( np.power(zr,((nu-1.0)/nu))/(nu/(1.0-nu)) )
			zr = np.power(zr,(1.0/(1.0-nu)))-np.power(phi,(1.0/(1.0-nu)))
			zr = np.exp(-np.power(zr, (1.0-nu)) )
			zr=kmin2*np.power(1.0-zr, (1.0/(1.0-gamma2)))
			kappa2[i]=zr
	return np.array(kappa2).real
#############################################
def SampleConditionalTheta(N, g, theta1):
	theta2=[None]*N
	if g==1:
		theta2=[i for i in theta1]
	elif g==0:
		for i in range(N):
			theta2[i]=2.0*np.pi*rnd.random_sample()
	else:
		twoPI= 2*np.pi
		sigma0=N/(4.0*np.pi)
		if sigma0>100.0:
			sigma0=100.0
		sigma=sigma0*(1.0/g-1.0)
		for i in range(N):
			l=np.sqrt(2.0)*sigma*erfinv( (-1.0+2.0*rnd.random_sample()) * erf(N/(2*np.sqrt(2)*sigma)) )
			theta2[i]=(theta1[i]+twoPI*l/N)%(twoPI)
	return theta2
#############################################
def ChangeVariablesFromS1ToH2(N, kappa, R, kmin):
	r=[None]*N
	for i in range(N):
		r[i]=R-2*np.log(kappa[i]/kmin)
		if r[i]<0.0:
			r[i]=0.0
	return r
#############################################
def PrintCoordinates(r,theta,kappa,name):
	file = open(name,"w")
	N=len(r)
	for i in range(N):
		print>>file, i, r[i], theta[i], kappa[i]
	file.close()
#############################################
def CreateNetworks(kappa, theta, T, kbar):
	N=len(kappa)
	twoPI= 2*np.pi
	links=[]
	edge=0
	invT=1.0/T
	for i in range(N-1):
		for j in range(i+1,N):
			dTheta=N/(twoPI)*np.abs(np.pi-np.abs(np.pi-np.abs(theta[i]-theta[j])))
			mu=np.sin(T*np.pi)/(twoPI*T*kbar)
			r=dTheta/(mu*kappa[i]*kappa[j])
			if rnd.random_sample()<(1.0/(1.0+np.power(r,invT))):
				links.append((i,j))
				edge+=1
	print "number of links:\t",edge
	return links
#############################################
def PrintNetwork(links,name):
	file = open(name,"w")
	for i in links:
		print>> file, i[0],i[1]
#############################################
def ReadLinks(name):
	links=[]
	edge=0
	file = open(name,"r")
	for row in file:
		if row and not "#" in row:
			i = int(row.split()[0])
			j = int(row.split()[1])
			links.append((i,j))
			edge+=1
	file.close()
	print "number of links:\t",edge
	return links
#############################################
def ReadCoordinates(name):
	coords = []   
	file = open(name,"r")
	for row in file:
		if row and not "#" in row:
			coords.append( (int(row.split()[0]) , float(row.split()[1]) , float(row.split()[2])) )

	file.close()
	coords.sort(key=lambda x:x[0])
	return coords
#############################################
def PlotNetwork(links, r, theta, name):

	fig = plt.figure()
	fig.set_size_inches(5,5)
	plt.rc('text', usetex=True)
	plt.rc('font', size=22, **{'family':'sans-serif','sans-serif':['Helvetica']})
	plt.rcParams['xtick.major.pad'] = 8
	plt.rcParams['ytick.major.pad'] = 8

	ax1 = fig.add_subplot(1,1,1,projection='polar')



	for x in links:
		i = x[0]
		j = x[1]
		ax1.plot([theta[i],theta[j]], [r[i],r[j]],'-',color = 'maroon', linewidth = 0.01, alpha=0.2)

	ax1.plot(theta, r,'o', color = 'orange', markeredgecolor='orange', markersize = 1.5, alpha=0.83)



	ax1.grid(False)

	ax1.spines['polar'].set_visible(False)

	ax1.set_yticklabels([])
	ax1.set_xticklabels([])

	#ax1.set_ylabel('Relative size of the structural set', fontsize=18)
	#ax1.set_xlabel('Relabeling probability', fontsize=18)
	#ax1.locator_params(nbins=1)
	#ax1.set_xlim(-0.05, 1.05)
	#ax1.set_ylim(0.15, 0.55)

	#ax1.legend(loc='upper left', numpoints=1, markerscale=1.0,  prop={'size':11})

	fig.tight_layout()
	fig.savefig(name+'.pdf', bbox_inches = 'tight')
	plt.close(fig)
#############################################
def CreateHyperbolicNetwork(seed,N,gamma,kbar,T, Fast=None, Print=None, Plot=None):

	np.random.seed(seed)

	kmin=CalculateKmin(kbar, gamma)
	print "kmin: ",kmin

	C=CalculateC(kbar, T, gamma)
	print "C: ",C

	R=CalculateR(N,C)
	print "R: ",R

	kappa=SampleKappa(N, kmin, gamma)
	theta=SampleTheta(N)

	r=ChangeVariablesFromS1ToH2(N, kappa[:], R, kmin)

	if Fast:
		PrintCoordinates(r[:],theta[:],kappa[:],"coords.txt")
		subprocess.call(["./HYP_C++.out","coords.txt",str(seed),str(N),str(kbar),str(T)])
		links=ReadLinks("links.coords.txt")
		if Plot:
			PlotNetwork(links[:],r[:],theta[:],"layer")
	else:
		links=CreateNetworks(kappa,theta,T,kbar)
		if Print:
			PrintCoordinates(r[:],theta[:],kappa[:],"coords.txt")
			PrintNetwork(links,"el.txt")
		if Plot:
			PlotNetwork(links,r,theta,"layer")

	return links
#############################################
def Create_GMM_Network (seed,N,gamma,kbar,T,g,nu,Fast=False,Print=False,Plot=False,GMCC=False):

	gamma1=gamma[0]
	gamma2=gamma[1]
	kbar1=kbar[0]
	kbar2=kbar[1]
	T1=T[0]
	T2=T[1]

	np.random.seed(seed)

	kmin1=CalculateKmin(kbar1, gamma1)
	print "kmin1: ",kmin1

	C1=CalculateC(kbar1, T1, gamma1)
	print "C1: ",C1

	R1=CalculateR(N,C1)
	print "R1: ",R1

	kmin2=CalculateKmin(kbar2, gamma2)
	print "kmin2: ",kmin2

	C2=CalculateC(kbar2, T2, gamma2)
	print "C2: ",C2

	R2=CalculateR(N,C2)
	print "R2: ",R2


	kappa1=SampleKappa(N, kmin1, gamma1)

	kappa2=SampleConditionalKappa(N, nu, kappa1[:], kmin1, gamma1, kmin2, gamma2)

	theta1=SampleTheta(N)

	theta2=SampleConditionalTheta(N, g, theta1[:])


	r1=ChangeVariablesFromS1ToH2(N, kappa1[:], R1, kmin1)
	r2=ChangeVariablesFromS1ToH2(N, kappa2[:], R2, kmin2)


	if Fast:
		PrintCoordinates(r1[:],theta1[:],kappa1[:],"coords1.txt")
		subprocess.call(["./HYP_CPP.out","coords1.txt",str(seed+1),str(N),str(kbar1),str(T1)])
		links1=ReadLinks("links.coords1.txt")

		PrintCoordinates(r2[:],theta2[:],kappa2[:],"coords2.txt")
		subprocess.call(["./HYP_CPP.out","coords2.txt",str(seed+2),str(N),str(kbar2),str(T2)])
		links2=ReadLinks("links.coords2.txt")

		if Plot:
			PlotNetwork(links1[:],r1[:],theta1[:],"layer1")
			PlotNetwork(links2[:],r2[:],theta2[:],"layer2")
	else:
		links1=CreateNetworks(kappa1,theta1,T1,kbar1)
		links2=CreateNetworks(kappa2,theta2,T2,kbar2)
		if Print:
			PrintCoordinates(r1[:],theta1[:],kappa1[:],"coords1.txt")
			PrintNetwork(links1,"el1.txt")
			PrintCoordinates(r2[:],theta2[:],kappa2[:],"coords2.txt")
			PrintNetwork(links2,"el2.txt")
		if Plot:
			PlotNetwork(links1[:],r1[:],theta1[:],"layer1")
			PlotNetwork(links2[:],r2[:],theta2[:],"layer2")


	if GMCC:
		subprocess.call(["./GMCC.out","links.coords1.txt","coords1.txt","links.coords2.txt","coords2.txt"])

	return links1, links2
#############################################


