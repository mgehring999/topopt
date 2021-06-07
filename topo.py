import pyomo.environ as pyo
import numpy as np
import matplotlib.pyplot as plt
from solidspy.preprocesor import rect_grid
from solidspy.uelutil import uel4nquad

def FE(nelx,nely,x,penal):
	KE = localK()
	ndof = 2*(nelx+1)*(nely+1)

	K = np.zeros((ndof,ndof))
	F = np.zeros((ndof,1))
	U = np.zeros((ndof,1))

	for elx in range(1,nelx):
		for ely in range(1,nely):
			n1 = (nely+1)*(elx-1)+ely
			n2 = (nely+1)* elx   +ely
			edof = [2*n1-1, 2*n1, 2*n2-1, 2*n2, 2*n2+1, 2*n2+2, 2*n1+1, 2*n1+2]
			kvec = np.resize(KE,(8,))
			K[edof,edof] += x[elx,ely]**penal*kvec

	#print(K)
	F[3] = -1
	fixeddofs =list(range(1,2*(nely+1),2))
	fixeddofs.append(2*(nelx+1)*(nely+1))
	alldofs = list(range(1,2*(nelx+1)*(nely+1)))
	freedofs = list(set(alldofs) - set(fixeddofs))
	nfreedof = len(freedofs)

	U[freedofs] = np.linalg.solve(np.resize(K[freedofs,freedofs],(nfreedof,nfreedof)),F[freedofs])

	return U

def localK():
	nu = 0.3
	E = 1.
	k = [9999,1/2.-nu/6,1/8.+nu/8,-1/4-nu/12, -1/8.+3*nu/8.,
		-1/4+nu/12, -1/8-nu/8, nu/6, 1/8-3*nu/8]

	KE = [ [k[1], k[2], k[3], k[4], k[5], k[6], k[7], k[8]],
	[k[2], k[1], k[8], k[7], k[6], k[5], k[4], k[3]],
	[k[3], k[8], k[1], k[6], k[7], k[4], k[5], k[2]],
	[k[4], k[7], k[6], k[1], k[8], k[3], k[2], k[5]],
	[k[5], k[6], k[7], k[8], k[1], k[2], k[3], k[4]],
	[k[6], k[5], k[4], k[3], k[2], k[1], k[8], k[7]],
	[k[7], k[4], k[5], k[2], k[3], k[8], k[1], k[6]],
	[k[8], k[3], k[2], k[5], k[4], k[7], k[6], k[1]]]

	return E/(1-nu**2)*np.array(KE)





#mesh = rect_grid(1,1,5,5)
#print(mesh)
#coord = np.array([[-1, -1], [1, -1], [1, 1], [-1, 1]])
#kl = uel4nquad(coord,0.3,1)
#print(kl)

#x = np.array([[.5]*10]*10)
#FE(10,10,x,3)
