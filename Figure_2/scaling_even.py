#This file prints Sq and S of bipartition as a function of the ratio of hoppings

import numpy as np
from mpmath import mp

def Cij_0(eigvec, Np):
	Cij=eigvec[:,0:Np]*eigvec[:,0:Np].T
	return(Cij)
	
def FreeFermions(subsystem,C):
	C=mp.matrix([[C[x,y] for x in subsystem] for y in subsystem])
	C_eigval=mp.eigh(C, eigvals_only=True)
	EH_eigval=mp.matrix([mp.log(mp.fdiv(mp.fsub(mp.mpf(1.0),x),x)) for x in C_eigval])
	S=mp.re(mp.fsum([mp.log(mp.mpf(1.0)+mp.exp(-x))+mp.fdiv(x,mp.exp(x)+mp.mpf(1.0)) for x in EH_eigval]))
	return(S)


# def FreeFermions(eigvec, subsystem, FermiVector):
	# r=range(FermiVector)
	# Cij=mp.matrix([[mp.fsum([eigvec[i,k]*eigvec[j,k] for k in r]) for i in subsystem] for j in subsystem])
	# C_eigval=mp.eigsy(Cij, eigvals_only=True)
	# EH_eigval=mp.matrix([mp.log(mp.fdiv(mp.fsub(mp.mpf(1.0),x),x)) for x in C_eigval])
	# S=mp.re(mp.fsum([mp.log(mp.mpf(1.0)+mp.exp(-x))+mp.fdiv(x,mp.exp(x)+mp.mpf(1.0)) for x in EH_eigval]))
	# return(S)
	
ClosedLoop = False

#t_vec= np.arange(0.01, 1.5, 0.01)#vector: contains hopping amplitude ratios	

#vector of all possible lenghts for the D part
L_vec=[80,96] #input possible lengths of the chain
window=4
for L in L_vec:
	Np = int(L / 2)
	l=[int(L/4), int(L/2), int(3*L/4), int(L)]
	A = list(range(l[0]))
	B = list(range(l[0],l[1]))
	D = list(range(l[1],l[2]))
	C = list(range(l[2],l[3]))
	
	
	outputfile="./Sq_function_of_hopping_L="+str(L)+".dat"#file in which I will write the results
	f1=open(outputfile,"w")
	
	 
	point1= 1-window/L
	point2= 1+window/L
	v_vec= np.concatenate((np.linspace(0.01,point1-0.01,20),np.linspace(point1,point2,100),np.linspace(point2+0.01,1.5,20)))
			
	mp.dps=L*4 #sets the digits of the decimal numbers
	
	END = L - 1

	if ClosedLoop == True:
		END = L
	#buildHamiltonian
	H=mp.matrix(L)

	for v in v_vec:
		print("Currently working on...")
		print("L="+str(L)) 
		print(v)
		print("................................")
		#the following cicle takes care of the staggered hopping
		for i in range(END):
			j = (i + 1) % L
			if i % 2 == 0:
				H[i, j] = -v
				H[j, i] = H[i, j]
			elif j % 2 == 0:
				H[i, j] = -1
				H[j, i] = H[i, j]
						

		#find its eigenvalues and vectors
		eigval, eigvec=mp.eigsy(H)
		eigval, eigvec=mp.eig_sort(eigval,eigvec)
		
		Cij=Cij_0(eigvec, Np)
			
		#From now on I will use the Free Fermions technique to calculate the entanglement entropies of the subsistems in units of log2
			
			
		SB=FreeFermions(B,Cij)/mp.log(mp.mpf(2.0))
			
		SAB=FreeFermions(A+B,Cij)/mp.log(mp.mpf(2.0))
			
		SBC=FreeFermions(B+C,Cij)/mp.log(mp.mpf(2.0))
						
		SABC=FreeFermions(D,Cij)/mp.log(mp.mpf(2.0)) 
			
			
			
		#In the end I calculate Sqtopo ad proposed by Wen
		Sq=(SAB+SBC-SB-SABC)

		f1.write ("%.5f,%.40f\n" %(v, Sq))#print to file

	f1.close()

