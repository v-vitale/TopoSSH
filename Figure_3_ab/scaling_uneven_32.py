#This file prints Sq and S of bipartition as a function of the ratio of hoppings

import numpy as np
from mpmath import mp


def FreeFermions(eigvec, subsystem, FermiVector):
	r=range(FermiVector)
	Cij=mp.matrix([[mp.fsum([eigvec[i,k]*eigvec[j,k] for k in r]) for i in subsystem] for j in subsystem])
	C_eigval=mp.eigsy(Cij, eigvals_only=True)
	EH_eigval=mp.matrix([mp.log(mp.fdiv(mp.fsub(mp.mpf(1.0),x),x)) for x in C_eigval])
	S=mp.re(mp.fsum([mp.log(mp.mpf(1.0)+mp.exp(-x))+mp.fdiv(x,mp.exp(x)+mp.mpf(1.0)) for x in EH_eigval]))
	return(S)
	
ClosedLoop = False

window=4

#vector of all possible lenghts for the D part
LA_vec=[32]
for LA in LA_vec:
	LD_vec=[96 128 160 192] 
	for LD in LD_vec:
		print("Currently working on...")
		print("LA="+str(LA)+" in ") 
		print(LA_vec)
		
		print("LD="+str(LD)+" in ") 
		print(LD_vec)
		print("................................")
		
		iteration=0 #this counts the iteration on the lengths
		L=3*LA+LD
		point1= 1-window/L
		point2= 1+window/L
		t_vec= np.concatenate((np.linspace(0.01,point1-0.01,30),np.linspace(point1,point2,101),np.linspace(point2+0.01,1.5,30)))
		
		mp.dps = L*3.5  
		Np = int(L / 2)
		A = list(range(LA))
		B = list(range(LA, 2*LA))
		D = list(range(2*LA, 2*LA+LD))
		C = list(range(2*LA+LD, L))

		outputfile="./Sq_function_of_hopping_A="+str(LA)+"_D="+str(LD)+".dat"#file in which I will write the results
		
		D_alternative=A+B+C
		if len(D)>len(D_alternative):
			D=D_alternative

		f1=open(outputfile,"w")
			
		mp.dps=L*4 #sets the digits of the decimal numbers
		
		END = L - 1

		if ClosedLoop == True:
			END = L

		for t in t_vec:
			#buildHamiltonian
			H=mp.matrix(L)
			#the following cicle takes care of the staggered hopping
			for i in range(END):
				j = (i + 1) % L
				if i % 2 == 0:
					H[i, j] = -t
					H[j, i] = H[i, j]
				elif j % 2 == 0:
					H[i, j] = -1
					H[j, i] = H[i, j]
							

			#find its eigenvalues and vectors
			eigval, eigvec=mp.eigsy(H)
			eigval, eigvec=mp.eig_sort(eigval,eigvec)
				
			#From now on I will use the Free Fermions technique to calculate the entanglement entropies of the subsistems in units of log2
				
				
			SB=FreeFermions(eigvec, B, Np)/mp.log(mp.mpf(2.0))
				
			SAB=FreeFermions(eigvec, A+B, Np)/mp.log(mp.mpf(2.0))
				
			SBC=FreeFermions(eigvec, B+C, Np)/mp.log(mp.mpf(2.0))
							
			SABC=FreeFermions(eigvec, D, Np)/mp.log(mp.mpf(2.0)) 
				
				
				
			#In the end I calculate Sqtopo ad proposed by Wen
			Sq=(SAB+SBC-SB-SABC)

			f1.write ("%.5f,%.40f\n" %(t, Sq))#print to file

		f1.close()

