import numpy as np
from mpmath import mp
import matplotlib.pyplot as plt


def Cij_0(eigvec, Np): #computes correlation matrix on ground state
	Cij=eigvec[:,0:Np]*eigvec[:,0:Np].T
	return(Cij)

def FreeFermions(subsystem,C_t):#implements free fermion technique by peschel
	C=mp.matrix([[C_t[x,y] for x in subsystem] for y in subsystem])
	C_eigval=mp.eigh(C, eigvals_only=True)
	EH_eigval=mp.matrix([mp.log(mp.fdiv(mp.fsub(mp.mpf(1.0),x),x)) for x in C_eigval])
	S=mp.re(mp.fsum([mp.log(mp.mpf(1.0)+mp.exp(-x))+mp.fdiv(x,mp.exp(x)+mp.mpf(1.0)) for x in EH_eigval]))
	return(S)
	


def rn(sigma): #draws a random number between -sigma, sigma
    ans = (np.random.random()-0.5)*sigma
    return (ans)


#########################################################################################################
#                                           PARAMETERS                                                  #
#########################################################################################################

L = 32 #number of sites

#disorder parameter ---> WE WILL HAVE THE POINTS ON THE DISORDER AXIS =[Win, Win+dW, Win+2dW, ... ,Wfin-dW]
Win=0 #first amplitude of disorder
Wfin=6 #last amplitude of disorder
dW=.2 #resolution of disorder amplitude

W_vec = np.arange(Win,Wfin,dW) #all disorder amplitude, in a vector

#intracell hopping parameter ---> WE WILL HAVE THE POINTS ON THE RATIO AXIS =[tin, tin+dt, tin+2dt, ... ,tfin-dt]
tin=0.001 #first intracell hopping amplitude
tfin=2 #last intracell hopping amplitudes
dt=.1 #resolution of intracell hopping amplitude

t_vec = np.arange(tin,tfin,dt) #all intracell hopping amplitudes, in a vector

repetitions=50 #number of repetitions

#########################################################################################################
#                                          END PARAMETERS                                               #
#########################################################################################################



#OPEN FILE
outputfile = "./BUMP_SD_length="+str(L)+"_W="+str(Win)+"-"+str(Wfin)+"_t="+str(tin)+"-"+str(tfin)+"_rep="+str(repetitions)+".dat"
f1 = open(outputfile, "w")






mp.dps = L*3.5  #number of digits
Np = int(L / 2) #number of particles

l=[int(L/4), int(L/2), int(3*L/4), int(L)] #where are the cuts located

A = list(range(l[0])) #sites in part A (WEN NOTATION)
B = list(range(l[0],l[1])) #sites in part B (WEN NOTATION)
D = list(range(l[1],l[2])) #sites in part D (WEN NOTATION)
C = list(range(l[2],l[3])) #sites in part C (WEN NOTATION)
	


H = mp.matrix(L) #1body Hamiltonian initialisation

END=L-1 #SET TO L FOR PBC

for t in t_vec: #LOOP ON HOPPING AMPLITUDES
	for s in W_vec: #LOOP ON DISORDER AMPLITUDES
		Sq_vec=mp.matrix(1,repetitions) #vector in which I will write all SD corresponding to given s and t
		
		for r in range(repetitions): #LOOP ON REPETITIONS
			print(t,s,r) #if an error occurs this will tell me where the loops have stopped
			
			#BUILD THE HAMILTONIAN
			for i in range(END):
				j = (i + 1) % L #j is nearest neighbor of i
				if i % 2 == 0:#intracell case
					H[i, j] = -t + rn(s)
					H[j, i] = H[i, j]
				elif j % 2 == 0: #intercell case
					H[i, j] = -1 + rn(s/2)
					H[j, i] = H[i, j]

            # find its eigenvalues and vectors

			eigval, eigvec = mp.eigh(H)
			eigval, eigvec = mp.eig_sort(eigval, eigvec)
			
			#build correlation matrix
			Cij=Cij_0(eigvec,Np)
			
			#use free fermions to obtain the entanglement entropies of the
			SB=FreeFermions(B,Cij)/mp.log(mp.mpf(2.0))
			
			SAB=FreeFermions(A+B,Cij)/mp.log(mp.mpf(2.0))
			
			SBC=FreeFermions(B+C,Cij)/mp.log(mp.mpf(2.0))
						
			SABC=FreeFermions(D,Cij)/mp.log(mp.mpf(2.0)) 

            # In the end I calculate Sqtopo ad proposed by Wen
			Sq_vec[r] = (SAB + SBC - SB - SABC)
		#average of vector above will give me the average of several realisation of disorder with given (s,t)	
		Sq=np.mean(Sq_vec)
		
		# print to file
		if s==W_vec[-1]: #if all s have been calculated go accapo to start with new line and s=Win
			f1.write("%.60f\n" % (Sq))
		else: #otherwise if s in not the last of its vector continue on the same line
			f1.write("%.60f," % (Sq))  

f1.close()