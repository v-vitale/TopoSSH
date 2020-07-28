import numpy as np
from matplotlib.ticker import FixedLocator
import matplotlib.pyplot as plt
from mpmath import mp

def tilde(n):
	out=0
	if n%2==0:
		out=int(n/2)
	else:
		out=int((n+1)/2)
	return (out)
	
def bar(n):
	out=0
	if n%2==0:
		out=int(n/2)
	else:
		out=int((n-1)/2)
	return (out)
	
def s(x):
	if x==0:
		out=0
	else:
		out=x*np.log(x)
	return(out)
	
def FreeFermions(eigvec, subsystem, FermiVector):
    r = range(FermiVector)
    Cij = mp.matrix([[mp.fsum([eigvec[i, k] * eigvec[j, k] for k in r]) for i in subsystem] for j in subsystem])
    C_eigval = mp.eigsy(Cij, eigvals_only=True)
    EH_eigval = mp.matrix([mp.log(mp.fdiv(mp.fsub(mp.mpf(1.0), x), x)) for x in C_eigval])
    S = mp.re(mp.fsum([mp.log(mp.mpf(1.0) + mp.exp(-x)) + mp.fdiv(x, mp.exp(x) + mp.mpf(1.0)) for x in EH_eigval]))
    return (S)


minLen = 8
maxLen = 67
pace = 4
r_vec=[1.5]
L_vec=range(minLen,maxLen,pace)
ClosedLoop=False

for r in r_vec:
	outputfile = "./SD_of_L_r="+mp.nstr(r,2)+".dat"
	f1 = open(outputfile, "w")

	Sq=np.zeros(len(L_vec))
	iteration=0
	R=r
	for L in L_vec:  # L is chain length
		print(R,L)
		mp.dps = L * 4.0  # sets the digits of the decimal numbers

		N=int(L/2)
		Np=N
		A = list(range(int(L / 4)))
		B = list(range(int(L / 4), int(L / 2)))
		D = list(range(int(L / 2), int(3 * L / 4)))
		C = list(range(int(3 * L / 4), L))
		if len(A + B + C) < len(D):
			D = A + B + C
			
		H = mp.matrix(L)
		END = L - 1

		if ClosedLoop:
			END = L


		for i in range(END):
			j = (i + 1) % L
			if i % 2 == 0:
				H[i, j] = -R
				H[j, i] = H[i, j]
			elif j % 2 == 0:
				H[i, j] = -1
				H[j, i] = H[i, j]

        # find its eigenvalues and vectors
		eigval, eigvec = mp.eigsy(H)
		eigval, eigvec = mp.eig_sort(eigval, eigvec)

        # From now on I will use the Free Fermions technique to calculate the entanglement entropies of the subsistems in units of log2

		SB = FreeFermions(eigvec, B, Np) / mp.log(mp.mpf(2.0))

		SAB = FreeFermions(eigvec, A + B, Np) / mp.log(mp.mpf(2.0))

		SBC= FreeFermions(eigvec, B + C, Np) / mp.log(mp.mpf(2.0))

		SABC = FreeFermions(eigvec, D, Np) / mp.log(mp.mpf(2.0))

        # In the end I calculate Sqtopo ad proposed by Wen
		Sq[iteration] = (SAB + SBC- SB - SABC)
		# print(tm, tn, td)
		# print(bm, bn, bd)
		print(Sq[iteration])

		f1.write("%i,%.300f\n" % (L, Sq[iteration]))  # print to file
		iteration += 1

	f1.close()