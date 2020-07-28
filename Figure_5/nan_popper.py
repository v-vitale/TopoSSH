from mpmath import mp
import numpy as np

tin=0.1
dt=0.1
tfin=2

Win=0 #first amplitude of disorder                                                                                                                                                                         
Wfin=6 #last amplitude of disorder                                                                                                                                                                         
dW=.2 #resolution of disorder amplitude                                                                                                                                                                    

L=64

t_vec = [mp.mpf(0.001)]+list(mp.arange(tin,tfin,dt))

W_vec = mp.arange(Win,Wfin,dW)


outputfile="./BUMP_SD_length="+str(L)+".dat"
f1=open(outputfile,"w")


for t in t_vec:
  
  for s in W_vec:
    print(t,s)
    inputfile="./BUMP_SD_length="+str(L)+"_t="+mp.nstr(t,3)+"_s="+mp.nstr(s,3)+".dat"
    ND=np.genfromtxt(inputfile)#,dtype=float, max_rows=250)
    ND = ND[~np.isnan(ND)]
    value=np.mean(ND)
    st=np.std(ND)
    f1.write("%.10f,%.10f,%.10f,%.10f,%i\n" % (t,s,value,st,len(ND)))

f1.close()
