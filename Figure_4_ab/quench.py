import numpy as np
from mpmath import mp
import matplotlib.pyplot as plt


def Construct_K(Phi,Psi,Omega,tt):
    
    Cos=mp.matrix([[mp.cos(eigvalst[i]*tt) if i==j else mp.mpf(0.0) for j in range(L)]for i in range(L)])
    Sin=mp.matrix([[mp.sin(eigvalst[i]*tt) if i==j else mp.mpf(0.0) for j in range(L)]for i in range(L)])
    
    OmegaCos=Cos*Omega*Cos
    OmegaSin=Sin*Omega*Sin

    Kmatrix=-Phi*(OmegaCos+OmegaSin.T)*Psi.T
    return(Kmatrix)

def Construct_Dp(Phi,Psi,Omega,tt):
    
    Cos=mp.matrix([[mp.cos(eigvalst[i]*tt) if i==j else mp.mpf(0.0) for j in range(L)]for i in range(L)])
    Sin=mp.matrix([[mp.sin(eigvalst[i]*tt) if i==j else mp.mpf(0.0) for j in range(L)]for i in range(L)])
    
    OmegaCosSin=Cos*Omega*Sin

    Dp=-Phi*(OmegaCosSin.T-OmegaCosSin)*Phi.T
    return(Dp)

def Construct_Dm(Phi,Psi,Omega,tt):
    
    Cos=mp.matrix([[mp.cos(eigvalst[i]*tt) if i==j else mp.mpf(0.0) for j in range(L)]for i in range(L)])
    Sin=mp.matrix([[mp.sin(eigvalst[i]*tt) if i==j else mp.mpf(0.0) for j in range(L)]for i in range(L)])
    
    OmegaSinCos=Sin*Omega*Cos

    Dm=-Psi*(OmegaSinCos-OmegaSinCos.T)*Psi.T
    return(Dm)

def SvN(Dm,Dp,Kpart):

    dim=2*len(Kpart)
    gamma=mp.matrix(dim,dim)
    KpartH=Kpart.H
    for m in list(range(0,dim,2)):
        for n in list(range(0,dim,2)):
            row=int((m+1)/2)
            col=int((n+1)/2)
            gamma[m,n]=mp.j*Dp[row,col]
            gamma[m,n+1]=mp.j*Kpart[row,col]
            gamma[m+1,n]=-mp.j*KpartH[row,col]
            gamma[m+1,n+1]=mp.j*Dm[row,col]
            
    eigvalgamma=mp.eigh(gamma ,eigvals_only = True)
    Spart=mp.fsum([-((mp.mpf(1.0)+x)/mp.mpf(2))*mp.log((mp.mpf(1.0)+x)/mp.mpf(2)) for x in eigvalgamma])

    return(Spart)

######################################################################################################
Lrange=[32]


for L in Lrange:
    mp.dps=3.0*L
    l=int(L/4)
    la=list(range(l))
    lb=list(range(la[-1]+1,la[-1]+l+1))
    ld=list(range(lb[-1]+1,L-l))
    lc=list(range(ld[-1]+1,ld[-1]+l+1))
    iter=0
    #Parameters at t=0-
    w=mp.mpf(0.1)
    delta0=mp.mpf(0.0)
    #Parameters at t=0+
    v=1.5
    deltat=mp.mpf(0.0)
    tsteps=mp.arange(0.0,50.01,0.5)
        
    Sqtopo=[]
    Sab=[]
    Sbc=[]
    Sb=[]
    Sabc=[]
        
    pvp=np.zeros(L)
    pvm=np.zeros(L)
        
    for i in list(range(L)):
        pvp[i]=(i+1)%L 
        pvm[(i+1)%L]=i
            
            
#        outputfile="./SqtopoOBC_mpdps_"+str(mp.dps)+"_L_"+str(L)+"_l_"+str(l)+"_mu_"+mp.nstr(chemical0)+"_muf_"+mp.nstr(chemicalt)+"_d_"+mp.nstr(delta0)+"_dt_"+mp.nstr(deltat)+".dat"
#        f1=open(outputfile,"w")
    
    #### t=0- #####
    ClosedLoop = False  # Set true for a closed ring topology, False for a open loop topology

    A0 = mp.matrix(L)
    At = mp.matrix(L)

    END = L-1
    if ClosedLoop == True:
        END = L

    for i in range(END):
        j = (i + 1) % L
        if i % 2 == 0:
            A0[i, j] = -w
            A0[j, i] = A0[i, j]
        elif j % 2 == 0:
            A0[i, j] = -1
            A0[j, i] = A0[i, j]
    
    B0=mp.matrix([[delta0 if i==j-1 else -delta0 if i==j+1 else mp.mpf(0.0) for j in list(range(L))] for i in list(range(L))])        
        
    EP0, eigvalsi, EF0  = mp.svd_r(A0+B0, compute_uv=True)
    eigvecf0=EF0.T
    eigvecp0=EP0
        
    #### t>0 #####
     
    for i in range(END):
        j = (i + 1) % L
        if i % 2 == 0:
            At[i, j] = -v
            At[j, i] = At[i, j]
        elif j % 2 == 0:
            At[i, j] = -1
            At[j, i] = At[i, j]
            
    Bt=mp.matrix([[deltat if i==j-1 else -deltat if i==j+1 else mp.mpf(0.0) for j in list(range(L))] for i in list(range(L))])        
                    
    EPt, eigvalst, EFt  = mp.svd_r(At+Bt, compute_uv=True)
    eigvecft=EFt.T
    eigvecpt=EPt
        
    KK0 =eigvecf0*eigvecp0.T
    Omega=(eigvecft.T)*KK0*eigvecpt  
        
    for tt in tsteps: 
        print("timestep -> ",tt)
                   
        KKt=Construct_K(eigvecft,eigvecpt,Omega,tt)
        Dmt=Construct_Dm(eigvecft,eigvecpt,Omega,tt)
        Dpt=Construct_Dp(eigvecft,eigvecpt,Omega,tt)
            
            
        Kabt = mp.matrix([[KKt[i,j] for j in la+lb] for i in la+lb])
        Dmabt = mp.matrix([[Dmt[i,j] for j in la+lb] for i in la+lb])
        Dpabt = mp.matrix([[Dpt[i,j] for j in la+lb] for i in la+lb])
            
        Kbct = mp.matrix([[KKt[i,j] for j in lb+lc] for i in lb+lc])
        Dmbct = mp.matrix([[Dmt[i,j] for j in lb+lc] for i in lb+lc])
        Dpbct = mp.matrix([[Dpt[i,j] for j in lb+lc] for i in lb+lc])
                        
        Kbt = mp.matrix([[KKt[i,j] for j in lb] for i in lb])
        Dmbt = mp.matrix([[Dmt[i,j] for j in lb] for i in lb])
        Dpbt = mp.matrix([[Dpt[i,j] for j in lb] for i in lb])
            
        Kabct = mp.matrix([[KKt[i,j] for j in la+lb+lc] for i in la+lb+lc])
        Dmabct = mp.matrix([[Dmt[i,j] for j in la+lb+lc] for i in la+lb+lc])
        Dpabct = mp.matrix([[Dpt[i,j] for j in la+lb+lc] for i in la+lb+lc])
            
        Sab.append(SvN(Dmabt,Dpabt,Kabt))
        Sbc.append(SvN(Dmbct,Dpbct,Kbct))
        Sb.append(SvN(Dmbt,Dpbt,Kbt))
        Sabc.append(SvN(Dmabct,Dpabct,Kabct))    
            
            
        print('step '+str(iter+1))
        temp=(Sab[iter]+Sbc[iter]-Sb[iter]-Sabc[iter])/mp.log(2)
        Sqtopo.append(mp.re(temp))
        print("--- Sqtopo ---")
        print(temp)
        print("======================================================")
 #           f1.write("%.10f,%.10f\n" % (tsteps[iter],Sqtopo[iter]))
            
        iter+=1
        
 #      f1.close()
    print('end')
    #plt.plot(tsteps,Sqtopo)
    #plt.title("$\Delta=0.5$")
    #plt.xlabel("t")
    #plt.ylabel("$S^q_{topo} $")
        
    
        
    #plt.legend(loc="upper right")
    #plt.show()
