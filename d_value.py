import numpy as np 
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import math
import os
pi = math.pi

Nc = 8
U = 7
coefp = 1.0
ess = [0.25,0.0,-0.25,-0.5,-0.75,-1.0]

n1s = [0.88,0.9,0.92,0.94,0.95,0.97]
n2s = [0.91,0.9,0.89,0.88,0.87,0.87] #V0.1_d0.9

#n1s = [0.83,0.85,0.87,0.89,0.91,0.93]
#n2s = [0.86,0.85,0.84,0.83,0.82,0.81] ##V0.1_d0.85

#n1s = [0.93,0.96]
#n2s = [0.89,0.87]  ## V0.5_d0.9

Vs = [0.5,1.0,1.5,2.0,2.5,3.0]
Vs = [0.1]
d = 0.9
T = 0.06

orbitals = [0,1,2]

num = [1]

mode = 'spm'
mode = 'd'

X1s,X2s,Ys = [],[],[]
for iV in range(0,len(Vs)):
    V = Vs[iV]
    for ie in range(0,len(ess)):
        es = ess[ie]

        if mode == 'd':
            dataname= '../Nc'+str(Nc)+'/coefp'+str(coefp)+'/d'+str(d)+'/es'+str(es)+'/V01_'+str(V)+'_V12_'+str(V)+'/leading_Evec_vs_K_T'+str(T)+'.txt'
            if os.path.exists(dataname):
                print(dataname)
                X1s.append(1-n1s[ie]);X2s.append(1-n2s[ie])
                a = np.loadtxt(dataname,skiprows=0,unpack=True)

                print(round(a[4,1]/pi,2)) ### check lowest frequency
                
                for m in range(0,len(orbitals)):
                    M = orbitals[m]
                    ys = 0
                    for i in range(0,len(num)):
                        n = num[i]-1
                        N = num[i]
                        k1s = a[2,512*M+n*8:512*M+n*8+8]  ##first K
                        k2s = a[3,512*M+n*8:512*M+n*8+8]  ##second K
                        print(a[0:2,512*M+n*8:512*M+n*8+8])
                        Vec = list(a[5,512*M+n*8:512*M+n*8+8])
                        K1,K2 = [],[]

                        for j in range(0,len(k1s)):
                            k1 = k1s[j]
                            if k1-pi>1e-2:
                                k = k1-2.*pi
                            else:
                                k = k1
                            K1.append(k)

                        for j in range(0,len(k2s)):
                            k2 = k2s[j]
                            if k2-pi>1e-3:
                                k2 = k2-2.*pi
                            else:
                                k2 = k2
                            K2.append(k2)  ### shift K points to (-pi,pi)

                        K1.extend([-pi,-pi,-pi,0,pi])
                        K2.extend([pi,0,-pi,-pi,-pi])
                        Vec.extend([Vec[5],Vec[4],Vec[5],Vec[1],Vec[5]])  ## add else K points
                        #print(K1)

                        for i in range(0,len(K1)):
                            ys += abs((math.cos(K1[i])-math.cos(K2[i]))*Vec[i])
                        Ys.append(ys)
                        
    fig = plt.figure()
    ax1 = fig.add_subplot(111)    
    ax1.fill_between([-1,1],[-2,-2],[10,10],color='powderblue',alpha=0.8)             
    ax1.plot(X1s, Ys[1::3], marker='o',color = 'steelblue',markersize=10,label='$IL$')
    ax1.set_xlabel('$\delta_{I}$', fontsize=28)
    ax1.set_ylim(-0.1,4.1)
    ax1.tick_params(axis='y',labelsize=18,labelcolor='steelblue')
    ax1.tick_params(axis='x',labelsize=18)
    ax1.set_xlim(0.0,0.125)
    ax1.set_xticks([0.00,0.02,0.04,0.06,0.08,0.1,0.12])
    ax1.set_ylabel('$D_{IL}$', fontsize=28, color='steelblue')
    ax1.plot([0,0.2],[0,0],color='k',linestyle='--')

    ax2 = ax1.twinx() 
    ax2.plot(X1s, Ys[0::3], marker='s',color = 'goldenrod',markersize=10,label='$OL$')
    ax2.set_ylim(-0.1,4.1)
    ax2.tick_params(axis='y',labelsize=18,labelcolor='goldenrod')
    ax2.set_ylabel('$D_{OL}$', fontsize=28,color='goldenrod')
    #plt.title('dvalue_leadingd_'+str(N)+'_Nc'+str(Nc)+'_U'+str(U)+'_T'+str(T),fontsize=14)
    ax1.legend(loc='center left',bbox_to_anchor=(0.1,0.6),fontsize=24)
    ax2.legend(loc='center left',bbox_to_anchor=(0.1,0.4),fontsize=24)


    ax3 = ax1.twiny()
    ax3.plot(X1s, Ys[0::3], marker='s',color = 'steelblue',markersize=10,alpha=0.0)
    ax3.set_xlim(0.15,0.085)
    ax3.set_xticks([0.15,0.14,0.13,0.12,0.11,0.1,0.09])
    ax3.tick_params(axis='x',labelsize=18)
    ax3.set_xlabel('$\delta_{O}$', fontsize=28)
  
    #plt.title('dvalue_leadingd_'+str(N)+'_Nc'+str(Nc)+'_U'+str(U)+'_T'+str(T),fontsize=14)
    #ax1.legend(loc='lower left',fontsize=18)
    #ax2.legend(loc='lower right',fontsize=18)
    

    
    
    plt.savefig('dvalue_leadingd_'+str(N)+'_Nc'+str(Nc)+'_U'+str(U)+'_d'+str(d)+'_T'+str(T)+'_V'+str(V)+'.pdf',bbox_inches='tight')
    #plt.savefig('dvalue_leadingd_'+str(N)+'_Nc'+str(Nc)+'_U'+str(U)+'_T'+str(T)+'_square.pdf',bbox_inches='tight')
    
    #if M == 0 or M == 2:
    #    plt.savefig('dvalue_leadingd_'+str(N)+'_Nc'+str(Nc)+'_U'+str(U)+'_V'+str(V)+'_d'+str(d)+'_T'+str(T)+'_Outer.pdf',bbox_inches='tight')
    #elif M == 1:
    #    plt.savefig('dvalue_leadingd_'+str(N)+'_Nc'+str(Nc)+'_U'+str(U)+'_V'+str(V)+'_d'+str(d)+'_T'+str(T)+'_Inner.pdf',bbox_inches='tight')


