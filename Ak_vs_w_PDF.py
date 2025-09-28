import matplotlib.pyplot as plt
from numpy import loadtxt
import os
import sys
import numpy as np

Nc = 8
U = 7
ts = 1.0
tp = -0.15
ess = [0.25,0.0,-0.25,-0.5,-0.75,-1.0]

n1s = [0.88,0.9,0.92,0.94,0.95,0.97]
n2s = [0.91,0.9,0.89,0.88,0.87,0.87]

V = 0.1
#Vs = [3.0]
T = 0.06
ds = [0.9]

bands = [0,1]
Ks = [2,4]

mode = 'K'

s = ['o','s','^','v','p','h','D','8','<','>','H','o','s','^','v','p','h','D','8','<','>','H']
colors = ['C0', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9', 'b','k','tomato','gold','teal','cyan']
#############################################
fig = plt.subplots(figsize=(8,6),constrained_layout=True)
for iid in range(0,len(ds)): 
	d = ds[iid] 
	for ie in range(0,len(ess)): 
		es = ess[ie]
		for ib in range(len(bands)):
			b = bands[ib]
			title1 = 'Nc'+str(Nc)+'_U'+str(U)+'_Gk_d'+str(d)+'_T'+str(T)+'_Outer'
			title2 = 'Nc'+str(Nc)+'_U'+str(U)+'_Gk_d'+str(d)+'_T'+str(T)+'_Inner'

			if mode == 'K':

				for ik in range(0,len(Ks)):
					dataname = './Nc'+str(Nc)+'_U'+str(U)+'_d'+str(d)+'/OmegaMaxEnt_final_result'+'/optimal_spectral_function_Gw_K'\
					+str(Ks[ik])+'_d'+str(d)+'_es'+str(es)+'_V'+str(V)+'_T'+str(T)+'_'+str(b)+'.dat.dat'
                                
					if os.path.isfile(dataname):
						print(dataname)
						x1,y1=loadtxt(dataname,unpack=True)
						y1 = y1 - ie*1.0
						#print(y1)

					if os.path.isfile(dataname) and b == 0 :
						plt.figure(1)
						if Ks[ik] == 4:
							plt.plot(x1,y1,c=colors[ie],linewidth=2.5,linestyle='-',label='$\delta_{I}=$'+str(round(1-n1s[ie],2))+', $\delta_{O}=$'+str(round(1-n2s[ie],2)))
						else:
							plt.plot(x1,y1,c=colors[ie],linewidth=2.5,linestyle='--',alpha=0.6)
 

						plt.plot([0,0],[-7,2.0],color='gray',linewidth=1.0,linestyle='--')
						plt.plot([-10.5,10.5],[0,0],color='gray',linewidth=1.0,linestyle='--')
						plt.plot([-10.5,10.5],[-1.0,-1.0],color='gray',linewidth=1.0,linestyle='--')
						plt.plot([-10.5,10.5],[-2.0,-2.0],color='gray',linewidth=1.0,linestyle='--')
						plt.plot([-10.5,10.5],[-3.0,-3.0],color='gray',linewidth=1.0,linestyle='--')
						plt.plot([-10.5,10.5],[-4.0,-4.0],color='gray',linewidth=1.0,linestyle='--')
						plt.plot([-10.5,10.5],[-5.0,-5.0],color='gray',linewidth=1.0,linestyle='--')
						#plt.plot([-10.5,10.5],[-6.0,-6.0],color='gray',linewidth=0.8,linestyle='--')

						#plt.title(title1)
						#plt.legend(loc='upper right',fontsize=14)
                	        	        #plt.fill_between([-4,4],0,20,color='gainsboro')
						plt.xlabel("$\\omega$",fontsize=28)
						plt.ylabel('$A_{OL}(K,\\omega)$',fontsize=28)
                               			#plt.ylim(-0.05,3.5)

						plt.xticks(fontsize=22)#,alpha=0);plt.tick_params(bottom=False)
						plt.yticks(fontsize=0,alpha=0);plt.tick_params(left=False)
						plt.tick_params(left=False)
						plt.text(-6,1.2,'(e)',fontsize=36)
						#plt.plot([0,0],[0,5],color='r',linestyle='--')
						#plt.plot([-10.5,10.5],[0,0],color='r',linestyle='--')
						plt.ylim(-5.05,2.0)
						plt.xlim(-7.5,7.5)
						plt.savefig('Aw_K_Nc'+str(Nc)+'_U'+str(U)+'_d'+str(d)+'_V'+str(V)+'_ts'+str(ts)+'_tp'+str(tp)+'_T'+str(T)+'_Outer.pdf',bbox_inches='tight')
                                
					if os.path.isfile(dataname) and b == 1:
						fig = plt.figure(figsize=(8,6),constrained_layout=True)
						plt.figure(2)
						if Ks[ik] == 4:
							plt.plot(x1,y1,c=colors[ie],linewidth=2.5,linestyle='-',label='AN, $\delta_{I}=$'+str(round(1-n1s[ie],2))+', $\delta_{O}=$'+str(round(1-n2s[ie],2)))
						else:
							plt.plot(x1,y1,c=colors[ie],linewidth=2.5,linestyle='--',alpha=0.6)
				
						plt.plot([0,0],[-7,2.0],color='gray',linewidth=1.0,linestyle='--')
						plt.plot([-10.5,10.5],[0,0],color='gray',linewidth=1.0,linestyle='--')
						plt.plot([-10.5,10.5],[-1.0,-1.0],color='gray',linewidth=1.0,linestyle='--')
						plt.plot([-10.5,10.5],[-2.0,-2.0],color='gray',linewidth=1.0,linestyle='--')
						plt.plot([-10.5,10.5],[-3.0,-3.0],color='gray',linewidth=1.0,linestyle='--')
						plt.plot([-10.5,10.5],[-4.0,-4.0],color='gray',linewidth=1.0,linestyle='--')
						plt.plot([-10.5,10.5],[-5.0,-5.0],color='gray',linewidth=1.0,linestyle='--')
						#plt.plot([-10.5,10.5],[-6.0,-6.0],color='gray',linewidth=0.8,linestyle='--')
	
						#plt.title(title1)
                                                #plt.legend(loc='upper right',fontsize=14)
						plt.fill_between([-10,10],[-10,-10],[5,5],color='powderblue',alpha=0.1)
						plt.xlabel("$\\omega$",fontsize=28,alpha=0.0)
						plt.ylabel('$A_{IL}(K,\\omega)$',fontsize=28)
                                                #plt.ylim(-0.05,3.5)

						plt.xticks(fontsize=22,alpha=0);plt.tick_params(bottom=False)
						plt.yticks(fontsize=0,alpha=0);plt.tick_params(left=False)
						plt.tick_params(left=False)
						plt.text(-6,1.2,'(b)',fontsize=36)
                                                #plt.plot([0,0],[0,5],color='r',linestyle='--')
                                                #plt.plot([-10.5,10.5],[0,0],color='r',linestyle='--')
						plt.ylim(-5.05,2.0)
						plt.xlim(-7.5,7.5)

						plt.savefig('Aw_K_Nc'+str(Nc)+'_U'+str(U)+'_d'+str(d)+'_V'+str(V)+'_ts'+str(ts)+'_tp'+str(tp)+'_T'+str(T)+'_Inner.pdf',bbox_inches='tight')
			
