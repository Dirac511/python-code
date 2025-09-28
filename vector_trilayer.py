import numpy as np 
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import math
import os
pi = math.pi

Nc = 8
U = 7
coefp = 1.0

ess = [-1.0]
Vs = [1.5]
ds = [0.9]
Ts = [0.06]

orbitals = [0,1,2]

num = [1]

mode = 'spm'
#mode = 'd'

for iid in range(0,len(ds)):
    d = ds[iid]
    for iV in range(0,len(Vs)):
        V = Vs[iV]
        for ie in range(0,len(ess)):
            es = ess[ie]
            for iT in range(0,len(Ts)):
                T = Ts[iT]

                if mode == 'd':
                    if es == 0.0:
                        dataname= '../Nc'+str(Nc)+'/coefp1.0/d'+str(d)+'/V01_'+str(V)+'_V12_'+str(V)+'/leading_Evec_vs_K_T'+str(T)+'.txt'
                    else:
                        dataname= '../Nc'+str(Nc)+'/coefp1.0/d'+str(d)+'/es'+str(es)+'/V01_'+str(V)+'_V12_'+str(V)+'/leading_Evec_vs_K_T'+str(T)+'.txt'

                    if os.path.exists(dataname):
                        print(dataname)
                        a = np.loadtxt(dataname,skiprows=0,unpack=True)

                        print(round(a[4,1]/pi,2)) ### check lowest frequency

                        for m in range(0,len(orbitals)):
                            M = orbitals[m]
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

                                X,Y = np.meshgrid(np.linspace(-pi,pi,300),np.linspace(-pi,pi,300))
                                Z = griddata((K1,K2),Vec,(X,Y),method='cubic')
                                #print(K1)
                                #print(K2)
                                #print(Vec)


                                fig,ax = plt.subplots()
                                contour = ax.contourf(X, Y, Z, levels= np.arange(-0.41,0.41,0.001), cmap='RdBu')
                                cb = fig.colorbar(contour,ax=ax)
                                ticks = [-0.4,-0.2,0.0,0.2,0.4]
                                cb.set_ticks(ticks)
                                cb.ax.tick_params(labelsize=22)
                                cb.mappable.set_clim(-0.4,0.4)

                                #ax.set_title('leadingd_'+str(N)+'_Nc'+str(Nc)+'_U'+str(U)+'_V'+str(V)+'_d'+str(d)+'_T'+str(T),fontsize=12)
                                x = [-pi,0,pi]
                                ax.set_xticks(x);ax.set_xticklabels([r'$-\pi$','0',r'$\pi$'],fontsize=26)#,alpha=0)
                                ax.set_yticks(x);ax.set_yticklabels([r'$-\pi$','0',r'$\pi$'],fontsize=26)#,alpha=0)
                                ax.set_xlabel('$K_x$',fontsize=30)#,alpha=0);ax.tick_params(bottom=False)
                                ax.set_ylabel('$K_y$',fontsize=30)#,alpha=0);ax.tick_params(left=False)
                                ax.scatter(K1,K2,color='g',s=100,zorder=12)
                                ax.set_xlim(min(x),max(x))
                                ax.set_ylim(min(x),max(x))

                                if M == 0:
                                    ax.set_title('$OL$_'+str(N),fontsize=40)
                                    plt.savefig('dwave_'+str(N)+'_Nc'+str(Nc)+'_U'+str(U)+'_es'+str(es)+'_V'+str(V)+'_d'+str(d)+'_T'+str(T)+'_lower.pdf',bbox_inches='tight')
                                elif M == 1:
                                    ax.set_title('$IL$_'+str(N),fontsize=40)
                                    plt.savefig('dwave_'+str(N)+'_Nc'+str(Nc)+'_U'+str(U)+'_es'+str(es)+'_V'+str(V)+'_d'+str(d)+'_T'+str(T)+'_inner.pdf',bbox_inches='tight')
                                elif M == 2:
                                    ax.set_title('$OL$_'+str(N),fontsize=40)
                                    plt.savefig('dwave_'+str(N)+'_Nc'+str(Nc)+'_U'+str(U)+'_es'+str(es)+'_V'+str(V)+'_d'+str(d)+'_T'+str(T)+'_upper.pdf',bbox_inches='tight')

                elif mode == 'spm':
                    if es == 0.0:
                        dataname = '../Nc' + str(Nc) + '/coefp1.0/d' + str(d) + '/V01_' + str(V) + '_V12_' + str(V) + '/check_leadingspm_Evec_vs_K_T' + str(T) + '.txt'
                    else:
                        dataname = '../Nc' + str(Nc) + '/coefp1.0/d' + str(d) + '/es' + str(es) + '/V01_' + str(V) + '_V12_' + str(V) + '/check_leadingspm_Evec_vs_K_T' + str(T) + '.txt'
                    if os.path.exists(dataname):
                        print(dataname)
                        a = np.loadtxt(dataname, skiprows=0, unpack=True)

                        print(round(a[3, 1] / pi, 2))  ### check lowest frequency

                        for m in range(0, len(orbitals)):
                            M = orbitals[m]
                            for i in range(0, len(num)):
                                n = num[i] - 1
                                N = num[i]
                                k1s = a[0, n*8:n*8+8]  ##first K
                                k2s = a[1, n*8:n*8+8]  ##second K
                                print(a[0:2, n*8:n*8+8])
                                Vec = list(a[3+M, n*8:n*8+8])
                                #print(Vec)
                                K1, K2 = [], []

                                for j in range(0, len(k1s)):
                                    k1 = k1s[j]
                                    if k1 - pi > 1e-2:
                                        k = k1 - 2. * pi
                                    else:
                                        k = k1
                                    K1.append(k)

                                for j in range(0, len(k2s)):
                                    k2 = k2s[j]
                                    if k2 - pi > 1e-3:
                                        k2 = k2 - 2. * pi
                                    else:
                                        k2 = k2
                                    K2.append(k2)  ### shift K points to (-pi,pi)

                                K1.extend([-pi, -pi, -pi, 0, pi])
                                K2.extend([pi, 0, -pi, -pi, -pi])
                                Vec.extend([Vec[5], Vec[4], Vec[5], Vec[1], Vec[5]])  ## add else K points
                                # print(K1)

                                X, Y = np.meshgrid(np.linspace(-pi, pi, 300), np.linspace(-pi, pi, 300))
                                Z = griddata((K1, K2), Vec, (X, Y), method='cubic')
                                # print(K1)
                                # print(K2)
                                # print(Vec)

                                fig, ax = plt.subplots()
                                contour = ax.contourf(X, Y, Z, levels=np.arange(-0.31, 0.31, 0.001), cmap='RdBu')
                                cb = fig.colorbar(contour, ax=ax)
                                ticks = [-0.3, -0.15, 0.0, 0.15,0.3]
                                cb.set_ticks(ticks)
                                cb.ax.tick_params(labelsize=22)
                                cb.mappable.set_clim(-0.3, 0.3)

                                # ax.set_title('leadingd_'+str(N)+'_Nc'+str(Nc)+'_U'+str(U)+'_V'+str(V)+'_d'+str(d)+'_T'+str(T),fontsize=12)
                                x = [-pi, 0, pi]
                                ax.set_xticks(x);
                                ax.set_xticklabels([r'$-\pi$', '0', r'$\pi$'], fontsize=26)  # ,alpha=0)
                                ax.set_yticks(x);
                                ax.set_yticklabels([r'$-\pi$', '0', r'$\pi$'], fontsize=26)  # ,alpha=0)
                                ax.set_xlabel('$K_x$', fontsize=30)  # ,alpha=0);ax.tick_params(bottom=False)
                                ax.set_ylabel('$K_y$', fontsize=30)  # ,alpha=0);ax.tick_params(left=False)
                                ax.scatter(K1, K2, color='g', s=100, zorder=12)
                                ax.set_xlim(min(x), max(x))
                                ax.set_ylim(min(x), max(x))

                                if M == 0:
                                    ax.set_title('$Bonding$_'+str(N), fontsize=40)
                                    plt.savefig('spmwave_' + str(N) + '_Nc' + str(Nc) + '_U' + str(U) + '_es'+str(es)+'_V' + str(V) + '_d' + str(d) + '_T' + str(T) + '_bonding.pdf',bbox_inches='tight')
                                elif M == 1:
                                    ax.set_title('$Non-bonding$_' + str(N), fontsize=40)
                                    plt.savefig('spmwave_' + str(N) + '_Nc' + str(Nc) + '_U' + str(U) + '_es'+str(es)+'_V' + str(V) + '_d' + str(d) + '_T' + str(T) + '_nonbonding.pdf',bbox_inches='tight')
                                elif M == 2:
                                    ax.set_title('$Anti-bonding$_' + str(N), fontsize=40)
                                    plt.savefig('spmwave_' + str(N) + '_Nc' + str(Nc) + '_U' + str(U) + '_es'+str(es)+'_V' + str(V) + '_d' + str(d) + '_T' + str(T) + '_antibonding.pdf',bbox_inches='tight')
