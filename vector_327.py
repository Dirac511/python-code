import numpy as np 
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import math
import os
pi = math.pi

Nc = 8
U = 7
coefp = 1.0

ess = [0.0]
V1s = [0.0]
V2s = [0.0]
ds = [2.7]
Ts = [0.08]

Nc = 4
nOrb = 4


orbitals = [0,1,2,3]

num = [2,3,4]

mode = 'band'
#mode = 'orbital'

for iid in range(0,len(ds)):
    d = ds[iid]
    for iV in range(0,len(V1s)):
        V1 = V1s[iV]
        for ie in range(0,len(ess)):
            es = ess[ie]
            for iT in range(0,len(Ts)):
                T = Ts[iT]

                if mode == 'orbital':
                    if es == 0.0:
                        dataname= './leading_Evec_vs_K_T'+str(T)+'.txt'
                    else:
                        dataname= '../Nc'+str(Nc)+'/coefp1.0/d'+str(d)+'/es'+str(es)+'/V01_'+str(V1)+'_V12_'+str(V2s[iV])+'/leading_Evec_vs_K_T'+str(T)+'.txt'

                    if os.path.exists(dataname):
                        print(dataname)
                        a = np.loadtxt(dataname,skiprows=0,unpack=True)

                        #print(round(a[3,1]/pi,2)) ### check lowest frequency

                        for m in range(0,len(orbitals)):
                            M = orbitals[m]
                            for i in range(0,len(num)):
                                n = num[i]-1
                                N = num[i]
                                k1s = a[2,16*Nc*M*(nOrb+1)+n*Nc:16*Nc*M*(nOrb+1)+n*Nc+Nc]  ##first K
                                k2s = a[3,16*Nc*M*(nOrb+1)+n*Nc:16*Nc*M*(nOrb+1)+n*Nc+Nc]  ##second K
                                print(a[0:2,16*Nc*M*(nOrb+1)+n*Nc:16*Nc*M*(nOrb+1)+n*Nc+Nc])
                                Vec = list(a[5,16*Nc*M*(nOrb+1)+n*Nc:16*Nc*M*(nOrb+1)+n*Nc+Nc])
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

                                if Nc == 8:
                                    K1.extend([-pi,-pi,-pi,0,pi])
                                    K2.extend([pi,0,-pi,-pi,-pi])
                                    Vec.extend([Vec[5],Vec[4],Vec[5],Vec[1],Vec[5]])  ## add else K points
                                elif Nc == 4:
                                    K1.extend([-pi,-pi,-pi,0,pi])
                                    K2.extend([pi,0,-pi,-pi,-pi])
                                    Vec.extend([Vec[3],Vec[2],Vec[3],Vec[1],Vec[3]])  ## add else K points

                                X,Y = np.meshgrid(np.linspace(-pi,pi,300),np.linspace(-pi,pi,300))
                                Z = griddata((K1,K2),Vec,(X,Y),method='cubic')
                                #print(K1)
                                #print(K2)
                                #print(Vec)


                                fig, ax = plt.subplots()
                                t = abs(max(Z.flatten().tolist(),key=abs))
                                contour = ax.contourf(X, Y, Z, levels=np.arange(-t-t/100, t+t/100, t/100), cmap='RdBu')
                                cb = fig.colorbar(contour, ax=ax)
                                ticks = [-t, -t/2, 0.0, t/2, t]
                                cb.set_ticks(ticks)
                                cb.ax.tick_params(labelsize=22)
                                cb.mappable.set_clim(-t-t/100, t+t/100)

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
                                    ax.set_title('$X0$_'+str(N),fontsize=40)
                                    plt.savefig('orbital_'+str(N)+'_Nc'+str(Nc)+'_U'+str(U)+'_d'+str(d)+'_T'+str(T)+'_lowerx.pdf',bbox_inches='tight')
                                elif M == 1:
                                    ax.set_title('$Z1$_'+str(N),fontsize=40)
                                    plt.savefig('orbital_'+str(N)+'_Nc'+str(Nc)+'_U'+str(U)+'_d'+str(d)+'_T'+str(T)+'_lowerz.pdf',bbox_inches='tight')
                                elif M == 2:
                                    ax.set_title('$X2$_'+str(N),fontsize=40)
                                    plt.savefig('orbital_'+str(N)+'_Nc'+str(Nc)+'_U'+str(U)+'_d'+str(d)+'_T'+str(T)+'_upperx.pdf',bbox_inches='tight')
                                elif M == 3:
                                    ax.set_title('$Z3$_'+str(N),fontsize=40)
                                    plt.savefig('orbital_'+str(N)+'_Nc'+str(Nc)+'_U'+str(U)+'_d'+str(d)+'_T'+str(T)+'_upperz.pdf',bbox_inches='tight')
                elif mode == 'band':
                    if es == 0.0:
                        dataname = './leadingspm_Evec_vs_K_T' + str(T) + '.txt'
                    else:
                        dataname = '../Nc' + str(Nc) + '/coefp1.0/d' + str(d) + '/es' + str(es) + '/V01_' + str(V1) + '_V12_' + str(V2s[iV]) + '/check_leadingspm_Evec_vs_K_T' + str(T) + '.txt'
                    if os.path.exists(dataname):
                        print(dataname)
                        a = np.loadtxt(dataname, skiprows=0, unpack=True)
                        #print(a)

                        print(round(a[2,0] / pi, 2))  ### check lowest frequency

                        for m in range(0, len(orbitals)):
                            M = orbitals[m]
                            for i in range(0, len(num)):
                                n = num[i] - 1
                                N = num[i]
                                k1s = a[0, n*Nc:n*Nc+Nc]  ##first K
                                k2s = a[1, n*Nc:n*Nc+Nc]  ##second K
                                print(a[0:2, n*Nc:n*Nc+Nc])
                                Vec = list(a[3+M, n*Nc:n*Nc+Nc])
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

                                if Nc == 8:
                                    K1.extend([-pi,-pi,-pi,0,pi])
                                    K2.extend([pi,0,-pi,-pi,-pi])
                                    Vec.extend([Vec[5],Vec[4],Vec[5],Vec[1],Vec[5]])  ## add else K points
                                elif Nc == 4:
                                    K1.extend([-pi,-pi,-pi,0,pi])
                                    K2.extend([pi,0,-pi,-pi,-pi])
                                    Vec.extend([Vec[3],Vec[2],Vec[3],Vec[1],Vec[3]])  ## add else K points

                                X, Y = np.meshgrid(np.linspace(-pi, pi, 300), np.linspace(-pi, pi, 300))
                                Z = griddata((K1, K2), Vec, (X, Y), method='cubic')
                                # print(K1)
                                # print(K2)
                                #print(Z)

                                fig, ax = plt.subplots()
                                t = abs(max(Z.flatten().tolist(),key=abs))
                                contour = ax.contourf(X, Y, Z, levels=np.arange(-t-t/100, t+t/100, t/100), cmap='RdBu')
                                cb = fig.colorbar(contour, ax=ax)
                                ticks = [-t, -t/2, 0.0, t/2, t]
                                cb.set_ticks(ticks)
                                cb.ax.tick_params(labelsize=22)
                                cb.mappable.set_clim(-t-t/100, t+t/100)

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
                                    ax.set_title('$Bondingd_{x^2-y^2}$_'+str(N), fontsize=40)
                                    plt.savefig('band_' + str(N) + '_Nc' + str(Nc) + '_U' + str(U) + '_d' + str(d) + '_T' + str(T) + '_bondingx.pdf',bbox_inches='tight')
                                elif M == 1:
                                    ax.set_title('$Bondingd_{z^2}$_' + str(N), fontsize=40)
                                    plt.savefig('band_' + str(N) + '_Nc' + str(Nc) + '_U' + str(U) + '_d' + str(d) + '_T' + str(T) + '_bondingz.pdf',bbox_inches='tight')
                                elif M == 2:
                                    ax.set_title('$Anti-bondingd_{x^2-y^2}$_' + str(N), fontsize=40)
                                    plt.savefig('band_' + str(N) + '_Nc' + str(Nc) + '_U' + str(U) + '_d' + str(d) + '_T' + str(T) + '_antibondingx.pdf',bbox_inches='tight')
                                elif M == 3:
                                    ax.set_title('$Anti-bondingd_{z^2}$_' + str(N), fontsize=40)
                                    plt.savefig('band_' + str(N) + '_Nc' + str(Nc) + '_U' + str(U) + '_d' + str(d) + '_T' + str(T) + '_antibondingz.pdf',bbox_inches='tight')
