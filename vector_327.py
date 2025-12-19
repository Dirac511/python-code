import numpy as np 
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import math
import os
pi = math.pi

Nc = 8
U = 3


ds = [2.4]
Ts = [0.08]

Nc = 4
nOrb = 4

index = ['00','01','02','03','10','11','12','13','20','21','22','23','30','31','32','33']
        #  0 ,  1 ,  2 ,  3 ,  4 ,  5 ,  6,,  7 ,  8 ,  9 , 10 , 11 , 12 , 13 , 14 , 15

orbitals = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]
#orbitals = [7]

num = [1,2,3,4,5,6]

mode = 'band'
mode = 'orbital'

for iid in range(0,len(ds)):
    d = ds[iid]
    for iT in range(0,len(Ts)):
        T = Ts[iT]

        if mode == 'orbital':

            dataname= './leading_Evec_vs_K_T'+str(T)+'.txt'

            if os.path.exists(dataname):
                print(dataname)
                a = np.loadtxt(dataname,skiprows=0,unpack=True)

                #print(round(a[3,1]/pi,2)) ### check lowest frequency
                if len(orbitals)==16:
                     fig, axes = plt.subplots(4,4,figsize=(10,10))
                     #plt.subplots_adjust(wspace=-0.6)
                     for i in range(0,len(num)):
                         n = num[i]-1
                         N = num[i]
                         R = []
                         for idx, m in enumerate(orbitals):
                             M = orbitals[m]
                             k1s = a[2,32*Nc*M+n*Nc:32*Nc*M+n*Nc+Nc]  ##first K
                             k2s = a[3,32*Nc*M+n*Nc:32*Nc*M+n*Nc+Nc]  ##second K
                             print(a[0:2,32*Nc*M+n*Nc:32*Nc*M+n*Nc+Nc])
                             Vec = list(a[5,32*Nc*M+n*Nc:32*Nc*M+n*Nc+Nc])

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
                             R.append(Z)
                             #print(K1)
                             #print(K2)
                             #print(Vec)
                         r = np.concatenate([z.flatten() for z in R])
                         t = np.max(np.abs(r))
                    
                         for idx,M in enumerate(orbitals):
                             u = idx // 4
                             v = idx % 4
                             ax = axes[u,v]
                             Z = R[idx]

                             ax.set_title(str(u)+str(v)+'_['+str(N)+']',fontsize=14)

                             contour = ax.contourf(X, Y, Z, levels=np.arange(-t-t/100, t+t/100, t/100), cmap='RdBu')
                             ax.set_xticks([])
                             ax.set_yticks([])
                             x = [-pi,0,pi]                          
                             if idx == 15:
                                  ax.set_xticks(x);ax.set_xticklabels([r'$-\pi$','0',r'$\pi$'],fontsize=14)#,alpha=0)
                                  ax.set_yticks(x);ax.set_yticklabels([r'$-\pi$','0',r'$\pi$'],fontsize=14)#,alpha=0)
                              #ax.set_xlabel('$K_x$',fontsize=20)#,alpha=0);ax.tick_params(bottom=False)
                              #ax.set_ylabel('$K_y$',fontsize=20)#,alpha=0);ax.tick_params(left=False)
                             ax.scatter(K1,K2,color='g',s=100,zorder=12)
                             ax.set_xlim(min(x),max(x))
                             ax.set_ylim(min(x),max(x))

                         fig.subplots_adjust(bottom=0.12)
                         cb = fig.add_axes([0.15, 0.05, 0.7, 0.03])
                         cb = fig.colorbar(contour, cax=cb,orientation='horizontal')
                         ticks = [-t, -t/2, 0.0, t/2, t]
                         cb.set_ticks(ticks)
                         cb.ax.tick_params(labelsize=14)
                         cb.mappable.set_clim(-t-t/100, t+t/100)

                         plt.savefig('orbital_'+str(N)+'_Nc'+str(Nc)+'_U'+str(U)+'_d'+str(d)+'_T'+str(T)+'_total'+'.png',bbox_inches='tight')
                        
                else:

                    for m in range(0,len(orbitals)):
                        M = orbitals[m]
                        for i in range(0,len(num)):
                            n = num[i]-1
                            N = num[i]
                            k1s = a[2,32*Nc*M+n*Nc:32*Nc*M+n*Nc+Nc]  ##first K
                            k2s = a[3,32*Nc*M+n*Nc:32*Nc*M+n*Nc+Nc]  ##second K
                            print(a[0:2,32*Nc*M+n*Nc:32*Nc*M+n*Nc+Nc])
                            Vec = list(a[5,32*Nc*M+n*Nc:32*Nc*M+n*Nc+Nc])

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
                          
                            fig,ax = plt.subplots()
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

                            for i in range(len(index)):
                                if M == i:
                                    ax.set_title(str(index[i])+'_['+str(N)+']',fontsize=40)
                                    plt.savefig('orbital_'+str(N)+'_Nc'+str(Nc)+'_U'+str(U)+'_d'+str(d)+'_T'+str(T)+'_'+str(index[i])+'.png',bbox_inches='tight')

######################################################################################################################

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
                            ax.set_title('$Anti-bondingd_{x^2-y^2}$_' + str(N), fontsize=40)
                            plt.savefig('band_' + str(N) + '_Nc' + str(Nc) + '_U' + str(U) + '_d' + str(d) + '_T' + str(T) + '_antibondingx.pdf',bbox_inches='tight')
                        elif M == 2:
                            ax.set_title('$bondingd_{z^2}$_' + str(N), fontsize=40)
                            plt.savefig('band_' + str(N) + '_Nc' + str(Nc) + '_U' + str(U) + '_d' + str(d) + '_T' + str(T) + '_bondingz.pdf',bbox_inches='tight')

                        elif M == 3:
                            ax.set_title('$Anti-bondingd_{z^2}$_' + str(N), fontsize=40)
                            plt.savefig('band_' + str(N) + '_Nc' + str(Nc) + '_U' + str(U) + '_d' + str(d) + '_T' + str(T) + '_antibondingz.pdf',bbox_inches='tight')
