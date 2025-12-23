import matplotlib
matplotlib.use('Agg')

import subprocess
import os
import sys
import time
import shutil

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import math

from pylab import *
import scipy
from scipy import *
from scipy import interpolate
from numpy import *
from sympy import symbols, Matrix, init_printing, simplify,lambdify
from scipy.linalg import eigh

import re
import mmap
from linecache import getline

# use the common py analysis toolbox valid for all DQMC projects

import read_data as data
import domains as domain
import utility as util

Ms = ['o','s','^','v','p','h','D','8','<','>','H','o','s','^','v','p','h','D','8','<','>','H', \
      'o','s','^','v','p','h','D','8','<','>','H','o','s','^','v','p','h','D','8','<','>','H']
#colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
colors = ['C0', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9', 'C0',\
          'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8',\
          'C0', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9', 'C0',\
          'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8']

######################################################
def Hamitonlian(kx,ky,chemicalpotential):
    m = chemicalpotential
    Hx  = ex + 2.*t11x*(cos(kx)+cos(ky)) + 4.*t11xy*cos(kx)*cos(ky) + 2.*t11xx*(cos(2*kx)+cos(2*ky)) - m
    Hz  = ez + 2.*t22x*(cos(kx)+cos(ky)) + 4.*t22xy*cos(kx)*cos(ky) + 2.*t22xx*(cos(2*kx)+cos(2*ky)) - m
    V   = 2.*t12x*(cos(kx)-cos(ky)) + 2.*t12xx*(cos(2*kx)-cos(2*ky))

    Hxp = s110 + 2.*s11x*(cos(kx)+cos(ky)) + 4. * s11xy*cos(kx)*cos(ky) + 2.*s11xx*(cos(2*kx)+cos(2*ky))
    Hzp = s220 + 2.*s22x*(cos(kx)+cos(ky)) + 4. * s22xy*cos(kx)*cos(ky) + 2.*s22xx*(cos(2*kx)+cos(2*ky))
    Vp  = 2.*s12x*(cos(kx)-cos(ky)) + 2.*s12xx*(cos(2.*kx)-cos(2.*ky))
    
    H = np.array([[Hx, V, Hxp, Vp],
                  [V, Hz, Vp, Hzp],
                  [Hxp, Vp, Hx, V],
                  [Vp, Hzp, V, Hz]], dtype=complex)
    return H

def Matrix(kx,ky,chemicalpotential):
    m = chemicalpotential
    H = Hamitonlian(kx,ky,m)
    eigvals, eigvecs = np.linalg.eigh(H)  # eigvecs[:, i] to eigvals[i]
    mat = eigvecs
    matinv = np.linalg.inv(mat)
    return matinv, mat

def findmuforNonband(density,number):
    n = density
    N = number
    a1 = -4; b1 = 4
    m = (a1 + b1) / 2.0
    d = 0.0;d1=0.0;d2=0.0
    X = []
    while n >= 0.0:
        for kx in arange(-pi,pi,pi/N):
            for ky in arange(-pi,pi,pi/N):
                H = Hamitonlian(kx,ky,m)
                
                eigvals = np.linalg.eigvalsh(H)
                
                x = sum(1 for e in eigvals if e <= 0)
                d += x * 2  

        n1 = d / (2 * N * 2 * N )
        print('n1:',n1)
        print('m:',m)
        if abs(n1 - n) < 0.0005:
            break

        if n1 > n:
            b1 = m
            m = (a1 + b1) / 2.0
            d = 0.0
        elif n1 < n:
            a1 = m
            m = (a1 + b1) / 2.0
            d = 0.0

    print('density =', n1)
    print('mu =', m)
    return n1,m,N


def noninteractingband(chemicalpotential,number):
    m = chemicalpotential
    N = number
    ye0=[];ye1=[];ye2=[];ye3=[]
    ye00=[];ye01=[];ye02=[]
    ye10=[];ye11=[];ye12=[]
    ye20=[];ye21=[];ye22=[]
    ye30=[];ye31=[];ye32=[]
    for kx in arange(0,pi+pi/N,pi/N):
        ky = 0
        H = Hamitonlian(kx,ky,m)

        eigvals = np.linalg.eigvalsh(H)
        
        ye00.append(eigvals[0]);ye10.append(eigvals[1]);ye20.append(eigvals[2]);ye30.append(eigvals[3])

    for ky in arange(pi/N,pi+pi/N,pi/N):
        kx = pi
        H = Hamitonlian(kx,ky,m)

        eigvals = np.linalg.eigvalsh(H)
        
        ye01.append(eigvals[0]);ye11.append(eigvals[1]);ye21.append(eigvals[2]);ye31.append(eigvals[3])

    for kx in arange(-pi+pi/N,0,pi/N):
        ky = kx
        H = Hamitonlian(kx,ky,m)

        eigvals = np.linalg.eigvalsh(H)
        
        ye02.append(eigvals[0]);ye12.append(eigvals[1]);ye22.append(eigvals[2]);ye32.append(eigvals[3])

    ye0 = ye00+ye01+ye02;ye1 = ye10+ye11+ye12;ye2 = ye20+ye21+ye22;ye3 = ye30+ye31+ye32;
    return ye0,ye1,ye2,ye3

######################################################################3
pi = math.pi
Nc = 8
Ts = [0.08]
ds = [2.4]

#model = 'tensile'
#model = 'compress'
model = 'compress_full'
#model = 'maier'

if model == 'compress':
    ex = 1.0493
    ez = 0.3918
    t11x = -0.4748
    t11xy = 0.0775
    t11xx = 0.0
    t22x =  -0.0779
    t22xy = -0.0147
    t22xx = 0.0
    t12x = 0.2051
    t12xx = 0.0
    s110 = 0.0083
    s11x = 0.0
    s11xy = 0.0
    s11xx = 0.0
    s220 =  -0.6174
    s22x =  0.0
    s22xy = 0.0
    s22xx = 0.0
    s12x =  -0.0277
    s12xx = 0.0

if model == 'compress_full':
    ex = 1.0953
    ez = 0.4378
    t11x = -0.4748
    t11xy = 0.0775
    t11xx = -0.0594
    t22x =  -0.0779
    t22xy = -0.0147
    t22xx = -0.0144
    t12x = 0.2051
    t12xx = 0.0264
    s110 = 0.0083
    s11x = -0.001
    s11xy = 0.0022
    s11xx = -0.0032
    s220 =  -0.6174
    s22x =  0.0122
    s22xy = 0.0068
    s22xx = 0.0022
    s12x =  -0.0277
    s12xx = 0.0002

if model == 'tensile':
    ex = 0.5071
    ez = 0.376
    t11x = -0.4101
    t11xy = 0.0
    t11xx = 0.0
    t22x =  -0.1128
    t22xy = 0.0
    t22xx = 0.0
    t12x = 0.2251
    t12xx = 0.0
    s110 = 0.0
    s11x = 0.0
    s11xy = 0.0
    s11xx = 0.0
    s220 =  -0.5772
    s22x =  0.0
    s22xy = 0.0
    s22xx = 0.0
    s12x =  0.0
    s12xx = 0.0

if model == 'maier':
    ex = 0.506
    ez = 0.0
    t11x = -0.515
    t11xy = 0.0
    t11xx = 0.0
    t22x =  -0.11
    t22xy = 0.0
    t22xx = 0.0
    t12x = 0.243
    t12xx = 0.0
    s110 = 0.0
    s11x = 0.0
    s11xy = 0.0
    s11xx = 0.0
    s220 =  -0.666
    s22x =  0.0
    s22xy = 0.0
    s22xx = 0.0
    s12x =  0.0
    s12xx = 0.0

print(ex)

########################################################
nvals = [];muvals =[]
SigmaRekzx0 = []; SigmaRekzxpi = []
SigmaImkzx0 = []; SigmaImkzxpi = []
SigmaRekzz0 = []; SigmaRekzzpi = []
SigmaImkzz0 = []; SigmaImkzzpi = []
Zkzx0 = []; Zkzxpi = []
Zkzz0 = []; Zkzzpi = []

for iid in range(0,len(ds)):
    d = ds[iid]
    nvals = [];muvals =[]
    SigmaRekx0 = []; SigmaRekxpi = []
    SigmaRekz0 = []; SigmaRekzpi = []
    Zkx0 = []; Zkxpi = []
    Zkz0 = []; Zkzpi = []
    n,m,N = findmuforNonband(d,100)
    for iT in range(len(Ts)):
        T = Ts[iT]
        # To record Sigma(w=0) at kz=0 , 2*pi/3 and 4*pi/3 vs T
        fname0 = './data/Nc' + str(Nc) + '/'+ str(model) + '/Sigma_vs_n'+str(d)+'_Re.txt'
        fname1 = './data/Nc' + str(Nc) + '/'+ str(model) + '/Sigma_vs_n'+str(d)+'_Zk.txt'

        util.del_existing_file(fname0);util.del_existing_file(fname1)

        if model == 'compress_full':
            dataname = '../Nc'+str(Nc)+'/compress/Full_TB/d'+str(d)+'/T='+str(T)+'/dca_sp.hdf5'
        else:
            dataname = '../Nc'+str(Nc)+'/'+str(model)+'/d'+str(d)+'/T='+str(T)+'/dca_sp.hdf5'

        if os.path.isfile(dataname):

            print(dataname)
            Ks = []
            Rvecs, Kvecs, qchannel, iQ = domain.Get_K(dataname)
            Kxs = [0, pi, pi]
            Kys = [0, 0, pi]
            iK00   = domain.K_2_iK(Kxs[0], Kys[0], Kvecs)
            iKpi0  = domain.K_2_iK(Kxs[1], Kys[1], Kvecs)
            #iKpi2pi2 = domain.K_2_iK(pi/2, pi/2, Kvecs)
            iKpipi = domain.K_2_iK(Kxs[2], Kys[2], Kvecs)
            #print(iK00, iKpipi, iKpi0)
            Ks.append(iK00);Ks.append(iKpi0);Ks.append(iKpipi)
            print(Ks)
            Nwsp, wn = domain.Get_w_sp(dataname)
            Nwsp_2 = int(Nwsp/2)

            sigmaRe, sigmaIm = data.Get_Sigma(dataname)
            mu = data.Get_mu(dataname)
            #print(sigmaIm[Nwsp_2:Nwsp_2+3,iKpi0,0,2,0,2])
            #print("G.shape=", sigmaIm.shape)

            print(wn[Nwsp_2]/pi)
            for i in range(len(Ks)):
                Uinv,U = Matrix(Kxs[i],Kys[i],m)
                #print(U)
                SigmaReM = np.matmul(Uinv,np.matmul(sigmaRe[Nwsp_2,Ks[i],0,:,0,:],U))
                SigmaImM = np.matmul(Uinv,np.matmul(sigmaIm[Nwsp_2,Ks[i],0,:,0,:],U))

                Zz0  = 1/(1-SigmaImM[0][0]/wn[Nwsp_2])
                Zzpi = 1/(1-SigmaImM[1][1]/wn[Nwsp_2])
                Zx0  = 1/(1-SigmaImM[2][2]/wn[Nwsp_2])
                Zxpi = 1/(1-SigmaImM[3][3]/wn[Nwsp_2])

                SigmaRekz0.append(SigmaReM[0][0]); SigmaRekzpi.append(SigmaReM[1][1])
                SigmaRekx0.append(SigmaReM[2][2]); SigmaRekxpi.append(SigmaReM[3][3])

                Zkz0.append(Zz0);Zkzpi.append(Zzpi)
                Zkx0.append(Zx0);Zkxpi.append(Zxpi)
                nvals.append(d);muvals.append(mu)

    if len(SigmaRekz0)>0 :
        util.write_data_5cols(fname0, nvals, SigmaRekz0, SigmaRekzpi, SigmaRekx0, SigmaRekxpi)
        util.write_data_6cols(fname1, nvals, Zkz0, Zkzpi, Zkx0, Zkxpi, muvals)

##################################################

for iid in range(len(ds)):
    d = ds[iid]
    dname0 = './data/Nc' + str(Nc) + '/'+ str(model) + '/Sigma_vs_n'+str(d)+'_Re.txt'
    dname1 = './data/Nc' + str(Nc) + '/'+ str(model) + '/Sigma_vs_n'+str(d)+'_Zk.txt'

    a = loadtxt(dname0,unpack=True)
    b = loadtxt(dname1,unpack=True)

    #print(a)
    ns = a[0,:]
    #print(ns)
    mus = b[5,:]

    sigmaRez0  = a[1,:];Zz0  = b[1,:]
    sigmaRezpi = a[2,:];Zzpi = b[2,:]
    sigmaRex0  = a[3,:];Zx0  = b[3,:]
    sigmaRexpi = a[4,:];Zxpi = b[4,:]
    #print(mus[1])

    Kxs = [0,pi,pi]
    Kys = [0,0, pi]
    Y0 = [];Y1 = [];Y2 = [];Y3 = []

    for i in range(len(ns)):
        kx = Kxs[i];ky=Kys[i]
        H = Hamitonlian(kx,ky,mus[i])

        E = np.linalg.eigvalsh(H)
        #print(E)

        Y0.append((E[0]+sigmaRez0[i])*Zz0[i]);Y1.append((E[1]+sigmaRezpi[i])*Zzpi[i])
        Y2.append((E[2]+sigmaRex0[i])*Zx0[i]);Y3.append((E[3]+sigmaRexpi[i])*Zxpi[i])

    Y0.append(Y0[0]);Y1.append(Y1[0]);Y2.append(Y2[0]);Y3.append(Y3[0])

    X =[0,pi,2*pi,3*pi]
    plt.plot(X,Y0,color=colors[iid],marker=Ms[iid],markersize=8)#,label=r'$BondingZ,n$'+str(d))
    plt.plot(X,Y1,color=colors[iid+1],marker=Ms[iid],markersize=8)#,label=r'$BondingZ,n$'+str(d))
    plt.plot(X,Y2,color=colors[iid+2],marker=Ms[iid],markersize=8)#,label=r'$BondingZ,n$'+str(d))
    plt.plot(X,Y3,color=colors[iid+3],marker=Ms[iid],markersize=8)#,label=r'$BondingZ,n$'+str(d))

################# non interacting #######################################################


#n,m,N = findmuforNonband(d,100)
xe = np.linspace(0,3*pi,3*N)
ye0,ye1,ye2,ye3 = noninteractingband(m,N)
plt.plot(xe,ye0,linestyle='--',color=colors[iid])#,label='E0')
plt.plot(xe,ye1,linestyle='--',color=colors[iid+1])#,label='E1')
plt.plot(xe,ye2,linestyle='--',color=colors[iid+2])#,label='E2')
plt.plot(xe,ye3,linestyle='--',color=colors[iid+3])#,label='E3')

plt.plot([0,3*pi],[0,0],linestyle='--',color='gray')
    
#plt.text(1.5,-0.3,'$T=0.125$eV',fontsize=20)
plt.title(str(model)+'_T'+str(T)+'_Nc'+str(Nc)+'_d'+str(d),fontsize=18)
#plt.legend(loc='best',fontsize=18)
plt.ylim(-0.8,3.5)
X =[0,pi,2*pi,3*pi]

plt.xlim(0,3*pi)
plt.ylabel('$Z\epsilon^*(k)$',fontsize=28)
plt.grid(axis='x')
plt.xticks(X,['$\Gamma$','$X$','$M$','$\Gamma$'],fontsize=20)
plt.yticks(fontsize=20)
plt.savefig(str(model)+'_Nc'+str(Nc)+'_T'+str(T)+'_d'+str(d)+'.pdf',bbox_inches='tight')
