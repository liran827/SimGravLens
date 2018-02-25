import numpy as np
import random as rand
import sys
import time
from scipy.interpolate import interp1d
from scipy.integrate import quad
from scipy.optimize import  brentq
import cosmo
print cosmo.H0
from scipy.optimize import fsolve
from math import *
import halo as ho

rad2min=180./np.pi*60.
rad2sec=180./np.pi*3600.

def MER_SIS(theta,zl,zs):

    Dl=cosmo.Da(0.,zl)
    Ds=cosmo.Da(0.,zs)
    Dls=cosmo.Da(zl,zs)


    Dfact=Dls/Ds/Dl

    M=(theta/rad2sec)**2.*cosmo.VC**2./4./cosmo.G/Dfact
    R=Dl*theta/rad2sec
    Vdisp=np.sqrt(M/R/np.pi*cosmo.G)
    return [M,Vdisp,R]

def ERing_SIS(sigma,zl,zs):
    Dl=cosmo.Da(0.,zl)
    Ds=cosmo.Da(0.,zs)
    Dls=cosmo.Da(zl,zs)


    Dfact=Dls/Ds/Dl
    theta=4.*np.pi*sigma**2./cosmo.VC**2.*Dls/Ds
    RE=theta*Dl

    return [RE,theta]

def ERing_PM(mass,zl,zs):
    Dl=cosmo.Da(0.,zl)
    Ds=cosmo.Da(0.,zs)
    Dls=cosmo.Da(zl,zs)
    Dfact=Dls/Ds/Dl


    theta=np.sqrt(4.*mass*cosmo.G*Dfact)/cosmo.VC
    RE=Dl*theta

    return [RE,theta]


def ERing_ENFW_generator():
    #The function is to calculate ERing-M table using Schaller et al. 2014
    #out put a M-Ring size distribution
    table=np.loadtxt("table.schaller")
    table=table[5:,:]

    rs1=table[:,2]
    dt1=table[:,4]
    
    rs2=table[:,5]
    dt2=table[:,6]

    R200=table[:,1]

    zl=1.0
    zs=2.0

    fout="ENFW_ringsize_table_"+"zl"+str(zl)+"_zs"+str(zs)
    for i in range(rs1.shape[0]):
        M200=Mproj_SNFW(R200[i],dt1[i],rs1[i],dt2[i],rs2[i])
        RE,theta=FindRing_SNFW(dt1[i],rs1[i],dt2[i],rs2[i])
        print >>fout,"%f %f %f"%(M200,RE,theta)
    return 0

def tsol_nfw_ring(r,para):
    dt,rs,rhoc,Dfact=para
    #print r,Mproj_nfw(r,dt,rs,rhoc)
    return Dfact*4.*cosmo.G*Mproj_nfw(r,dt,rs,rhoc)/cosmo.VC**2.-r*r

def ERing_NFW(mass,zl,zs): 
    Dl=cosmo.Da(0.,zl)
    Ds=cosmo.Da(0.,zs)
    Dls=cosmo.Da(zl,zs)
    Dfact=Dl*Dls/Ds # not the same Dfact as those in above function
    [c200,rs,r200,rhos]=ho.nfw_para(mass,zl)
    dt=1.
    para=[dt,rs,rhos,Dfact]
    RE=fsolve(tsol_nfw_ring,1e-6*r200,args=para)
    theta=RE/Dl
    return [RE,theta]


#NFW profile properties from Wright 2000

def Mproj_nfw(r,dt,rs,rhoc):
    #projected mass with the radius r for a nfw halo
    x=r/rs
    SS=np.pi*r**2.
    if(x<1):
        return SS* 4./x/x*rs*dt*rhoc * ( \
            2./sqrt(1-x**2.)*np.arctanh( sqrt( (1-x)/(1+x)) ) +log(x/2.)  )
    elif(x==1):
        return SS*4.*rs*dt*rhoc*(1+log(0.5))
    else:
        return SS*4./x/x*rs*dt*rhoc * ( \
            2./sqrt(x**2.-1)*np.arctan( sqrt( (x-1)/(1+x)) ) + log(x/2.)  )

def Mproj_SNFW(r,para): 
    dt1,rs1,dt2,rs2=para
    rhoc=cosmo.rho_crit(0)         
    return Mproj_nfw(r,dt1,rs1,rhoc)+Mproj_nfw(r,dt2,rs2,rhoc)




