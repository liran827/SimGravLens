import numpy as np
import scipy.integrate as integrate
import math 


VC=300000.
H0=67.8
G=4.302e-9
Omega_L0=0.728
Omega_M0=0.272

print "H0=",H0
print "Do you need to set H0=100?"

def set_para(h0=67.8,OL=0.728,OM=0.272):
    global H0, Omega_L0, Omega_M0
    H0=h0
    Omega_L0=OL
    Omega_M0=OM
    print "H0=",H0
    print "Omega_L0=",Omega_L0
    print "Omega_M0=",Omega_M0
    return 0

def reset_para():
    global H0, Omega_L0, Omega_M0,VC, G
    VC=300000.
    H0=67.8
    G=4.302e-9
    Omega_L0=0.728
    Omega_M0=0.272

    print "H0=",H0
    print "Do you need to set H0=100?"


    return 0



def rho_crit(z):
    fac=Omega_L0 + Omega_M0* np.power(1.0+z , 3.)
    xH=H0*np.sqrt(fac)
    rho=3.0*xH*xH/(8.0 * math.pi * G)
    return rho

# Physical coordinate, R500
def R500(m,z):
    r500=np.power( m/(4.* math.pi /3. * rho_crit(z)*500.) , 1./3. )
    return r500



#Angular Diameter distance
def IntD(a):
    return VC/H0 * np.power( ( Omega_M0*a + Omega_L0*np.power(a , 4.) ) , -0.5)

def Da(z1,z2):
    a1=1./(1.+z1)
    a2=1./(1.+z2)
    DD=integrate.quad(IntD, a2, a1)
    return DD[0]*a2

def Dc(z1,z2):
    a1=1./(1.+z1)
    a2=1./(1.+z2)
    DD=integrate.quad(IntD, a2, a1)
    return DD[0]

def dDcom_dz(z):
    a=1./(1+z)
    return VC/H0 * np.power( ( Omega_M0*a + Omega_L0*np.power(a , 4.) ) , -0.5)*a*a




#Angular distance between to position

def angdist(ra1,dec1,ra2,dec2,degree=1):
    conv=180./np.pi
    if(degree > 0):
	ra1=ra1/conv
	ra2=ra2/conv
	dec1=dec1/conv
	dec2=dec2/conv

    res=np.arccos(np.sin(dec1)*np.sin(dec2) + np.cos(dec1)*np.cos(dec2)*np.cos(ra1-ra2) )

    if(degree > 0):
        res=res*conv

    return np.abs(res)


def sigma_crit(zl,zs):
    VC= 300000.
    G=4.2994e-9
    return np.power(VC,2.)/(4.*np.pi*G)*Da(0.,zs)/Da(0.,zl)/Da(zl,zs)/1e12



def sigma_crit_comv(Dl,Ds,zl):
    VC= 300000.
    G=4.2994e-9
    Dls=Ds-Dl
    return np.power(VC,2.)/(4.*np.pi*G)*Ds/Dl/Dls/(1+zl)/1e12


def sigma_crit_phy(zl , zs, Dl=0, Ds=0):
    #Dl,Ds,Dls are all comoving distance, returns sigma_crit in comoving coordinate
    VC= 300000.
    G=4.2994e-9
    a_l=1./(1.+zl)
    a_s=1./(1.+zs)
    if(Dl==0):
        Dl=Da(0.,zl)
        Ds=Da(0.,zs)
    Dls=Ds-Dl/a_l*a_s
    return np.power(VC,2.)/(4.*np.pi*G)*Ds/Dl/Dls/1e12

    

def resample(Nsub):
    ind=np.floor(np.random.rand(Nsub)*Nsub)
    ind=ind.astype(int)
    weight=np.zeros(Nsub)
    for i in range(Nsub):
        weight[ind[i]]+=1
    return weight

    return VC/H0 * np.power( ( Omega_M0*a + Omega_L0*np.power(a , 4.) ) , -0.5)*a*a
