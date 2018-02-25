# The programe is to sample a PIEMD halo + severl PIEMD/NFW subhaloes to a mesh.

import matplotlib.pyplot as plt
from math import *
import numpy as np
import random as rand
import sys
import time
from scipy.interpolate import interp1d
from scipy.integrate import quad
from scipy.optimize import  brentq
import cosmo

def mc(x):
    return log(1+x)-x/(1+x) 

def dens_tnfw(r,mass,rs,rt):
    ct=rt/rs
    mi=ct*ct/(2.* np.power(ct*ct+1,3)*(1+ct)*2.*ct*ct) \
            *( (ct*ct+1)*ct * (ct*(ct+1)-ct*ct*(ct-1)*(2+3*ct)-2*np.power(ct,4.)) \
            +ct*(ct+1)* 2.*ct*ct*(2*(3*ct*ct-1)*atan(1.) \
            +ct*(ct*ct-3)*np.log(ct*ct*np.power(1+ct,2)/(2*ct*ct)) ) )
    rhos=mass/(4.*np.pi*np.power(rs,3.)*mi)/1e12

    tau=rt/rs
    x=r/rs

    if(x<1):
        Fx=1/np.sqrt(1-x*x)*np.arctanh(np.sqrt(1-x*x))
    else:
        Fx=1/np.sqrt(x*x-1)*np.arctan(np.sqrt(x*x-1))
    Lx=np.log(x/(tau+np.sqrt(tau*tau+x*x)))

    sigma_r=4*rhos*rs* np.power(tau,4.)/(4.*np.power(tau*tau+1,3)) * (2*(tau*tau+1)/(x*x-1)*(1-Fx) + 8*Fx + (np.power(tau,4.)-1)/tau/tau/(tau*tau+x*x)  \
            -np.pi*(4*(tau*tau+x*x) + tau*tau +1)/np.power(tau*tau+x*x,1.5) \
            +( tau*tau*(np.power(tau,4.)-1) + (tau*tau+x*x)*(3*np.power(tau,4.)-6*tau*tau-1) )*Lx  \
            /np.power(tau,3.)/np.power(tau*tau+x*x,1.5) )


    return sigma_r

def dens_tnfw_array(r,mass,rs,rt):
    ct=rt/rs
    mi=ct*ct/(2.* np.power(ct*ct+1,3)*(1+ct)*2.*ct*ct) \
            *( (ct*ct+1)*ct * (ct*(ct+1)-ct*ct*(ct-1)*(2+3*ct)-2*np.power(ct,4.)) \
            +ct*(ct+1)* 2.*ct*ct*(2*(3*ct*ct-1)*atan(1.) \
            +ct*(ct*ct-3)*np.log(ct*ct*np.power(1+ct,2)/(2*ct*ct)) ) )
    rhos=mass/(4.*np.pi*np.power(rs,3.)*mi)/1e12

    tau=rt/rs
    x=r/rs
    
    Fx=np.zeros_like(x)
    ind,=np.where(x<1)
    Fx[ind]=1/np.sqrt(1-x[ind]*x[ind])*np.arctanh(np.sqrt(1-x[ind]*x[ind]))
    ind,=np.where(x>=1)
    Fx[ind]=1/np.sqrt(x[ind]*x[ind]-1)*np.arctan(np.sqrt(x[ind]*x[ind]-1))

    Lx=np.log(x/(tau+np.sqrt(tau*tau+x*x)))

    sigma_r=4*rhos*rs* np.power(tau,4.)/(4.*np.power(tau*tau+1,3)) * (2*(tau*tau+1)/(x*x-1)*(1-Fx) + 8*Fx + (np.power(tau,4.)-1)/tau/tau/(tau*tau+x*x)  \
            -np.pi*(4*(tau*tau+x*x) + tau*tau +1)/np.power(tau*tau+x*x,1.5) \
            +( tau*tau*(np.power(tau,4.)-1) + (tau*tau+x*x)*(3*np.power(tau,4.)-6*tau*tau-1) )*Lx  \
            /np.power(tau,3.)/np.power(tau*tau+x*x,1.5) )


    return sigma_r


def dens_tnfw_in(r,mass,rs,rt):
    ct=rt/rs
    mi=ct*ct/(2.* np.power(ct*ct+1,3)*(1+ct)*2.*ct*ct) \
            *( (ct*ct+1)*ct * (ct*(ct+1)-ct*ct*(ct-1)*(2+3*ct)-2*np.power(ct,4.)) \
            +ct*(ct+1)* 2.*ct*ct*(2*(3*ct*ct-1)*atan(1.) \
            +ct*(ct*ct-3)*np.log(ct*ct*np.power(1+ct,2)/(2*ct*ct)) ) )
    rhos=mass/(4.*np.pi*np.power(rs,3.)*mi)/1e12

    tau=rt/rs
    x=r/rs

    if(x<1):
        Fx=1/np.sqrt(1-x*x)*np.arctanh(np.sqrt(1-x*x))
    else:
        Fx=1/np.sqrt(x*x-1)*np.arctan(np.sqrt(x*x-1))
    Lx=np.log(x/(tau+np.sqrt(tau*tau+x*x)))

    sigma_in=4*rhos*rs*np.power(tau,4.)/2./np.power(tau*tau+1,3)/x/x * \
          ( 2*(tau*tau+4*x*x-3)*Fx + 1./tau*(np.pi*(3*tau*tau-1)+2*tau*(tau*tau-3)*np.log(tau)) \
            + 1./np.power(tau,3.)/np.sqrt(tau*tau+x*x)*(-tau*tau*tau*np.pi*(4*x*x+3*tau*tau-1) + (2*np.power(tau,4)*(tau*tau-3.)+np.power(x,2.)*(3*np.power(tau,4)-6*tau*tau-1))*Lx))  
    return sigma_in




def mass2concentration(m,z):
    #return concentration defined as r200_crit/rs
    return 4.67/(1. + z) * pow(m/1.e14 , -0.11)

def nfw_para(m200,z):
    c200=mass2concentration(m200,z)
    r200=np.power(m200/(4.*np.pi/3.*200.*cosmo.rho_crit(z)),1./3.)
    rs=r200/c200
    rhos=cosmo.rho_crit(z)*200./3.* pow(c200 , 3.)/(log(1. + c200)- c200/(1.+c200) );

    return [c200,rs,r200,rhos]

def nfw_para2(m200,c200,z):
    r200=np.power(m200/(4.*np.pi/3.*200.*cosmo.rho_crit(z)),1./3.)
    rs=r200/c200
    rhos=cosmo.rho_crit(z)*200./3.* pow(c200 , 3.)/(log(1. + c200)- c200/(1.+c200) );

    return [c200,rs,r200,rhos]

def profile_tnfw_2d(r,rs,rt):
    tau=rt/rs
    x=r/rs
    if(x<1):
        Fx=1/np.sqrt(1-x*x)*np.arctanh(np.sqrt(1-x*x))
    else:
        Fx=1/np.sqrt(x*x-1)*np.arctan(np.sqrt(x*x-1))
    Lx=np.log(x/(tau+np.sqrt(tau*tau+x*x)))
    sigma_r=1./(4.*np.power(tau*tau+1,3)) * (2*(tau*tau+1)/(x*x-1)*(1-Fx) + 8*Fx + (np.power(tau,4.)-1)/tau/tau/(tau*tau+x*x)  \
            -np.pi*(4*(tau*tau+x*x) + tau*tau +1)/np.power(tau*tau+x*x,1.5) \
            +( tau*tau*(np.power(tau,4.)-1) + (tau*tau+x*x)*(3*np.power(tau,4.)-6*tau*tau-1) )*Lx  \
            /np.power(tau,3.)/np.power(tau*tau+x*x,1.5) )
    return sigma_r

def M2disp(M,R):
    #mass with a projected radius R, with SIS model
    return sqrt(M*cosmo.G/np.pi/R)

def dens_sis(R,sigma):
    return sigma**2./2./cosmo.G/R/1e12





    
