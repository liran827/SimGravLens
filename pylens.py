import numpy as np
from matplotlib.pyplot import *
import matplotlib.pyplot as plt
from math import *
import halo as ho
import cosmo

from ctypes import *
lenslib=CDLL("./clens/liblens.so")
construct_lens_plane=lenslib.construct_lens_plane
mapping_source=lenslib.mapping_source


from numpy.ctypeslib import ndpointer
#the input is the smoothed density file from sph smoothing process
#all units includes h



def read_sdens(file_base, boxsize, nnn,mp):
    x1 = np.fromfile(file_base+'_posx1.bin',dtype=np.float32)
    x1 = x1.reshape((nnn,nnn))#/apr*Da(0.2)
    x2 = np.fromfile(file_base+'_posx2.bin',dtype=np.float32)
    x2 = x2.reshape((nnn,nnn))#/apr*Da(0.2)

    sdens_sph = np.fromfile(file_base+'_sdens.bin',dtype=np.float32)
    sdens_sph = sdens_sph.reshape((nnn,nnn))


    return [x1,x2,sdens_sph*mp]

def SIS(r,vdisp,rlim):
    if(r>rlim):
        return 0
    return pow(vdisp,2.)/2./cosmo.G/r/1000.

def SIS_array(r,vdisp,rlim):
    res=np.zeros_like(r)
    ind,=np.where(r<rlim)

    res[ind]=np.power(vdisp,2.)/2./cosmo.G/r[ind]/1000.
    return res

def SIS_mesh(vdisp,xi,yi,nnn,cellsize):
    #xi,yi the position of the lens, in kpc at its lens plane
    #return density in unit of 1e10*Msun/kpc^2, with or without h factor
    hN=nnn/2.
    cx1=hN-0.5+xi#//center of big halo
    cy1=hN-0.5+yi

    rx=(np.arange(nnn)-cx1)*cellsize
    ry=(np.arange(nnn)-cy1)*cellsize

    xmesh,ymesh=np.meshgrid(rx,ry)
    rr=np.sqrt(np.power(xmesh,2.)+np.power(ymesh,2.))
    rr,=np.array(rr.reshape((1,nnn*nnn),order='C'),dtype='float64')
    sis_dens= SIS_array(rr,vdisp,200.)/1e10;
    return sis_dens

def NFW_mesh(mass,z,xi,yi,nnn,cellsize):
    #xi,yi the position of the lens, in kpc at its lens plane
    #cellsize is the physical size of the grid,kpc
    #nnn is the meshgrid

    hN=nnn/2.
    cx1=hN-0.5+xi#//center of big halo
    cy1=hN-0.5+yi

    rx=(np.arange(nnn)-cx1)*cellsize
    ry=(np.arange(nnn)-cy1)*cellsize

    xmesh,ymesh=np.meshgrid(rx,ry)
    rr=np.sqrt(np.power(xmesh,2.)+np.power(ymesh,2.))
    rr,=np.array(rr.reshape((1,nnn*nnn),order='C'),dtype='float64')/1000. # in Mpc

    c200,rs,r200,rhos=ho.nfw_para(mass,z) # in unit of Mpc

    return ho.dens_tnfw_array(rr,mass,rs,r200)/1e4





def DFL_angle(sigma2d, z_lens, z_source, lensNumberSide,lensSize):
    #output dpdx dpdy

    dpdx=np.zeros(lensNumberSide*lensNumberSide)
    dpdy=np.zeros(lensNumberSide*lensNumberSide)

    lensAngularDistance = cosmo.Da(0.,z_lens)*1000.#unit kpc h^-1
    sourceAngularDistance = cosmo.Da(0.,z_source)*1000.#unit kpc h^-1
    source_lensAngularDistance = sourceAngularDistance - lensAngularDistance*(1+z_lens)/(1+z_source)

    print lensAngularDistance, sourceAngularDistance, source_lensAngularDistance
    construct_lens_plane(sigma2d.ctypes.data_as(POINTER(c_double)),c_double(lensAngularDistance),c_double(source_lensAngularDistance),c_double(sourceAngularDistance),c_double(lensSize),c_int(lensNumberSide), dpdx.ctypes.data_as(POINTER(c_double)), dpdy.ctypes.data_as(POINTER(c_double)) )

    return [dpdx,dpdy]


def lens_mag(xmesh,ymesh,dpdx,dpdy,nnn):
    #input deflection angle map dpdx, dpdy (nnnxnnn 1d array), in rad
    #xmesh,ymesh are coordinates in unit of arcsec


    ph_1=dpdx.reshape(nnn,nnn)
    ph_2=dpdy.reshape(nnn,nnn)



    ph11=np.zeros((nnn,nnn))
    ph22=np.zeros((nnn,nnn))
    ph12=np.zeros((nnn,nnn))

    ph11[:-1,:-1]=(np.diff(ph_1,axis=1)/np.diff(xmesh,axis=1))[:-1,:]#unit arcsec
    ph22[:-1,:-1]=(np.diff(ph_2,axis=0)/np.diff(ymesh,axis=0))[:,:-1]
    ph12[:-1,:-1]=(np.diff(ph_1,axis=0))[:,:-1]/np.diff(ymesh,axis=0)[:,:-1]


    #kappa=sigma2d.reshape(nnn,nnn)/ (sigma_crit / 1e4)#xxx unit

    kappa=0.5*(ph11+ph22)
    gm1=0.5*(ph11-ph22)
    gm2=ph12

    gm_mod=np.sqrt( np.power(gm1,2.)+np.power(gm2,2.) )

    mu=1./( np.power(1-kappa,2.)- np.power( gm_mod , 2.)    )

    sign=np.zeros_like(mu)
    sign[0:-1,:]=(mu[0:-1,:]*mu[1:,:])
    ind1=np.where(sign<0)
    sign[:,0:-1]=(mu[:,0:-1]*mu[:,1:])
    ind2=np.where(sign<0)

    ind=[np.append(ind1[0],ind2[0]),np.append(ind1[1],ind2[1])]



    return mu,ind
