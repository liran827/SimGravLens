#resimulate analytical model
import numpy as np
import matplotlib.pyplot as plt
import cosmo
cosmo.set_para(h0=100.) #set H0=100
import lensinglib as ll
import time

try:
    import astropy.io.fits as pyfits
except:
    import pyfits
#filenumber
Num=35

#source parameters
n_source=2 #number of source galaxy
spar=np.zeros((2,7))
#spar_data=np.loadtxt("sersic_source_param.dat")
#spar=spar_data[0:3,:]

# I0, xc, yc, r0, phi, ellp, ns
spar[0,0]=30. #I0
spar[0,1]=-0.077527 #xc
spar[0,2]=-0.0997  #yc
spar[0,3]=0.5 #R0
spar[0,4]=3.54 #phi
spar[0,5]=0.9 #ellp
spar[0,6]=0.5 #ns

spar[1,0]=100. #I0
spar[1,1]=-0.159
spar[1,2]=0.119
spar[1,3]=0.2 #R0
spar[1,4]=112.19 #phi
spar[1,5]=0.3 #ellp
spar[1,6]=1.0 #ns



#lens param
bsie,xi,yi,pa,q=1.6,  0.317     , -0.101     ,  0.1,  0.8

#foreground lens light
p_phot=[1000., xi, yi, 1.3, pa*180./np.pi, \
q,2.2,0.0001]
#p_phot= I0, xi,yi, r0,pa,q, sersic_n,I_back

z_lens=0.5869
z_source=2.4504
z_sub=0.3
msub=1e11
subx=55.87-100.
suby=96.37-100.

np.savez("./mock/param_Run"+str(Num).zfill(3)+".npz",
        spar=spar.reshape(-1),lpar=[bsie,xi,yi,pa,q],p_phot=p_phot,
        z_lens=z_lens,z_source=z_source,z_sub=z_sub,
        msub=msub,subx=subx,suby=suby)
