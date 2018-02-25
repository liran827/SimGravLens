#resimulate analytical model
import numpy as np
import matplotlib.pyplot as plt
from mock_lib import *
import cosmo
import lensinglib as ll
import time
try:
    import astropy.io.fits as pyfits
except:
    import pyfits
from scipy.interpolate import interp2d

#=============================
# import function from c lib
#=============================
#import ctype
from ctypes import *
lenslib=CDLL("./clens/liblens.so")
construct_lens_plane=lenslib.construct_lens_plane
mapping_source=lenslib.mapping_source
mapping_source_pix=lenslib.mapping_source_pix
from pylens import *

#==============================
# Defined FUNCTIONS
#==============================

def interp_grid(xp,yp,zp,xdev,ydev,gridsize,intp_range=200):
    nnn=np.sqrt(xp.shape[0])

    x_2d=xp.reshape(nnn,nnn)
    y_2d=xp.reshape(nnn,nnn)

    x_2d_temp=x_2d[nnn/2-200:nnn/2+200,nnn/2-200:nnn/2+200]
    y_2d_temp=y_2d[nnn/2-200:nnn/2+200,nnn/2-200:nnn/2+200]



    ix=np.floor(nx) + nnn/2
    iy=np.floor(ny) + nnn/2

    ind=np.where()




#=================================================================


#set observation parameters
expo=200.#200s expo time
nu_0=0.03 #ADU/s/pixel
sec2rad=1./180.*np.pi/3600.

#simulation configuration
lens_scale=13

OUTPUT=True
ADD_FOREGROUND=True
ADD_NOISE=True


if(lens_scale<10):
    print "warning: scale<5"


for i in [35]:


    RunNum=str(i).zfill(3)

    if(i<30):
        print "for i<30, please use sim_analytical.py instead"
        sys.exit()
        #if filenumber < 30, using the old one plane ray-tracing code
        #if filenumber >30, using this two plane code

    base="./mock/"
    para=np.load(base+"param_Run"+str(RunNum).zfill(3)+".npz")


    #subhalo positions, image center as 0,0
    subx=para["subx"]+0.
    suby=para["suby"]+0.

    msub=para["msub"]+0.
    z_lens=para["z_lens"]+0
    z_source=para['z_source']+0
    z_sub=para['z_sub']+0


    if(lens_scale<10):#if run a test
        msub=1e8
        OUTPUT=False


    #main lens parameters, SIE model
    bsie,xi,yi,pa,q = para['lpar'].astype("float32")
    #pa position angle, q: ellipticity



    #============================
    #      SET Lens plane
    #============================

    img_scale = 1 #if parameter >1, image resolution is increased
    nnn = 131*lens_scale #region to perform fft
    dpix_image = 0.04/img_scale # image resolution in arcsec
    boxsize = cosmo.Da(0.,z_lens)*dpix_image*sec2rad*nnn * 1000. #physical size kpc/h
    lensSize = sec2rad*dpix_image*nnn #angular size rad
    cellsize = boxsize/nnn #cellsize kpc/h


    #=============================
    #    MAKE LENS
    #=============================
    def rot_matrix(t):
        return np.array([[np.cos(t),np.sin(t)],[-np.sin(t),np.cos(t)]])

    print "make lens"
    print

    sigma_crit = cosmo.sigma_crit(z_lens,z_source) #critical density

    #set main lens mesh
    xmesh,ymesh = np.meshgrid( (np.arange(nnn)-nnn/2.+0.5)*dpix_image , (np.arange(nnn)-nnn/2.+0.5)*dpix_image)
    xp,yp = np.dot(rot_matrix(pa),np.array([xmesh.reshape(-1,order='C') - xi,ymesh.reshape(-1,order='C') - yi]) ) #arcsec

    #surface density for SIE
    sigma2d=np.zeros_like(xp)
    sigma2d=(sigma_crit/1e4) *np.sqrt(q)/2.*bsie/np.sqrt(np.power(xp,2.)+np.power(yp*q,2.))

    #calculate deflection angle
    #dpdx1,dpdy1= DFL_angle(sigma2d, z_lens, z_source, nnn,lensSize)


    lpar_sie=bsie,xi,yi,pa*180./np.pi,q,0.
    dpdx_sie,dpdy_sie,dpdx_sie2,dpdx_sie2=def_sie(xmesh.reshape(-1,order='C'),ymesh.reshape(-1,order='C'),lpar_sie) #arcsec
    dpdx1=-dpdy_sie/(180./np.pi*3600.) #rad
    dpdy1=-dpdx_sie/(180./np.pi*3600.)
    print "Using analytical SIE model for main lens"

    #set a subhalo mesh
    sub2d=np.zeros(nnn*nnn)
    sub2d=NFW_mesh(msub,z_sub,subx,suby,nnn,cellsize)  #surface density for a sub NFW halo
    #dpdx_a,dpdy_a= DFL_angle(sigma2d+sub2d, zsub, z_source, nnn,lensSize) #total deflection angle

    dpdx_s,dpdy_s= DFL_angle(sub2d, z_sub, z_source, nnn,lensSize)



    print "lens plane done"
    print


    #=====================================
    #Test model, comment when not needed
    #=====================================

    #lpar_sie=bsie,xi,yi,pa*180./np.pi,q,0.
    #dpdx_sie,dpdy_sie,dpdx_sie2,dpdx_sie2=def_sie(xmesh.reshape(-1,order='C'),ymesh.reshape(-1,order='C'),lpar_sie) #arcsec

    #plt.figure()
    #xx=(np.arange(nnn)-(nnn/2.-0.5) )*dpix_image
    #plt.plot(xx,180./np.pi*3600.*dpdy1.reshape(nnn,nnn)[nnn/2+10,:])
    #plt.plot(xx,-dpdx_sie.reshape(nnn,nnn)[nnn/2+10,:])
    #plt.xlim(xmax=5,xmin=-5)
    #diff=180./np.pi*3600.*dpdy1.reshape(nnn,nnn)[nnn/2+10,:]+dpdx_sie.reshape(nnn,nnn)[nnn/2+10,:]
    #ind=np.where((xx<5)*(xx>-5))
    #diff=diff[ind]
    #plt.figure()
    #plt.plot(xx[ind],diff)
    #plt.plot(xx,dpdy1.reshape(nnn,nnn)[nnn/2,:]*180./np.pi*3600.+dpdx_sie.reshape(nnn,nnn)[nnn/2,:],'--')

    #print
    #print "warning: already using analytical sie model"
    #print
    #dpdx1=-dpdy_sie/(180./np.pi*3600.) #rad
    #dpdy1=-dpdx_sie/(180./np.pi*3600.)


#==================================
#add deflection angle of two planes
#==================================

    if(z_sub<z_lens):
        Dfact=cosmo.Da(z_sub,z_lens)/cosmo.Da(0,z_lens)
        xmesh_dev,ymesh_dev=xmesh-Dfact*dpdx_s.reshape(nnn,nnn)*(180./np.pi*3600.),ymesh-Dfact*dpdy_s.reshape(nnn,nnn)*(180./np.pi*3600.)

        nx_dev_temp=np.round( xmesh[nnn/2-200:nnn/2+200,nnn/2-200:nnn/2+200]/dpix_image + (nnn/2 - 0.5) )
        ny_dev_temp=np.round( ymesh[nnn/2-200:nnn/2+200,nnn/2-200:nnn/2+200]/dpix_image + (nnn/2 - 0.5) )

        dpdx1_dev=dpdx1
        dpdy1_dev=dpdy1

        #only change the central
        dpdx1_dev.reshape(nnn,nnn)[nnn/2-200:nnn/2+200,nnn/2-200:nnn/2+200]=dpdx1.reshape(nnn,nnn)[ny_dev_temp.astype("int32"),nx_dev_temp.astype("int32")]
        dpdy1_dev.reshape(nnn,nnn)[nnn/2-200:nnn/2+200,nnn/2-200:nnn/2+200]=dpdy1.reshape(nnn,nnn)[ny_dev_temp.astype("int32"),nx_dev_temp.astype("int32")]

        dpdx_a=dpdx1_dev+dpdx_s
        dpdy_a=dpdy1_dev+dpdy_s


    else:
        print "only coding z_sub<z_lens case!!"
        sys.exit()



    #======================================
    #  MAKE SOURCE PLANE
    #======================================

    print "make source"

    #SourceSize=lensSize # in rad
    nnn_source=801 #odd
    dpix_source=0.01
    SourceSize=dpix_source*nnn_source*sec2rad

    n_source=para['spar'].shape[0]/7
    sourcemodel=['sersic']*n_source
    source_param=para['spar']
    n_x_source=nnn_source
    n_y_source=nnn_source
    source_xbase=np.outer(np.ones(n_y_source), np.arange(n_x_source)-n_x_source/2.)*dpix_source
    source_ybase=np.outer(np.arange(n_y_source)-n_y_source/2., np.ones(n_x_source))*dpix_source
    I_source=source_model(source_xbase, source_ybase, source_param, sourcemodel, n_source=n_source)
    source2d=I_source.reshape((1,nnn_source*nnn_source),order='C')[0]

    print "source done"
    print
    #================================

    print "mapping start"
    print
    #IMAGE position setting

    NX=200
    NY=200
    x_start=nnn/2-NX/2
    x_end=nnn/2+NX/2
    y_start=nnn/2-NY/2
    y_end=nnn/2+NY/2

    image2d=np.zeros(NX*NY)
    image2d_a=np.zeros(NX*NY)
    image2d_b=np.zeros(NX*NY)

    # image with no sub, green function
    mapping_source_pix(dpdx1.ctypes.data_as(POINTER(c_double)), dpdy1.ctypes.data_as(POINTER(c_double)), \
    c_int(x_start), c_int(x_end), c_int(y_start), c_int(y_end), c_double(lensSize), c_int(nnn),\
    c_double(SourceSize),c_int(nnn_source), source2d.ctypes.data_as(POINTER(c_double)),image2d.ctypes.data_as(POINTER(c_double)))
    #image with sub
    mapping_source_pix(dpdx_a.ctypes.data_as(POINTER(c_double)), dpdy_a.ctypes.data_as(POINTER(c_double)), \
    c_int(x_start), c_int(x_end), c_int(y_start), c_int(y_end), c_double(lensSize), c_int(nnn),c_double(SourceSize),\
    c_int(nnn_source), source2d.ctypes.data_as(POINTER(c_double)),image2d_a.ctypes.data_as(POINTER(c_double)))
    #image with no sub, analytical
    #mapping_source_pix(dpdy2.ctypes.data_as(POINTER(c_double)), dpdx2.ctypes.data_as(POINTER(c_double)), c_int(x_start), c_int(x_end), c_int(y_start), c_int(y_end), c_double(lensSize), c_int(nnn),c_double(SourceSize),c_int(nnn_source), source2d.ctypes.data_as(POINTER(c_double)),image2d_a.ctypes.data_as(POINTER(c_double)))

    print "mapping finished"
    print


    #========================
    # Making PSF
    #========================
    from scipy.signal import fftconvolve
    PSF=pyfits.getdata("acs_F606W_lens_1x.fits")

    if(img_scale !=1):
        PSF_new=np.zeros([PSF.shape[0]*img_scale, PSF.shape[1]*img_scale])

        NX_PSF=PSF.shape[0]*img_scale
        NY_PSF=PSF.shape[0]*img_scale

        for i in range(NX_PSF):
            for j in range(NY_PSF):
                PSF_new[i,j]=PSF[i/img_scale,j/img_scale]

        PSF_new/=img_scale**2
    else:
        PSF_new=PSF


    #===========adding the lens galaxy==========

    p_phot=para['p_phot']

    if(ADD_FOREGROUND):
        n_x_lens=NX
        n_y_lens=NY
        lens_xbase=np.outer(np.ones(n_y_lens), np.arange(n_x_lens)-n_x_lens/2.)*dpix_image
        lens_ybase=np.outer(np.arange(n_y_lens)-n_y_lens/2., np.ones(n_x_lens))*dpix_image
        I_phot=sersic_phot(lens_xbase, lens_ybase, p_phot)
        I_phot*=image2d.max()/I_phot.max()*2.
        ind=np.where(I_phot<0)
        I_phot[ind[0],ind[1]]=0


        img_c=fftconvolve(I_phot[:NY,:NX]+image2d.reshape(NX,NY,order='C'),PSF_new,mode='same')
        img_c2=fftconvolve(I_phot[:NY,:NX]+image2d_a.reshape(NX,NY,order='C'),PSF_new,mode='same')
        if(ADD_NOISE):
            img_c_ob=np.random.poisson(img_c*expo)/expo
            img_c2_ob=np.random.poisson(img_c2*expo)/expo
            noise_map=img_c-img_c_ob
            noise_map2=img_c2-img_c2_ob

    #===========adding the noise===============
    if(ADD_NOISE):
        Npix=img_c.shape[0]
        nu=nu_0*expo  #electron per pix in expo s
        I_noise=(np.random.poisson(nu,Npix*Npix).reshape(Npix,Npix)-nu)/expo

        img_c_ob+=I_noise
        img_c2_ob+=I_noise
        noise_map+=I_noise
        noise_map2+=I_noise





    #==============================
    #PLOTTING
    #==============================
    plt.figure()
    plt.imshow(np.log10(I_phot),origin='lower')
    plt.colorbar()

    plt.figure()
    plt.imshow(np.log10(I_source),origin='lower')
    plt.colorbar()


    plt.figure()
    plt.imshow(np.log10((sigma2d+sub2d).reshape(nnn,nnn)[y_start:y_end,x_start:x_end]/(sigma_crit/1e4)),origin='lower')
    plt.colorbar()


    plt.figure()
    plt.imshow(image2d.reshape(NX,NY,order='C'),origin='lower')
    plt.colorbar()


    plt.figure()
    plt.imshow(image2d_a.reshape(NX,NY,order='C'),origin='lower')
    plt.colorbar()



    plt.figure()
    plt.imshow(img_c_ob,origin='lower',cmap='gray')
    plt.colorbar()


    plt.figure()
    plt.imshow(img_c2_ob,origin='lower',cmap='gray')
    plt.colorbar()
    #===============================================================
    # OUTPUT TO FITS FILES
    #===============================================================


    if(OUTPUT):
        outname="./mock/Sim_RUN"+RunNum+".fits"
        prihdr=pyfits.Header()
        prihdr.set("dpix_l",dpix_image)
        prihdr.set("dpix_s",dpix_source)
        prihdr.set("N_IMG",nnn)
        prihdr.set("N_source",source2d.shape[0])
        prihdr.set("zlens",z_lens)
        prihdr.set("zsource",z_source)
        prihdr.set("zsub",z_sub)
        prihdr.set("expo",expo)
        prihdr.set("sky_nu_0",nu_0)

        hdu0=pyfits.PrimaryHDU(header=prihdr)

        hdu1=pyfits.ImageHDU(img_c2_ob) #with sub
        hdu2=pyfits.ImageHDU(img_c_ob)  # with nosub
        #hdu3=pyfits.ImageHDU(dpdx1.reshape(nnn,nnn))
        #hdu4=pyfits.ImageHDU(dpdy1.reshape(nnn,nnn))
        hdu3=pyfits.ImageHDU(noise_map)
        hdu4=pyfits.ImageHDU(noise_map2)
        #hdu5=pyfits.ImageHDU(PSF_new)
        #hdu6=pyfits.ImageHDU(I_phot)
        #hdu7=pyfits.ImageHDU(source2d.reshape(nnn_source,nnn_source,order='C'))
        #hdu8=pyfits.ImageHDU((sigma2d+sub2d).reshape(nnn,nnn,order='C'))

        hdulist=pyfits.HDUList([hdu0,hdu1,hdu2,hdu3,hdu4])
        hdulist.writeto(outname,clobber=True)
