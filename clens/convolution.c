#include "convolution.h"
#include <stdlib.h>
#include <fftw3.h>
#include <math.h>

fftw_complex* isolated_gravity_green_function_u(int num_side, double length_side) 
{
    int m,n;
    int Nside;
    double *gf_r, u,v, g0, Lside;
    double scale_zeropoint;
    double lp, spt,PI, r2, r;
    fftw_complex* gf_k;
    fftw_plan pf;
    scale_zeropoint = (double)num_side;
    Nside = 2*num_side;
    Lside = 2*length_side;
    g0 = (Lside*Lside)/((double)Nside*Nside)*43007.1;
    gf_r = (double*)malloc((int)Nside*Nside*sizeof(double));
    gf_k = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*Nside*(Nside/2+1));
    for (m=0; m<=num_side; m++)
        for (n=0; n<=num_side; n++) {
            u = m*length_side/num_side;
            v = n*length_side/num_side;
            r2 = (u*u + v*v);
            gf_r[m*Nside+n] =  - g0*(u/r2);
        }
    for (m=num_side+1; m<Nside; m++)
        for (n=0; n<=num_side; n++)  {
            gf_r[m*Nside+n] = -gf_r[(Nside-m)*Nside+n];
        }
    for (m=0; m<=num_side; m++)
        for (n=num_side+1; n<Nside; n++)  {
            gf_r[m*Nside+n] = gf_r[m*Nside+(Nside-n)];
        }
    for (m=num_side+1; m<Nside; m++)
        for (n=num_side+1; n<Nside; n++)  {
            gf_r[m*Nside+n] = -gf_r[(Nside-m)*Nside+(Nside-n)];
        }
    gf_r[0] = 0.0;
    pf = fftw_plan_dft_r2c_2d(Nside, Nside, gf_r, gf_k, FFTW_ESTIMATE);
    fftw_execute(pf);
    fftw_destroy_plan(pf);
    free(gf_r);
    return gf_k;
}

fftw_complex* isolated_gravity_green_function_v(int num_side, double length_side) 
{
    int m,n;
    int Nside;
    double *gf_r, u,v, g0, Lside;
    double scale_zeropoint;
    double lp, spt,PI, r2, r;
    fftw_complex* gf_k;
    fftw_plan pf;
    scale_zeropoint = (double)num_side;
    Nside = 2*num_side;
    Lside = 2*length_side;
    g0 = (Lside*Lside)/((double)Nside*Nside)*43007.1;
    gf_r = (double*)malloc((int)Nside*Nside*sizeof(double));
    gf_k = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*Nside*(Nside/2+1));
    for (m=0; m<=num_side; m++)
        for (n=0; n<=num_side; n++) {
            u = m*length_side/num_side;
            v = n*length_side/num_side;
            r2 = (u*u + v*v);
            gf_r[m*Nside+n] = - g0*(v/r2);
        }
    for (m=num_side+1; m<Nside; m++)
        for (n=0; n<=num_side; n++)  {
            gf_r[m*Nside+n] = gf_r[(Nside-m)*Nside+n];
        }
    for (m=0; m<=num_side; m++)
        for (n=num_side+1; n<Nside; n++)  {
            gf_r[m*Nside+n] = -gf_r[m*Nside+(Nside-n)];
        }
    for (m=num_side+1; m<Nside; m++)
        for (n=num_side+1; n<Nside; n++)  {
            gf_r[m*Nside+n] = -gf_r[(Nside-m)*Nside+(Nside-n)];
        }
    gf_r[0] = 0.0;
    pf = fftw_plan_dft_r2c_2d(Nside, Nside, gf_r, gf_k, FFTW_ESTIMATE);
    fftw_execute(pf);
    fftw_destroy_plan(pf);
    free(gf_r);
    return gf_k;
}

void isolated_gravity_plane_u(double *dpot_du, double *sigma, int num_side, double length_side)
{
    int i,j;
    int Nside;
    double PI;
    double green, g0;
    double ku, kv, k2, k_scale;
    double ggc, fu, fv, ff, ohPi;
    double Re_delta, Im_delta, Re_gf, Im_gf;
    double *dlt;
    fftw_complex *iso_gf,*pot_k;
    fftw_plan pf, pb;

    double prof, lp, spt, lsize;
    PI = acos(-1.0);

    Nside = 2*num_side;
    ohPi = PI/(double)Nside;
    g0  = 1.0/((double)Nside*Nside); // normalize density field

    iso_gf = isolated_gravity_green_function_u(num_side,length_side);
    dlt    = (double*)malloc((int)Nside*Nside*sizeof(double));
    for (i=0; i<Nside; i++) {
        for (j=0; j<Nside; j++) {
            if (i<num_side && j<num_side)
                dlt[i*Nside + j] = sigma[i*num_side + j];
            else
                dlt[i*Nside + j] = 0.0;
        }
    }
    pot_k= (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*Nside*(Nside/2+1));
    pf  = fftw_plan_dft_r2c_2d(Nside, Nside, dlt, pot_k, FFTW_ESTIMATE);
    fftw_execute(pf);
    fftw_destroy_plan(pf);

    for (i=0; i<Nside; i++) {
        ku  = (double)i;
        if (i > Nside/2) ku = (double)(i-Nside);
        for (j=0; j<=Nside/2; j++) {
            kv = (double)j;
            k2 = (kv*kv + ku*ku);

            green = 1.0;
            Re_delta = g0*pot_k[i * (Nside/2+1) + j][0];
            Im_delta = g0*pot_k[i * (Nside/2+1) + j][1];
            Re_gf = iso_gf[i * (Nside/2+1) + j][0];
            Im_gf = iso_gf[i * (Nside/2+1) + j][1];

            pot_k[i*(Nside/2+1)+j][0]=green*(Re_delta*Re_gf-Im_delta*Im_gf);
            pot_k[i*(Nside/2+1)+j][1]=green*(Re_delta*Im_gf+Im_delta*Re_gf);
        }
    }
    pb = fftw_plan_dft_c2r_2d(Nside,Nside, pot_k, dlt, FFTW_ESTIMATE);
    fftw_execute(pb);
    fftw_destroy_plan(pb);
    fftw_free(pot_k);
    fftw_free(iso_gf);

    for (i=0; i<num_side; i++) {
        for (j=0; j<num_side; j++) {
            dpot_du[i*num_side+j] = dlt[i*Nside + j];
        }
    }
    free(dlt);
}

void isolated_gravity_plane_v(double *dpot_dv, double *sigma, int num_side, double length_side)
{
    int i,j;
    int Nside;
    double PI;
    double green, g0;
    double ku, kv, k2, k_scale;
    double ggc, fu, fv, ff, ohPi;
    double Re_delta, Im_delta, Re_gf, Im_gf;
    double *dlt;
    fftw_complex *iso_gf,*pot_k;
    fftw_plan pf, pb;
    double prof, lp, spt, lsize;
    PI = acos(-1.0);

    Nside = 2*num_side;
    ohPi = PI/(double)Nside;
    g0  = 1.0/((double)Nside*Nside); // normalize density field

    iso_gf = isolated_gravity_green_function_v(num_side,length_side);
    dlt    = (double*)malloc((int)Nside*Nside*sizeof(double));
    for (i=0; i<Nside; i++) {
        for (j=0; j<Nside; j++) {
            if (i<num_side && j<num_side)
                dlt[i*Nside + j]=sigma[i*num_side + j];
            else
                dlt[i*Nside + j]=0.0;
        }
    }
    pot_k= (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*Nside*(Nside/2+1));
    pf  = fftw_plan_dft_r2c_2d(Nside, Nside, dlt, pot_k, FFTW_ESTIMATE);
    fftw_execute(pf);
    fftw_destroy_plan(pf);

    for (i=0; i<Nside; i++) {
        ku  = (double)i;
        if (i > Nside/2) ku = (double)(i-Nside);
        for (j=0; j<=Nside/2; j++) {
            kv = (double)j;
            k2 = (kv*kv + ku*ku);

            green = 1.0;
            Re_delta = g0*pot_k[i * (Nside/2+1) + j][0];
            Im_delta = g0*pot_k[i * (Nside/2+1) + j][1];
            Re_gf = iso_gf[i * (Nside/2+1) + j][0];
            Im_gf = iso_gf[i * (Nside/2+1) + j][1];

            pot_k[i*(Nside/2+1)+j][0]=green*(Re_delta*Re_gf-Im_delta*Im_gf);
            pot_k[i*(Nside/2+1)+j][1]=green*(Re_delta*Im_gf+Im_delta*Re_gf);
        }
    }
    pb = fftw_plan_dft_c2r_2d(Nside,Nside, pot_k, dlt, FFTW_ESTIMATE);
    fftw_execute(pb);
    fftw_destroy_plan(pb);
    fftw_free(pot_k);
    fftw_free(iso_gf);

    for (i=0; i<num_side; i++) {
        for (j=0; j<num_side; j++) {
            dpot_dv[i*num_side+j] = dlt[i*Nside + j];
        }
    }
    free(dlt);
}

