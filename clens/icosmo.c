#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include "icosmo.h"


//*************************************************************************
//cosmology lib
//*************************************************************************
double IntD(double a , void *p){

	double result = VC/H0* pow( (Omega_M0*a + Omega_Lambda0*pow(a , 4. )) , -0.5 );
	return result;
}

double Dang(double z1,double z2){
	//The function returns the angular distance between z1, z2 (z2>z1)

	gsl_function FDa;
	double alpha = 0.;
	double a1, a2;
	double result, error;
	int key=1;
	a1=1./ (z1+1.);
	a2=1./ (z2+1.);
	FDa.function = &IntD;
	FDa.params = &alpha;

	gsl_integration_workspace* w = gsl_integration_workspace_alloc (1000);
	gsl_integration_qag(&FDa,a2,a1,0,1.e-3,1000,key,w,&result,&error);
	gsl_integration_workspace_free (w);

	//printf("%g,%g,%g\n",result,H0,error);
	return result*a2;

}

double Dcom(double z1,double z2){
	//The function returns the comoving distance between z1, z2 (z2>z1)

	gsl_function FDa;
	double alpha = 0.;
	double a1, a2;
	double result, error;
	int key=1;
	a1=1./ (z1+1.);
	a2=1./ (z2+1.);
	FDa.function = &IntD;
	FDa.params = &alpha;

	gsl_integration_workspace* w = gsl_integration_workspace_alloc (1000);
	gsl_integration_qag(&FDa,a2,a1,0.,1.e-3,1000,key,w,&result,&error);
	gsl_integration_workspace_free (w);

	//printf("%g,%g,%g\n",result,H0,error);
	return result;

}
