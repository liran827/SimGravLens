#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <math.h>
#include "lens.h"

//generate lens
//generate potential
//alpha
//source
//mapping

void grav2alpha(double *dpdx, double *dpdy, double source_lensAngularDistance, double sourceAngularDistance, int lensNumberSide)
{
    int i, j;
    double norm = (4*source_lensAngularDistance/sourceAngularDistance)/(LightSpeed*LightSpeed);
    for (i=0; i<lensNumberSide; i++) {
        for (j=0; j<lensNumberSide; j++) {
            dpdx[i*lensNumberSide+j] *= norm;
            dpdy[i*lensNumberSide+j] *= norm;
        }
    }
}

int construct_lens_plane(double* sigma2d, double lensAngularDistance, double source_lensAngularDistance, double sourceAngularDistance, double lensSize, int lensNumberSide, double* dpdx, double* dpdy)
{


    //static double *dpdx;
    //static double *dpdy;

    //dpdx = (double*)malloc( (int)sizeof(double)*lensNumberSide*lensNumberSide );
    isolated_gravity_plane_u(dpdx, sigma2d, lensNumberSide, lensSize*lensAngularDistance );

    //dpdy = (double*)malloc( (int)sizeof(double)*lensNumberSide*lensNumberSide );
    isolated_gravity_plane_v(dpdy, sigma2d, lensNumberSide, lensSize*lensAngularDistance );



    grav2alpha(dpdx, dpdy, source_lensAngularDistance, sourceAngularDistance, lensNumberSide);
    return 0;
}

double source_plane_gauss(const double theta_x, const double theta_y, const double center_x, const double center_y, const double ss)
{
	int nx, ny;
	double x, y, r;
	double PI = acos(-1.0);
	double sec2rad=PI/180./60./60.;
	double tx,ty;
	double cx,cy;



	//analytical source
    	cx =center_x * sec2rad;
    	cy =center_y * sec2rad;
	tx = (theta_x - cx); // rad
	ty = (theta_y - cy); // rad
	r=sqrt(tx*tx + ty*ty)/sec2rad;
	//printf("txty %f %f %f\n",tx,ty,r);
	if(r < 1) //r range in arcsec
		return exp(-r*r/2./(ss*ss)); //gaussian
	else
		return 0.0;

}

void mapping_source_pix(double* dpdx, double* dpdy, int x_start, int x_end, int y_start, int y_end, double lensSize, int lensNumberSide, double SourceSize, int SourceNumberSide, double* source2d, double* image2d)
{
	int nx, ny, Nx, Ny, Nx1, Ny1;
	int I0, I1, I2, I3, J0, J1, J2, J3;
	double x, y;
	double tx01, tx02, tx03, tx12, tx13;
	double ty01, ty02, ty03, ty12, ty13;
	double dx, dy, mu, lum, tr1, tr2, tr3, tr4;
    double SourceCellSize;
    int hN_Source;
	double min, ave, max;

    hN_Source=SourceNumberSide/2;


	Nx = x_end - x_start;
	Ny = y_end - y_start;

	Nx1 = Nx + 1;
	Ny1 = Ny + 1;


	static double *src_x;
	static double *src_y;

	src_x = (double*)malloc((int)sizeof(double)*Nx1*Ny1);
	src_y = (double*)malloc((int)sizeof(double)*Nx1*Ny1);

	dx = lensSize / lensNumberSide ;
	dy = lensSize / lensNumberSide ;

    SourceCellSize= SourceSize / SourceNumberSide;


    printf("hello,sourcecellsize=%e %e\n",SourceCellSize, SourceSize);

	for (nx=0; nx<Nx1; nx++) {
		x = (nx + x_start - 0.5*lensNumberSide + 0.5) * dx;
		for (ny=0; ny<Ny1; ny++) {
			y = (ny + y_start - 0.5*lensNumberSide + 0.5) * dy;
			src_x[nx*Ny1+ny] = x + dpdx[(nx+x_start)*lensNumberSide+(ny+y_start)];
			src_y[nx*Ny1+ny] = y + dpdy[(nx+x_start)*lensNumberSide+(ny+y_start)];
		}
	}


	//norm = dx * dy * ( sourceAngularDistance / lensAngularDistance );

	for (nx=0; nx<Nx; nx++){
		I0 = I3 = nx;
		I1 = I2 = nx + 1;

		for (ny=0; ny<Ny; ny++){
			J0 = J1 = ny;
			J2 = J3 = ny + 1;

      
			tx01 = src_x[I1*Ny1+J1] - src_x[I0*Ny1+J0] ;
			tx02 = src_x[I2*Ny1+J2] - src_x[I0*Ny1+J0] ;
			tx03 = src_x[I3*Ny1+J3] - src_x[I0*Ny1+J0] ;
			tx12 = src_x[I2*Ny1+J2] - src_x[I1*Ny1+J1] ;
			tx13 = src_x[I3*Ny1+J3] - src_x[I1*Ny1+J1] ;

			ty01 = src_y[I1*Ny1+J1] - src_y[I0*Ny1+J0] ;
			ty02 = src_y[I2*Ny1+J2] - src_y[I0*Ny1+J0] ;
			ty03 = src_y[I3*Ny1+J3] - src_y[I0*Ny1+J0] ;
			ty12 = src_y[I2*Ny1+J2] - src_y[I1*Ny1+J1] ;
			ty13 = src_y[I3*Ny1+J3] - src_y[I1*Ny1+J1] ;

			tr1 = (tx01*ty02 - tx02*ty01);
			if (tr1 < 0.0)
				tr1 = - tr1;

			tr2 = (tx03*ty02 - tx02*ty03);
			if (tr2 < 0.0)
				tr2 = - tr2;

			tr3 = (tx01*ty13 - tx13*ty01);
			if (tr3 < 0.0)
				tr3 = - tr3;

			tr4 = (tx13*ty12 - tx12*ty13);
			if (tr4 < 0.0)
				tr4 = - tr4;

			//mu = 4 * norm / ( tr1 + tr2 + tr3 + tr4 );

			//printf("lum= %d %d %f %f %f\n",I0*Ny1,J0,J1,J2,J3);
            //printf("src_x= %f \n", src_x[I0*Ny1+J0]/SourceCellSize );

			lum  = 0.0;
			lum += source2d[ (int)(src_x[I0*Ny1+J0]/SourceCellSize + hN_Source) * SourceNumberSide + (int)( src_y[I0*Ny1+J0] / SourceCellSize + hN_Source) ];
			lum += source2d[ (int)(src_x[I0*Ny1+J1]/SourceCellSize + hN_Source) * SourceNumberSide + (int)( src_y[I0*Ny1+J1] / SourceCellSize + hN_Source) ];
			lum += source2d[ (int)(src_x[I0*Ny1+J2]/SourceCellSize + hN_Source) * SourceNumberSide + (int)( src_y[I0*Ny1+J2] / SourceCellSize + hN_Source) ];
			lum += source2d[ (int)(src_x[I0*Ny1+J3]/SourceCellSize + hN_Source) * SourceNumberSide + (int)( src_y[I0*Ny1+J3] / SourceCellSize + hN_Source) ];
			lum /= 4.0;
			image2d[nx*Ny+ny] =  lum; //xy corrected LIRAN
			//image2d[ny*Nx+nx] =  lum; //xy corrected LIRAN
		}
	}
	free(src_x);
	free(src_y);

    printf("hello,sourcecellsize=%e\n",SourceCellSize);

	//return image2d;
}


void mapping_source(double* dpdx, double* dpdy, int x_start, int x_end, int y_start, int y_end,double center_x, double center_y, double ss, double lensSize, int lensNumberSide,  double* image2d)
{
	int nx, ny, Nx, Ny, Nx1, Ny1;
	int I0, I1, I2, I3, J0, J1, J2, J3;
	double x, y;
	double tx01, tx02, tx03, tx12, tx13;
	double ty01, ty02, ty03, ty12, ty13;
	double dx, dy, mu, lum, tr1, tr2, tr3, tr4;

	double min, ave, max;


	Nx = x_end - x_start;
	Ny = y_end - y_start;

	Nx1 = Nx + 1;
	Ny1 = Ny + 1;

	//double *image2d;
	//image2d = (double*)malloc( (int) sizeof(double) * (Nx) * (Ny) );


	static double *src_x;
	static double *src_y;

	src_x = (double*)malloc((int)sizeof(double)*Nx1*Ny1);
	src_y = (double*)malloc((int)sizeof(double)*Nx1*Ny1);

	dx = lensSize / lensNumberSide ;
	dy = lensSize / lensNumberSide ;

	for (nx=0; nx<Nx1; nx++) {
		x = (nx + x_start - 0.5*lensNumberSide + 0.5) * dx;
		for (ny=0; ny<Ny1; ny++) {
			y = (ny + y_start - 0.5*lensNumberSide + 0.5) * dy;
			src_x[nx*Ny1+ny] = x + dpdx[(nx+x_start)*lensNumberSide+(ny+y_start)];
			src_y[nx*Ny1+ny] = y + dpdy[(nx+x_start)*lensNumberSide+(ny+y_start)];
		}
	}


	//norm = dx * dy * ( sourceAngularDistance / lensAngularDistance );

	for (nx=0; nx<Nx; nx++){
		I0 = I3 = nx;
		I1 = I2 = nx + 1;

		for (ny=0; ny<Ny; ny++){
			J0 = J1 = ny;
			J2 = J3 = ny + 1;

			tx01 = src_x[I1*Ny1+J1] - src_x[I0*Ny1+J0] ;
			tx02 = src_x[I2*Ny1+J2] - src_x[I0*Ny1+J0] ;
			tx03 = src_x[I3*Ny1+J3] - src_x[I0*Ny1+J0] ;
			tx12 = src_x[I2*Ny1+J2] - src_x[I1*Ny1+J1] ;
			tx13 = src_x[I3*Ny1+J3] - src_x[I1*Ny1+J1] ;

			ty01 = src_y[I1*Ny1+J1] - src_y[I0*Ny1+J0] ;
			ty02 = src_y[I2*Ny1+J2] - src_y[I0*Ny1+J0] ;
			ty03 = src_y[I3*Ny1+J3] - src_y[I0*Ny1+J0] ;
			ty12 = src_y[I2*Ny1+J2] - src_y[I1*Ny1+J1] ;
			ty13 = src_y[I3*Ny1+J3] - src_y[I1*Ny1+J1] ;

			tr1 = (tx01*ty02 - tx02*ty01);
			if (tr1 < 0.0)
				tr1 = - tr1;

			tr2 = (tx03*ty02 - tx02*ty03);
			if (tr2 < 0.0)
				tr2 = - tr2;

			tr3 = (tx01*ty13 - tx13*ty01);
			if (tr3 < 0.0)
				tr3 = - tr3;

			tr4 = (tx13*ty12 - tx12*ty13);
			if (tr4 < 0.0)
				tr4 = - tr4;

			//mu = 4 * norm / ( tr1 + tr2 + tr3 + tr4 );

			lum  = 0.0;
			lum += source_plane_gauss(src_x[I0*Ny1+J0],  src_y[I0*Ny1+J0], center_x, center_y,ss);
			lum += source_plane_gauss(src_x[I1*Ny1+J1],  src_y[I1*Ny1+J1], center_x, center_y,ss);
			lum += source_plane_gauss(src_x[I2*Ny1+J2],  src_y[I2*Ny1+J2], center_x, center_y,ss);
			lum += source_plane_gauss(src_x[I3*Ny1+J3],  src_y[I3*Ny1+J3], center_x, center_y,ss);
			lum /= 4.0;
			image2d[nx*Ny+ny] =  lum;
			//image2d[ny*Nx+nx] =  lum;
			//printf("lum= %d %d %f %f %f\n",nx,ny,lum, center_y,ss);
		}
	}
	free(src_x);
	free(src_y);


	//return image2d;
}



typedef union
{
  float as_float;
  char  as_char[4];
} Float32;

float endian(float f)
{
    Float32 u32;
    char c;
    u32.as_float = f;
    c = u32.as_char[0];
    u32.as_char[0] = u32.as_char[3];
    u32.as_char[3] = c;

    c = u32.as_char[1];
    u32.as_char[1] = u32.as_char[2];
    u32.as_char[2] = c;
    return u32.as_float;
}


char *getstring(char *str,int num)
{
    printf("%s\n",str);
    return str;
}

char *reverse(char *str,int num)
{
    getstring(str,num);
    int half = num / 2;
    int i;
    char temp;
    for(i =0;i < half;++ i)
    {
        temp = str[num - 1 - i];
        str[num - 1 - i] = str[i];
        str[i] = temp;
    }
    printf("%s\n",str);
    return str;
}
