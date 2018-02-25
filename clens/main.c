#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "icosmo.h"




int construct_lens_plane(double* sigma2d, double lensAngularDistance, double source_lensAngularDistance, double sourceAngularDistance, double lensSize, int lensNumberSide, double* dpdx, double* dpdy);
//double* mapping_source(double* dpdx, double* dpdy, int x_start, int x_end, int y_start, int y_end,double center_x, double center_y, double ss, double lensSize, int lensNumberSide, double lensAngularDistance, double sourceAngularDistance);
void mapping_source(double* dpdx, double* dpdy, int x_start, int x_end, int y_start, int y_end,double center_x, double center_y, double ss, double lensSize, int lensNumberSide, double lensAngularDistance, double sourceAngularDistance, double* image2d);
double SIS(double r, double disp,double rlim);
void input_lens2(const char fname[], double *sigma,int lensNumberSide);
void density_mesh_physical_scale(double *sigma, int lensNumberSide, double lensSize, double lensAngularDistance);
float endian(float f);


double source_plane_gauss(const double theta_x, const double theta_y, const double center_x, const double center_y, const double ss);
void mapping_source_pix(double* dpdx, double* dpdy, int x_start, int x_end, int y_start, int y_end, double lensSize, int lensNumberSide, double SourceSize, int SourceNumberSide,double* source2d, double* image2d);


void output_map_to_fits(const char fname[], int nrows, int ncols, double* mesh);
char lensFilename[200];

int main(void) {

	int    lensNumberSide,SourceNumberSide;
	double lensSize,SourceSize;
	double lensAngularDistance;
	double sourceAngularDistance;
	double source_lensAngularDistance;
	double center[2];

	int x_start, x_end;
	int y_start, y_end;
	double center_x, center_y, ss;


	double *sigma2d;
	double *image2d;
    double *source2d;

	double* dpdx;
	double* dpdy;

	double z_lens=0.3;
	double z_source=1.0;

    int i=0,j=0;

	lensNumberSide = 800;
    SourceNumberSide=10000;



	lensAngularDistance = Dang(0.,z_lens)*1000.;//unit kpc  
	sourceAngularDistance = Dang(0.,z_source)*1000.;//unit kpc         
	source_lensAngularDistance = sourceAngularDistance - lensAngularDistance*(1+z_lens)/(1+z_source);
	lensSize = 100./lensAngularDistance;//in rad
    SourceSize= 1.5*lensSize;
	
    sigma2d = (double*)malloc((int)sizeof(double)*lensNumberSide*lensNumberSide);
	dpdx = (double*)malloc( (int)sizeof(double)*lensNumberSide*lensNumberSide );
	dpdy = (double*)malloc( (int)sizeof(double)*lensNumberSide*lensNumberSide );
	image2d = (double*)malloc( (int)sizeof(double)*lensNumberSide*lensNumberSide );

	x_start=50;
	x_end=lensNumberSide-1-x_start;
	y_start=50;
	y_end=lensNumberSide-1-y_start;
	center_x=0.;
	center_y=0.;
	ss=0.05;
    int hN_Source=SourceNumberSide/2.;


	//sprintf(lensFilename,"../../sph_openmp/output_files/sdens_Aq_sph_main90.bin");
	density_mesh_physical_scale(sigma2d, lensNumberSide, lensSize, lensAngularDistance);

	construct_lens_plane(sigma2d, lensAngularDistance,source_lensAngularDistance, sourceAngularDistance, lensSize, lensNumberSide,  dpdx, dpdy);

    // make a pixelized source
    double SourceCellSize=0.;
	source2d = (double*)malloc((int)sizeof(double)*SourceNumberSide*SourceNumberSide);
    SourceCellSize= SourceSize/SourceNumberSide; //in rad
    for(i=0;i<SourceNumberSide;i++)
    {
        for(j=0;j<SourceNumberSide;j++)
        {
            source2d[i*SourceNumberSide+j]= source_plane_gauss((i-hN_Source)*SourceCellSize, (j-hN_Source)*SourceCellSize, 0.,0.,ss);

        }
    }

    mapping_source_pix(dpdx, dpdy, x_start, x_end, y_start, y_end, lensSize, lensNumberSide, SourceSize,SourceNumberSide, source2d, image2d);
    printf("hello!\n");

	//mapping_source(dpdx, dpdy, x_start, x_end, y_start, y_end, center_x, center_y, ss, lensSize, lensNumberSide, lensAngularDistance, sourceAngularDistance,image2d); // for gaussian source with dispersion of ss
    //

	char fname [50];
	sprintf(fname,"result");
	printf("%s\n",fname);
	output_map_to_fits(fname,x_end-x_start,y_end-y_start,image2d); //output sigma map of the lens
	//output_map_to_fits(fname,lensNumberSide,lensNumberSide,image2d); //output sigma map of the lens
    //
	sprintf(fname,"source");
	printf("%s\n",fname);
	output_map_to_fits(fname,SourceNumberSide,SourceNumberSide,source2d); //output sigma map of the lens
	//output_map_to_fits(fname,lensNumberSide,lensNumberSide,image2d); //output sigma map of the lens


	sprintf(fname,"dpdx");
	printf("%s\n",fname);
	output_map_to_fits(fname,lensNumberSide,lensNumberSide,dpdx); //output sigma map of the lens
	
	sprintf(fname,"dpdy");
	printf("%s\n",fname);
	output_map_to_fits(fname,lensNumberSide,lensNumberSide,dpdy); //output sigma map of the lens


}



void density_mesh_physical_scale(double *sigma, int lensNumberSide, double lensSize, double lensAngularDistance)
{
	int m, n, hN, N;
	double box  = lensSize * lensAngularDistance; //physical size of the box
	double hbox = box/2.0;
	double cellsize=box/lensNumberSide; //physical size of each cell
	N = lensNumberSide;
	hN = N/2;
	double PI = acos(-1.0);
	// using mesh file
	//char lensFilename [100];
	
    //input_lens2(lensFilename, sigma,lensNumberSide); //sigma is 1d array
	
	   for(m=0;m<lensNumberSide*lensNumberSide;m++)
	   {
	   sigma[m]*=0.; //NEED DENSITY not PARTICLE NUMBER in unit of 10^10/kpc^2 h*h
	//sigma[m]/=pow(cellsize,2.0); //NEED DENSITY not PARTICLE NUMBER!
	}
	
	

	printf("sigma=%e\n",sigma[(hN-1)*N+(hN-1)]);
	printf("cellsize=%e, box=%e\n",cellsize,box);


	double disp1=350.;
	double cx1=hN-0.5;//center of big halo
	double cy1=hN-0.5;
	double r1=0.;
	for (m=0; m<lensNumberSide; m++) {
		for (n=0; n<lensNumberSide; n++) {
			r1=sqrt( pow( (m-cx1)*cellsize,2.) +  pow( (n-cy1)*cellsize,2.) );

			sigma[m*lensNumberSide + n]+= SIS(r1,disp1,200.)/1e10;
		}
	}

	char fname [50];
	sprintf(fname,"lenstemp_Aq_sub_test");

	output_map_to_fits(fname,lensNumberSide,lensNumberSide,sigma); //output sigma map of the lens

}

void output_map_to_fits(const char fbase[], int nrows, int ncols, double* mesh)
{
    int i,j, n, nd,nb;
    int Ndata;
    int Nblock;
    char header[36*80];
    char fname[190];
    char stNAXIS1[80], stNAXIS2[80], stDate[80], stTITLE[80], stAUTHOR[80];
    
    float Dblock[36*20];
    double temp;
    FILE *fout;

    memset(header,' ', 36*80);
    memset(Dblock,'0', 36*80);

    sprintf(fname,"%s.fits", fbase);
    if ( !(fout = fopen(fname, "w") ) ) {
        printf("cannot open the file `%s'\n", fname);
        exit(0);
    }
    
    sprintf(stNAXIS1, "NAXIS1  = %20d ", nrows );
    sprintf(stNAXIS2, "NAXIS2  = %20d ", ncols );
    
    ///////////////////////000000000111111111122222222223//////////
    ///////////////////////123456789012345678901234567890//////////
    memcpy (header      , "SIMPLE  =                    T ", 30 );
    memcpy (header +1*80, "BITPIX  =                  -32 ", 30 );
    memcpy (header +2*80, "NAXIS   =                    2 ", 30 );
    memcpy (header +3*80, stNAXIS1, 30 );
    memcpy (header +4*80, stNAXIS2, 30 ); 
    memcpy (header +8*80, "END", 3); 

    fwrite (header , 1 , 36*80*sizeof(char) , fout );
    
    Ndata = nrows * ncols * sizeof(float);
    Nblock= Ndata / (36*80*sizeof(char) );



    nd = 0;  
    for (n=Nblock; n>0; n--) {
        for (nb =0; nb<36*20; nb++) {
            Dblock[nb] = endian((float)mesh[nd]);
            nd ++;
        }
        fwrite (Dblock, 1 , 36*80*sizeof(char), fout );
    }

    memset(Dblock,'0', 36*80);
    for (nb =0; nb<36*20; nb++) {
        Dblock[nb] = endian((float)mesh[nd]);
        nd ++;
        if (nd >= nrows * ncols )
            break;
    }
    fwrite (Dblock, 1 , 36*80*sizeof(char), fout );

    fclose(fout);
}




double SIS(double r, double disp,double rlim)
{
    	double PI = acos(-1.0);
	if(r>rlim)
		return 0;
	return pow(disp,2.)/2./G/r/1000;
}

void input_lens2(const char fname[], double *sigma,int lensNumberSide)
{
	//READ LENS MESH FILE FROM fname, STORE the DENSITY in sigma

	float *temp;
	temp=(float*)malloc((int)sizeof(float)*lensNumberSide*lensNumberSide);
	FILE *fin;
	if(!(fin=fopen(fname,"r")))
	{
		printf("lens file cannot open\n");
		exit(0);
	}

	fread(temp, sizeof(float),lensNumberSide*lensNumberSide,fin);
	//printf("%d\n",lensNumberSide);
	int i;
	for (i=0;i<lensNumberSide*lensNumberSide;i++)
		sigma[i]=temp[i];
	fclose(fin);	

}


