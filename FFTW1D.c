/*
| Computes Fourier transform 1D, this is just learning how to use FFTW
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fftw3.h>

#define error(x)		{printf("\n\nError generating,creating or opening "x"\n\n");exit(-1);}
#define errorrc(x)		{printf("\n\nError reading %s\nMaybe file does not exist\n\n",x);exit(-1);}
#define errorwc(x)		{printf("\n\nError generating,creating or writing %s\n\n",x);exit(-1);}
#define CLOSEFILE(x)            {fclose(x); x = NULL;}
#define SQR(x)			((x)*(x))
#define FREEP(x)		{free(x); x = NULL;}
#define sfsg			{printf("\n\n So far, so good...");getchar();printf("\n\n");}

typedef long int LI;
typedef unsigned long int ULI;

/****** global variables ******/

double pi, sqrt2pi, tempr;
fftw_plan plan_uhatf, plan_uhatb, plan_phat0f, plan_phat0b, plan_S1b, plan_S1f, plan_dydtb, plan_dydtf, plan_chif;

int main(){

	LI i, N, L;
        extern fftw_plan plan_uhatf, plan_uhatb;
	extern double pi;
	double *uhat; /* arrays */
	double *k;
	char syscall[180], namex[60], namek[60], namey[60];
	FILE *fx, *fk, *fy;

	pi = 4.0*atan(1.0);

	// Grid size
	N = (LI) 256;

	// Allocating necessary arrays
	// Transform is going to be from u to uhat? Or is it in place?

  if( (uhat = (double*) fftw_malloc(sizeof(double) * N)) == NULL)
		error("vector uhat");

	/* array containing frequencies in each position to use with FFTW3 and its half-complex vectors */
	// What does this array do? Reorganize frequencies in normal order?
	/*
	    k[0]=0.0;
	    k[N/2]=(2.0*pi/L)*(double)(N/2);

		for(i=1;i<N/2;i++)
		{
		    k[i]=(2.0*pi/L)*(double)i;
		    k[N-i]=(2.0*pi/L)*(double)i;
	    }
	*/

	/** initialize FFTW **/
	plan_uhatf = fftw_plan_r2r_1d(N, uhat, uhat, FFTW_R2HC, FFTW_MEASURE);
	plan_uhatb = fftw_plan_r2r_1d(N, uhat, uhat, FFTW_HC2R, FFTW_MEASURE);


	/* writing initial uhat(t) in the buffer, also initializing uold */
	/* check if initializing uhat from a previous run - if not, start it as zero */

	/* rect function */
	/*
		L=20;
		for(i=0;i<L;i++)
				uhat[i]=1.0;
		for(i=L;i<N;i++)
				uhat[i]=0.0;
	*/

	/* sine function */
	/*
		L = 10;
		for(i=0;i<N;i++)
			uhat[i]=sin(i/(double)L);
	*/

	/* data for uhat array as sinc function */
	L = 10;
	uhat[0]=1./pi/(double)L;
	for(i=1;i<N;i++)
		uhat[i]=sin(i/(double)L)/pi/i;

 /* open output files */
 sprintf(namex,"rect_x.dat");
 if( (fx = fopen(namex,"w")) == NULL)
   errorwc(namex);
 sprintf(namek,"rect_k.dat");
 if( (fk = fopen(namek,"w")) == NULL)
	 errorwc(namek);
 sprintf(namey,"rect_y.dat");
 if( (fy = fopen(namey,"w")) == NULL)
	 errorwc(namey);

  for(i=0;i<N;i++)
    fprintf(fx, "%lf \n",uhat[i]);

  fftw_execute(plan_uhatf); /* forward transform uhat */

  for(i=0;i<N;i++)
	  fprintf(fk, "%lf \n",uhat[i]);

	fftw_execute(plan_uhatb); /* backward transform uhat */

	for(i=0;i<N;i++)
		fprintf(fy, "%lf \n",uhat[i]/(double)N); /* has to normalize */

  /* Closing new FFTW plans and freeing remaining memory */
  /* Closing some FFTW plans and freeing some memory */

  fftw_destroy_plan(plan_uhatb);
  fftw_destroy_plan(plan_uhatf);
  fftw_free(uhat);
  CLOSEFILE(fx);
  CLOSEFILE(fk);
  CLOSEFILE(fy);

  printf("\n\n Transforms are done \n\n\a");

	double *ux, *uy;
	double max, oneu;

	if( (ux = (double*) malloc(sizeof(double) * N)) == NULL)
		error("vector ux");
	if( (uy = (double*) malloc(sizeof(double) * N)) == NULL)
		error("vector uy");

  printf("\n\n Now checking results \n\n\a");

	/* open same files for reading */
	sprintf(namex,"rect_x.dat");
	if( (fx = fopen(namex,"r")) == NULL)
		errorwc(namex);
	sprintf(namey,"rect_y.dat");
	if( (fy = fopen(namey,"r")) == NULL)
		errorwc(namey);

	/* read one file in real space */
	for(i=0;i<N;i++){
	  fscanf(fx, "%lf", &ux[i]);
  }

	/* read other file in real space */
	for(i=0;i<N;i++){
	  fscanf(fy, "%lf", &uy[i]);
  }

	/* calculate difference between two files */
	for(i=0;i<N;i++){
	  ux[i] = fabs( ux[i] - uy[i] );
  }

	/* find maximum difference */
	max=0;
	for(i=0;i<N;i++){
		oneu = ux[i];
		if(oneu > max)
			max = oneu;
  }
	printf("Max difference is: %e \n",max);

	CLOSEFILE(fx);
  CLOSEFILE(fy);
	FREEP(ux);
	FREEP(uy);

	printf("\n\n And we are done \n\n\a");

return 0;
}
