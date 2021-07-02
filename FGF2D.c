/*
Generates 2D FGF, Out of place FFTW transform
gcc -o x FGF2D.c -lfftw3 -lm
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>
#include <time.h>

#define error(x)      {printf("\n\nError generating,creating or opening "x"\n\n");exit(-1);}
#define errorrc(x)    {printf("\n\nError reading %s\nMaybe file does not exist\n\n",x);exit(-1);}
#define errorwc(x)    {printf("\n\nError generating,creating or writing %s\n\n",x);exit(-1);}
#define CLOSEFILE(x)  {fclose(x); x = NULL;}
#define SQR(x)        ((x)*(x))
#define FREEP(x)      {free(x); x = NULL;}
#define sfsg          {printf("\n\n So far, so good...");getchar();printf("\n\n");}
#define RAND()        (2.*RCTE*rand()/DRAND-RCTE)
#define ABS2(a,b)     ( SQR(((double)(a))) + SQR(((double)(b))) )
#define CRD(i,j)      (N2*(i)+(j)) // coordinates in 2D array

typedef long int LI;
typedef unsigned long int ULI;
double DRAND=(double)RAND_MAX;
double RCTE=.5*sqrt(12.);

/****** global variables ******/

fftw_plan plan_uf, plan_ub;

int main(){

	LI i, j, N, N2;
  extern fftw_plan plan_uf, plan_ub;
	extern double pi;
	double *u0, *K;
	fftw_complex *v0, *Phat; /* arrays */
	int dim;
	double dx,sqdx,H,Ltot;
	char name1[60], name2[60];
	FILE *f1, *f2;

	time_t t;
  /* Intializes random number generator */
  srand((unsigned) time(&t));

	// Grid size
	N = (LI) 1<<9; // 1<<N = 2^N
	N2 = (int)(N/2)+1;

	Ltot = 1.;
	dx = Ltot/(double)N;
	sqdx = sqrt(dx);
	H = 1./3.; // Holder exponent
	dim = 2;

	// Allocating necessary arrays
	if( (K = (double*) malloc(sizeof(double) * N)) == NULL)
		error("vector K");
  if( (u0 = (double*) fftw_malloc(sizeof(double) * N * N)) == NULL)
		error("vector u0");
	// v0 = F( u0 ), fftw format
	if( (v0 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N * N2 )) == NULL)
		error("vector v0");
	if( (Phat = (fftw_complex*) malloc(sizeof(fftw_complex) * N * N2 )) == NULL)
		error("vector Phat");

	/** initialize FFTW **/
	plan_uf = fftw_plan_dft_r2c_2d(N, N, u0, v0, FFTW_MEASURE);
	plan_ub = fftw_plan_dft_c2r_2d(N, N, v0, u0, FFTW_MEASURE);
	/* or maybe, with the in-place transform
	plan_uf = fftw_plan_r2r_2d(N, N, u0, u0, FFTW_R2HC, FFTW_R2HC, FFTW_MEASURE);
	plan_uf = fftw_plan_r2r_2d(N, N, u0, u0, FFTW_HC2R, FFTW_HC2R, FFTW_MEASURE);
	*/

	// Array of frequencies in Fourier space
  K[0]=0.0;
  K[N/2]=(double)(N/2);
	for(i=1;i<N/2;i++){
	  K[i]=(double)i;
	  K[N-i]=-(double)i;
  }

	// Assign random vector dW (white noise increments)
	for(i=0;i<N*N;i++){
		u0[i] = RAND() * dx;
	}

	// Save initial array to a text file

	sprintf(name1,"FGF2D_In.dat");
  if( (f1 = fopen(name1,"w")) == NULL)
    errorwc(name1);

	for(i=0;i<N*N;i++)
    fprintf(f1, "%le \n",u0[i]);

	CLOSEFILE(f1);


	// Assign PH operator directly in Fourier space
	// Imag. components are zero
	// frequencies are i/Ltot, normalization is eps = 1/Ltot
	for(i=0;i<N;i++){
		for (j = 0; j < N2; j++){
			//Phat[CRD(i,j)] = pow( ABS2(i,j)+1., -.5*H-.25*dim );
			Phat[CRD(i,j)] = pow( ABS2(K[i],K[j])+1., -.5*H-.25*dim );
		}
	}

	fftw_execute(plan_uf);

	// multiplication directly in Fourier space
	// Phat has no imaginary component, hence
	// u_k = Re u_k * Re P_k + I * Im u_k * Re P_k
	for(i=0;i<N*N2;i++)
		v0[i] = v0[i] * Phat[i];

	fftw_execute(plan_ub);

	// Save final array to a text file
	sprintf(name2,"FGF2D_Out.dat");
  if( (f2 = fopen(name2,"w")) == NULL)
    errorwc(name2);

	for(i=0;i<N*N;i++)
    fprintf(f2, "%le \n",u0[i]);

	CLOSEFILE(f2);

  fftw_destroy_plan(plan_ub);
  fftw_destroy_plan(plan_uf);
  fftw_free(u0);
	fftw_free(v0);
	fftw_free(Phat);
	FREEP(K)

	printf("\n\n And we are done \n\n\a");

return 0;
}
