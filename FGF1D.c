/*
Generates 1D FGF
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
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

typedef long int LI;
typedef unsigned long int ULI;
double DRAND=(double)RAND_MAX;
double RCTE=.5*sqrt(12.);

/****** global variables ******/

fftw_plan plan_uf, plan_ub;

int main(){

	LI i, N;
  extern fftw_plan plan_uf, plan_ub;
	extern double pi;
	double *u0, *Phat; /* arrays */
	double dx,sqdx,H,Ltot;
	char name1[60], name2[60];
	FILE *f1, *f2;

	time_t t;
  /* Intializes random number generator */
  srand((unsigned) time(&t));

	// Grid size
	N = (LI) 1<<18; // 1<<N = 2^N

	Ltot = 1.;
	dx = Ltot/(double)N;
	sqdx = sqrt(dx);
	H = 1./3.; // Holder exponent

	// Allocating necessary arrays
  if( (u0 = (double*) fftw_malloc(sizeof(double) * N)) == NULL)
		error("vector u0");
	if( (Phat = (double*) malloc(sizeof(double) * N)) == NULL)
		error("vector Phat");

	/** initialize FFTW **/
	plan_uf = fftw_plan_r2r_1d(N, u0, u0, FFTW_R2HC, FFTW_MEASURE);
	plan_ub = fftw_plan_r2r_1d(N, u0, u0, FFTW_HC2R, FFTW_MEASURE);

	// Assign random vector dW (white noise increments)
	for(i=0;i<N;i++){
		u0[i] = RAND() * sqdx;
		//printf("%le\n",u0[i]);
	}

	// Save initial array to a text file
	sprintf(name1,"FGF1D_In.dat");
  if( (f1 = fopen(name1,"w")) == NULL)
    errorwc(name1);

	for(i=0;i<N;i++)
    fprintf(f1, "%le \n",u0[i]);

	CLOSEFILE(f1);

	// Assign PH operator directly in Fourier space
	// Imag. components are zero
	// frequencies are i/Ltot, normalization is eps = 1/Ltot
	for(i=0;i<N/2;i++)
		Phat[i] = pow( SQR((double)i)+1., -.5*H-.25 );

	fftw_execute(plan_uf);

	// multiplication directly in Fourier space
	// Phat has no imaginary component, hence
	// u_k = Re u_k * Re P_k + I * Im u_k * Re P_k

	/*
	No normalization is needed.
	If we discretize the expression for the continuous Fourier transform, we find
	that uhat(k) \approx dx uhat(k_j)
	That is, to compare a continuous and discrete Fourier transforms, we have
	to multiply by dx
	Instead, when working with rough velocity fields, like dW, and integrations
	over dW, as is the case for the action of the operator PH, no such
	normalization is needed
	P_H W = IF( |k|^{-H-1/2} F(W) )
	Where F and IF are continuous Fourier transforms and W is delta-correlated
	discretizing this, we find
	P_H W(x) \approx dk sum_j exp(2 pi i k_j x)|k_j|^{-H-1/2} sum_l w(dy_l) exp(-2 pi i k_j y_l)
	which does not need any normalization constant
	*/

	for(i=0;i<N/2+1;i++)
		u0[i] = u0[i] * Phat[i];
	for(i=1;i<(int)N/2;i++)
		u0[N-i] = u0[N-i] * Phat[i];

	fftw_execute(plan_ub);

	// Save final array to a text file
	sprintf(name2,"FGF1D_Out.dat");
  if( (f2 = fopen(name2,"w")) == NULL)
    errorwc(name2);

	for(i=0;i<N;i++)
    fprintf(f2, "%le \n",u0[i]);

	CLOSEFILE(f2);

  fftw_destroy_plan(plan_ub);
  fftw_destroy_plan(plan_uf);
  fftw_free(u0);
	FREEP(Phat)

	printf("\n\n And we are done \n\n\a");

return 0;
}
