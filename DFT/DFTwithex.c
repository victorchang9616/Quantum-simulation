#include <stdio.h>
#include <math.h>
#define NRANSI
#include "nr.h"
#include "nrutil.h"
#define rmax 5.0   // the rmax 
#define N 2200 // the total points from 0 to rmax
#define tol 1.0e-4
#define pi 3.141592


float dxsav, *xp, *xb1, *xb2, **yp;
int kmax, kount;
int nrhs;                                                     
float delr = rmax / (float)N;      // 0.01        
float *VH,  *u, E,*Vx;


void 
poisson(float x, float y[], float dydx[]) 
{
	nrhs++;
	dydx[1] = y[2];
	int i = x / delr;
	dydx[2] = -pow(u[i], 2) / x;
}
//Normalize u
void 
normu()
{
	float integral=0;
	int i=0;
	for(i=1; i<N; i++){
		integral+=delr*u[i]*u[i];
	}
	for(i=1; i<N; i++){
		u[i]/=sqrt(integral);
	}
}

// subroutine hartreecal uses u to calculate VH[]
void 
hartreecal() 
{
	int  nbad, nok;
	float eps = 1.0e-4, h1 = 0.1, hmin = 0.0, x1 = 0.001, *ystart;

	kmax = 100;
	dxsav = delr / kmax;
	ystart = vector(1, 2);
	xp = vector(1, kmax);
	yp = matrix(1, 2, 1, kmax);

	ystart[1] = 0;
	ystart[2] = 0.1;    

	for (int i = 0; i < N; i++) 
	{
		odeint(ystart, 2, x1 + i*delr, x1 + (i + 1)*delr, eps, h1, hmin, &nok, &nbad, poisson, rkqs);
		VH[i + 2] = yp[1][kount];
		ystart[1] = yp[1][kount];
		ystart[2] = yp[2][kount];
	}
	VH[1] = VH[2];
	VH[0] = VH[2];
	float qmax = 0.0;
	for (int i = 0; i < N; i++) 
	{
		qmax += pow(u[i],2.0) * delr;
	}

	float alpha = (qmax - VH[N]) / rmax;
	for (int i = 0; i <= N; i++) 
	{
		VH[i] = 2.0 * (alpha*i*delr + VH[i] )/ (i * delr); // 2?
	}
	free_matrix(yp, 1, 2, 1, kmax);
	free_vector(xp, 1, kmax);
	free_vector(ystart, 1, 2);

}

// function exchangecal uses u to calculate Vx with local density approx. 
void 
exchangecal() 
{
	
	for (int i = 1; i <= N; i++) 
	{		
		Vx[i] = -pow(3 * u[i] * u[i] / (2 * pow(pi * delr * i, 2) ), (1.0 / 3.0));
	}
}

// function u0 requires E, VH[], Vx[] to calculate u[], returns u[0]
float 
u0(float E0) 
{                                  
	int nbad, nok;
	float eps = 1.0e-4, h1 = 0.1, hmin = 0.0, *ystart;
	E = E0;
	kmax = 100;                                         
	dxsav = delr / kmax;

	ystart = vector(1, 2);
	xp = vector(1, kmax);
	yp = matrix(1, 2, 1, kmax);

	ystart[1] = 0;
	ystart[2] = 0.0012; 
	
	for (int i = 0; i<(N-1); i++) 
	{
		odeint(ystart, 2, rmax - i*delr, rmax - (i + 1)*delr, eps, h1, hmin, &nok, &nbad, derivs, rkqs);
		u[N - i - 1] = yp[1][kount];
		ystart[1] = yp[1][kount];
		ystart[2] = yp[2][kount];
	}
	u[N] = u[N - 1];
	u[0] = u[1];
	normu();
	free_matrix(yp, 1, 2, 1, kmax);
	free_vector(xp, 1, kmax);
	free_vector(ystart, 1, 2);

	return u[0];
}
void 
derivs(float x, float y[], float dydx[]) 
{
	nrhs=nrhs+1;
	dydx[1] = y[2];
	int i = x / delr;         
	dydx[2] = 2.0*(-E - 2.0 / x + VH[i] + Vx[i])*y[1];
}
int main()
{
	int nb;
	float zero, epsilon, Energy = 0.0, Energylocal;
	FILE *file = fopen("DFTwiexepi.txt", "w"); 
	FILE *file1 = fopen("DFTwiexhart.txt", "w"); 
	FILE *file2 = fopen("DFTwiexexch.txt", "w"); 
    xb1 = vector(1, 20);
	xb2 = vector(1, 20);
	u = vector(1, (N + 1));
	VH = vector(1, (N + 1));
	Vx = vector(1, (N + 1));
	// initialization
	E = -2.5;
	for (int i = 0; i < N + 1; i++) 
	{
		VH[i] = 0.0;
		Vx[i] = 0.0;
	}
	zbrak(u0, -5.0, 1.0, 30, xb1, xb2, &nb);      		
	epsilon = zbrent(u0, xb1[1], xb2[1], tol);
	printf("initial E=\n %10.6f\n", epsilon);
	zero = u0(epsilon);  
	E = epsilon;
	for (int i = 0; i < 250; i++) 
	{
		hartreecal();
		exchangecal();
		zbrak(u0, -5.0, 1.0, 30, xb1, xb2, &nb);                		
		epsilon = zbrent(u0, xb1[1], xb2[1], tol);
		zero = u0(epsilon);
		printf("for %dth ieration check zero =%10.6f\n", i, zero);
		fprintf(file, "%d %f\n",i,epsilon);
		printf("for %dth iteration of epsilon=%10.6f\n", i, epsilon);;
		Energylocal = 2.0 * epsilon;
		for (int i = 1; i <= N; i++) 
		{
			Energylocal -= delr * u[i] * u[i] * VH[i];
		}
		fprintf(file1, "%d %f\n",i,Energylocal);
		for (int i = 1; i <= N; i++) 
		{
			Energylocal += 0.5 * delr  *u[i] * u[i] * Vx[i];     
		}
		fprintf(file2, "%d %f\n",i,Energylocal);
		if (fabs(Energy - Energylocal) < tol) break;
		else 
		{
			E = epsilon;
			Energy = Energylocal;
		}
		printf("the Energy =%10.6f\n", Energy);		
		printf("\n");
	}
	printf("the result is: %10.6f\n" , Energylocal);
	fclose(file);
	fclose(file1);
	fclose(file2);
    return 0;
}

#undef NRANSI
