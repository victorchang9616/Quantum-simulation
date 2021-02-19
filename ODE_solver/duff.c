#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#define NRANSI
#include "nr.h"
#include "nrutil.h"                                                        // define parameters of duffing equation
#define N 2
#define F0 2.7
#define a 0.7
#define b 0.6
#define omega 2.7
#define gamma 0.2
#define pi 3.1415

  
int kmax,kount;
int nrhs; 
float dxsav,*xp,**yp;                                                     //count times

void derative(float x,float y[],float dydx[])      
{
	nrhs++;
                                                              //dervative functions of duffing equation
	dydx[1] = y[2];
//	dydx[2] = 1.5 * y[1]- 2.5 * y[1] * y[1] * y[1] -gamma * y[2] + F0 * cos(omega * x);
	dydx[2] = -gamma * y[2] + 2 * a * y[1] - 4 * b * pow(y[1],3.0)  + F0 * cos(omega * x);
}

int main(void)
{
	int n,nb,nk;
	float eps=1.0e-3,h1=0.1,hmin=0.0;
	float *ystart,x1,x2;
	FILE *file = fopen("duff.txt", "w");                     //open file
	kmax=1500;                                              //number of steps
	dxsav=(x2-x1)/kmax;
	ystart=vector(1,N);
	xp=vector(1,kmax);
	yp=matrix(1,N,1,kmax);
	                                                        //initial values of x and p
	ystart[1]=0.0;
	ystart[2]=0.0;
	for(n=0;n<8000;n++){                                    //repeat 8000 times(8000 datas)
	nrhs = 0;
	x1 = 2 * pi * n / omega;                                // n = # of period
	x2 = 2 * pi * (n+1) / omega;
	odeint(ystart,N,x1,x2,eps,h1,hmin,&nk,&nb,derative,rkqs);
	fprintf(file, "%f %f\n",yp[1][kount],yp[2][kount]);      //output data
	ystart[1]=yp[1][kount];
	ystart[2]=yp[2][kount];
	}
	fclose(file);                                            //close file
	free_matrix(yp,1,N,1,kmax);
	free_vector(xp,1,kmax);
	free_vector(ystart,1,N);
	return 0;
}
#undef NRANSI
