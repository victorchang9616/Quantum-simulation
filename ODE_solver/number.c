#include <stdio.h>
#include <math.h>
#define N 8000
int main(void)
{

	FILE *myFile= fopen("duff.txt", "r");
	FILE *outFile = fopen("number.txt", "w");
	float x[N], p[N];
	int i, j, k, l, n[128][128];
	float b, number;
	
	for (i = 0; i < N; i++)
	{
		fscanf(myFile, "%f", &x[i]);
		fscanf(myFile, "%f", &p[i]);
	}
	
	for (l = 2; l <= 128; l=l*2) 
	{
		number = 0;
       for (i = 0; i < l; i++) 
	   {
		for (j = 0; j < l; j++) n[i][j] = 0;
	   }
		b = 6.0 / l;
		for (int ni = 0; ni < N; ni++)	
		{
			for (j = 0; j < l; j++) 
			{
				for (k = 0; k < l; k++) 
				{
					if (x[ni] < ((j + 1) * b - 3.0) && x[ni] > (j * b - 3.0)) 
					{
						if (p[ni] < ((k + 1) * b - 3.0) && p[ni] > (k * b - 3.0)) n[j][k] = 1;
						
					}
				}
			}
		}
       for (i = 0; i < l; i++) 
	   {
		for (j = 0; j < l; j++) if (n[i][j] == 1) number++;
	   }
	  
	  fprintf(outFile, "%f %f\n", log(b), log(number));
	}
	

	fclose (myFile);
	fclose (outFile);
	return 0;
}
