// DLA Simulationi with simple visualization of result 
// written by Kevin E. Bassler, 1/22/2019
// version 2: visuslizes each particle addition
// version 3: visuslizes each step
//
// Note: The lattice has two outer rows and columns on each side that are "dead." 
// If a walker walks onto a dead site, then it "dies."  
#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define NRANSI
#include "nr.h"
#include "nrutil.h"

int L, ll;
long seed;
int **lattice;			// a pointer to the pointer to lattice

void initialize_lattice();
void addparticle();
void write_lattice();
void write_latticewithwalker(int i, int j);

int main() {

	int n, npart;

	printf("The lattice will be L x L. How large is L? (Should be odd integer.): ");
	scanf("%d", &L);
	printf("Input number of particles added to the seed: ");
	scanf("%d", &npart);
	printf("Input random number seed. (Should be large negative odd integer.): ");
	scanf("%ld", &seed);

	ll = (L - 1) / 2;
	lattice = imatrix(-(ll + 2), ll + 2, -(ll + 2), ll + 2);	// declares lattice as matrix

	initialize_lattice();		// initializes empty lattice with seed in middle
	//write_lattice();		// outputs simple visualization
	//getchar();			// waits for return before proceeding
	for (n = 1; n <= npart; ++n) {
		addparticle();		// adds new particle to seed
		//write_lattice();
		//getchar();
	}
	write_lattice();

	FILE *fp; //save y,x,lattice[y][x] into txt file
	FILE *fp1; //save y,x if lattice[y][x]==1 into txt file for plotting
	fp = fopen("dla.txt", "w");
  fp1 = fopen("dla1.txt", "w");
	int x, y;
	for (y = -ll - 2; y <= ll + 2; ++y) {
		for (x = -ll - 2; x <= ll + 2; ++x) {
			fprintf(fp, "%d %d %d\n", y, x, lattice[y][x]);
      if(lattice[y][x]==1) fprintf(fp1, "%d %d \n", y, x);
		}
	}
	fclose(fp);
	fclose(fp1);

	free_imatrix(lattice, -(ll + 2), ll + 2, -(ll + 2), ll + 2);	// frees memory
	return 0;
}

void write_lattice()
{
	int x, y;

	for (y = -ll - 2; y <= ll + 2; ++y) {
		for (x = -ll - 2; x <= ll + 2; ++x) {
			if (lattice[y][x] == 0)  printf("- ");		// empty site represented by -
			else if (lattice[y][x] == 1)  printf("o ");	// seed represented by o
		}
		printf("\n");
	}
	printf("\n");

	return;
}

void write_latticewithwalker(int i, int j)
{
	int x, y;

	for (y = -ll - 2; y <= ll + 2; ++y) {
		for (x = -ll - 2; x <= ll + 2; ++x) {
			if ((x == i) && (y == j)) printf("+ ");			// walker represented by +
			else if (lattice[y][x] == 0)  printf("- ");
			else if (lattice[y][x] == 1)  printf("o ");
		}
		printf("\n");
	}
	printf("\n");

	return;
}

void addparticle()
{
	int x, y, roll, side, site, done, nsum;

	done = 0;

	do {
		roll = 4 * L*ran2(&seed); // random int to determine walker's initial location: (x,y)
		//printf("%d\n",roll);
		side = roll / L;
		site = roll % L;
		if (side == 0) { x = -ll + site; y = ll; }
		else if (side == 1) { x = ll; y = ll - site; }
		else if (side == 2) { x = ll - site; y = -ll; }
		else if (side == 3) { x = -ll; y = -ll + site; }

		//write_latticewithwalker(x, y);	// visualize lattice and walker
		//getchar();

		while ((done == 0) && (abs(x) <= ll) && (abs(y) <= ll))  // do unless done or dead   
		{
			roll = 4 * ran2(&seed);		// random int to determine walk direction
			if (roll == 0) --x;
			else if (roll == 1) ++x;
			else if (roll == 2) --y;
			else if (roll == 3) ++y;
			nsum = lattice[y - 1][x] + lattice[y + 1][x] + lattice[y][x - 1] + lattice[y][x + 1];
			if (nsum > 0) done = 1;		// determine if next to seed
			//write_latticewithwalker(x, y);
			//getchar();
		}

	} while (done == 0);	// repeat if not done (happens only if dead)

	lattice[y][x] = 1;		// add walker to seed

	return;
}

void initialize_lattice()
{
	int x, y;

	for (y = -ll - 2; y <= ll + 2; ++y) {
		for (x = -ll - 2; x <= ll + 2; ++x) {
			lattice[y][x] = 0;		// initialize empty lattice
		}
	}
	lattice[0][0] = 1;		// start seed as single particle at center of lattice

	return;
}

#undef NRANSI
