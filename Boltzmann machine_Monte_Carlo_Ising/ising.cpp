#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <list> 

using namespace std;

// Global variables
                
int Lx, Ly, N;// size and number of sites of system 
int **spinmesh; // spin sites matrix                                                
double J = 1, T, Ttarget=4.0,BoltzmannFactors[17];
double eacu, macu, m2ac, *corree, *corrmm;                
int ncorrt = 10, Config;                 
list<double> elist, mlist;
void ComputeBoltzmannFactors() {
	for (int i = -8; i <= 8; i += 4) 
	{
		BoltzmannFactors[i + 8] = exp(-(i * J ) / T);
	}
}
// generate ramdom number from 0 to 1
inline double std_rand()
{
	return rand() / (RAND_MAX + 1.0);
}

  

// randomly choose +1/-1 for each site
void Initialize() {
	int i, j;
	spinmesh = new int*[Lx];
	for (i = 0; i < Lx; i++)
		spinmesh[i] = new int[Ly];
	for (i = 0; i < Lx; i++)
		for (j = 0; j < Ly; j++)
		{
		if (std_rand() < 0.5) spinmesh[i][j] = 1;
		else spinmesh[i][j] = -1;
		}			

}

void Initiallist() {

	eacu = macu = m2ac = 0;
	elist.clear();
	mlist.clear();
	if (corree != NULL) delete[] corree;
	if (corrmm != NULL) delete[] corrmm;
	corree = new double[ncorrt + 1];
	corrmm = new double[ncorrt + 1];
	for (int i = 0; i <= ncorrt; i++) 
	{
		corree[i] = 0.0;
		corrmm[i] = 0.0;
	}
	Config = 0;
}
/*
void MetropolisStep() {
	double rate;	
	int i, j, ileft, iright, jdown, jup, deltaE;
	i = int(Lx*std_rand());
	j = int(Ly*std_rand());
	if(i == 0) ileft = Lx - 1;
	else ileft = i - 1;
	
	if(i == Lx - 1) iright = 0;
	else iright = i + 1;

	if(j == 0) jdown = Ly - 1;
	else jdown = j - 1;

	if(j == Ly - 1) jup = 0;
	else jup = j + 1;
	deltaE = 2 * spinmesh[i][j] * (spinmesh[ileft][j] + spinmesh[iright][j] + spinmesh[i][jdown] + spinmesh[i][jup]);
	rate = BoltzmannFactors[deltaE + 8];
	if (std_rand() < ratio) 
	{
		spinmesh[i][j] = -spinmesh[i][j];
	}	
}
*/
void MCScreen() 
{
	for (int dice = 0; dice < N; dice++)
	{
	double rate;	
	int i, j, ileft, iright, jdown, jup, deltaE;
	i = int(Lx*std_rand());
	j = int(Ly*std_rand());
	if(i == 0) ileft = Lx - 1;
	else ileft = i - 1;
	
	if(i == Lx - 1) iright = 0;
	else iright = i + 1;

	if(j == 0) jdown = Ly - 1;
	else jdown = j - 1;

	if(j == Ly - 1) jup = 0;
	else jup = j + 1;
	deltaE = 2 * spinmesh[i][j] * (spinmesh[ileft][j] + spinmesh[iright][j] + spinmesh[i][jdown] + spinmesh[i][jup]);
	rate = BoltzmannFactors[deltaE + 8];
	if (std_rand() < rate) spinmesh[i][j] = -spinmesh[i][j];	
	}
}

double Magnetization() {
	int szsum = 0;
	double szpersite=0.0;
	for (int i = 0; i < Lx; i++)
	{
		for (int j = 0; j < Ly; j++) 
		{
			szsum += spinmesh[i][j];
		}
	}
	szpersite=szsum / double(N);
	return szpersite;
}

double Energy() {
	int Esum = 0, iright, jup;
	double Epersite=0.0;
	for (int i = 0; i < Lx; i++)
	{
		for (int j = 0; j < Ly; j++) 
		{
			// periodic boundary condition
			if(i == Lx - 1) iright = 0;
			else iright = i + 1;
			
			if(j == Ly - 1) jup = 0;
			else jup = j + 1;

			Esum += spinmesh[i][j] * (spinmesh[iright][j] + spinmesh[i][jup]);
		}
	}
	Epersite=(-J*Esum) / double(N);
	return Epersite;
}

double M2() {
	int E2sum = 0;
	double sz2persite=0.0;
	for (int i = 0; i < Lx; i++)
	{
		for (int j = 0; j < Ly; j++) 
		{
			E2sum += spinmesh[i][j] * spinmesh[i][j];
		}
	}
	sz2persite=E2sum / double(N);
	return sz2persite;
}

void scrambling() {
	double e = Energy();
	double m = Magnetization();
	double m2 = M2();	
	if (elist.size() == ncorrt) 
	{   	
		eacu += e;
		macu += m;
		m2ac += m2;
		corree[0] += e * e;
		corrmm[0] += m * m;
		Config++;
		list<double>::const_iterator ie = elist.begin(), im = mlist.begin();
		for (int i = 1; i <= ncorrt; i++) 
		{
			//corree[i] += *ie++ * e;
			corree[i] =corree[i] + *ie++ * e;
			//corrmm[i] += *im++ * m;
		    corrmm[i] =corrmm[i] + *im++ * m;
		}
		elist.pop_back();
		mlist.pop_back();
	}
	elist.push_front(e);
	mlist.push_front(m);
}

double tau_e, tau_m, tau_m2, evee, evem, evem2, stde, stdm;
/*
void Compute() {
	// energy correlation
	double av = eacu / Config;
	double c0 = corree[0] / Config - av * av;
	tau_e = 0;
	for (int i = 1; i <= ncorrt; i++)	tau_e += (corree[i] / Config - av * av) / c0;
	// magnetization correlation
	av = macu / Config;
	c0 = corrmm[0] / Config - av * av;
	tau_m = 0;
	for (int i = 1; i <= ncorrt; i++)	tau_m += (corrmm[i] / Config - av * av) / c0;
	// m^2 correlation
	av = corrmm[0] / Config;
	c0 = corrmm[0] / Config - av * av;
	tau_m2 = 0;
	for (int i = 1; i <= ncorrt; i++)	tau_m2 += (corrmm[i] / Config - av * av) / c0;
	// everage e, m, m^2
	evee = eacu / Config;
	evem = macu / Config;
	evem2 = m2ac / Config;
	// standard devation stde, stdm
	stde = sqrt(fabs(eacu * eacu / Config / Config - corree[0] / Config));
	stdm = sqrt(fabs(macu * macu / Config / Config - corrmm[0] / Config));
}
*/
int main(int argc, char *argv[]) {

	int trialsteps;
	Lx = 40;
	Ly = Lx;
	N = Lx * Ly;
	// Temperature from 0.5 to 4.0 with steps 0.1 
	double T1 = 0.5, T2 = 4.0;
	// trial steps 10000, MC steps 10000000
	int nT = 35, MCSweep = 10000000;
	// output M, M^2, E function of T 
	ofstream file("isingT.txt");
	// output E, M per MC sweep at T = 1.6
	ofstream file1("isingEM.txt");
	Initialize();	
	trialsteps = int(0.01 * MCSweep);
	for (int i = 0; i <= nT; i++) 
	{
		T = T1 + i * (T2 - T1) / double(nT);
		ComputeBoltzmannFactors();
		for (int i = 0; i < trialsteps; i++) {
			MCScreen();
//			if (T == Ttarget) file1 << Energy() << " " << Magnetization() << endl;
		}				
		Initiallist();		
		for (int i = 0; i < MCSweep; i++) {
			MCScreen(); // MC updated to a configuration of spin system with at time t 
			if (T == Ttarget) file1 << Energy() << " " << Magnetization() << endl;			
			scrambling();		
		}
		cout <<  T << "completed" <<endl ;
//		Compute();
	// energy correlation
	    double av = eacu / Config;
		double c0 = corree[0] / Config - av * av;
		tau_e = 0;
		for (int i = 1; i <= ncorrt; i++)	tau_e += (corree[i] / Config - av * av) / c0;
	// magnetization correlation
		av = macu / Config;
		c0 = corrmm[0] / Config - av * av;
		tau_m = 0;
		for (int i = 1; i <= ncorrt; i++)	tau_m += (corrmm[i] / Config - av * av) / c0;
	// m^2 correlation
		av = corrmm[0] / Config;
		c0 = corrmm[0] / Config - av * av;
		tau_m2 = 0;
		for (int i = 1; i <= ncorrt; i++)	tau_m2 += (corrmm[i] / Config - av * av) / c0;
	// everage e, m, m^2
		evee = eacu / Config;
		evem = macu / Config;
		evem2 = m2ac / Config;
	// standard devation stde, stdm
		stde = sqrt(fabs(eacu * eacu / Config / Config - corree[0] / Config));
		stdm = sqrt(fabs(macu * macu / Config / Config - corrmm[0] / Config));
		file << T << ' ' << tau_e << ' ' << tau_m << ' ' << tau_m2 << ' ' << evee << ' ' << evem << ' ' << evem2 << ' ' << stde << ' ' << stdm << endl;
	}
	file.close();
	return 0;
}
