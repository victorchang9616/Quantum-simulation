#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <list> 

using namespace std;

// for ferromagnetic ising model J=1
double J = 1, T, w[17];                 
int Lx, Ly, N, **s;                                                           

void ComputeBoltzmannFactors() {
	for (int i = -8; i <= 8; i += 4) {
		w[i + 8] = exp(-(i * J ) / T);
	}
}
// generate ramdom number from 0 to 1
inline double std_rand()
{
	return rand() / (RAND_MAX + 1.0);
}

double eAv, mAv, m2Av, *cee, *cmm;                
int nSave = 10, Count, steps = 0;                 
list<double> eSave, mSave;  

// randomly choose +1/-1 for each cite
void Initialize() {
	int i, j;
	s = new int*[Lx];
	for (i = 0; i < Lx; i++)
		s[i] = new int[Ly];
	for (i = 0; i < Lx; i++)
		for (j = 0; j < Ly; j++){
		if (std_rand() < 0.5) s[i][j] = 1;
		else s[i][j] = -1;
		}			
	steps = 0;
}

void InitializeData() {
	int i;
	eAv = mAv = m2Av = 0;
	eSave.clear();
	mSave.clear();
	if (cee != NULL) delete[] cee;
	if (cmm != NULL) delete[] cmm;
	cee = new double[nSave + 1];
	cmm = new double[nSave + 1];
	for (i = 0; i <= nSave; i++)
		cee[i] = cmm[i] = 0;
	Count = 0;
}

void MetropolisStep() {
	int i, j, iPrev, iNext, jPrev, jNext, sum, deltaE;
	double ratio;
	i = int(Lx*std_rand());
	j = int(Ly*std_rand());

	// periodic boundary conditions
	if(i == 0) iPrev = Lx - 1;
	else iPrev = i - 1;
	
	if(i == Lx - 1) iNext = 0;
	else iNext = i + 1;

	if(j == 0) jPrev = Ly - 1;
	else jPrev = j - 1;

	if(j == Ly - 1) jNext = 0;
	else jNext = j + 1;

	sum = s[iPrev][j] + s[iNext][j] + s[i][jPrev] + s[i][jNext];
	deltaE = 2 * s[i][j] * sum;
	ratio = w[deltaE + 8];
	if (std_rand() < ratio) {
		s[i][j] = -s[i][j];
	}	
}

void MonteCarloSweep() {
	int i;
	for (i = 0; i < N; i++)	MetropolisStep();
	steps++;
}

double Magnetization() {
	int sum = 0, i, j;
	for (i = 0; i < Lx; i++){
		for (j = 0; j < Ly; j++) {
			sum += s[i][j];
		}
	}
	return sum / double(N);
}

double Energy() {
	int sum = 0, i, j, iNext, jNext;
	for (i = 0; i < Lx; i++){
		for (j = 0; j < Ly; j++) {
			// periodic boundary condition
			if(i == Lx - 1) iNext = 0;
			else iNext = i + 1;
			
			if(j == Ly - 1) jNext = 0;
			else jNext = j + 1;

			sum += s[i][j] * (s[iNext][j] + s[i][jNext]);
		}
	}
	return (-J*sum) / double(N);
}

double Magnetization2() {
	int sum = 0, i, j;
	for (i = 0; i < Lx; i++){
		for (j = 0; j < Ly; j++) {
			sum += s[i][j] * s[i][j];
		}
	}
	return sum / double(N);
}

void Accumulate() {
	int i;
	double e = Energy();
	double m = Magnetization();
	double m2 = Magnetization2();	
	if (eSave.size() == nSave) {   	
		eAv += e;
		mAv += m;
		m2Av += m2;
		cee[0] += e * e;
		cmm[0] += m * m;
		Count++;
		list<double>::const_iterator ie = eSave.begin(), im = mSave.begin();
		for (i = 1; i <= nSave; i++) {
			cee[i] += *ie++ * e;
			cmm[i] += *im++ * m;
		}
		eSave.pop_back();
		mSave.pop_back();
	}
	eSave.push_front(e);
	mSave.push_front(m);
}

double tau_e, tau_m, tau_m2, evee, evem, evem2, stde, stdm;

void Compute() {
	int i;
	// energy correlation
	double av = eAv / Count;
	double c0 = cee[0] / Count - av * av;
	tau_e = 0;
	for (i = 1; i <= nSave; i++)	tau_e += (cee[i] / Count - av * av) / c0;
	// magnetization correlation
	av = mAv / Count;
	c0 = cmm[0] / Count - av * av;
	tau_m = 0;
	for (i = 1; i <= nSave; i++)	tau_m += (cmm[i] / Count - av * av) / c0;
	// m^2 correlation
	av = cmm[0] / Count;
	c0 = cmm[0] / Count - av * av;
	tau_m2 = 0;
	for (i = 1; i <= nSave; i++)	tau_m2 += (cmm[i] / Count - av * av) / c0;
	// everage e, m, m^2
	evee = eAv / Count;
	evem = mAv / Count;
	evem2 = m2Av / Count;
	// standard devation stde, stdm
	stde = sqrt(fabs(eAv * eAv / Count / Count - cee[0] / Count));
	stdm = sqrt(fabs(mAv * mAv / Count / Count - cmm[0] / Count));
}

int main(int argc, char *argv[]) {

	int i, s, trialsteps;
	Lx = 20;
	Ly = Lx;
	N = Lx * Ly;
	// Temperature from 1.5 to 4.0 with steps 0.1 
	double T1 = 1.5, T2 = 4.0;
	// trial steps 10000, MC steps 10000000
	int TSteps = 25, MCSteps = 10000000;
	// output M, M^2, E function of T 
	ofstream file("isingT.txt");
	// output E, M per MC sweep at T = 1.6
	ofstream file1("isingEM.txt");
	Initialize();	
	trialsteps = int(0.001 * MCSteps);
	for (i = 0; i <= TSteps; i++) {
		T = T1 + i * (T2 - T1) / double(TSteps);
		ComputeBoltzmannFactors();
		for (s = 0; s < trialsteps; s++) {
			MonteCarloSweep();
			if (T == 1.6) file1 << Energy() << " " << Magnetization() << endl;
		}				
		InitializeData();		
		for (s = 0; s < MCSteps; s++) {
			MonteCarloSweep();
			Accumulate();		
		}
		cout << "completed" << ((T - 1.5) * 100.0 / 2.5) << "% " <<endl ;
		Compute();
		file << T << ' ' << tau_e << ' ' << tau_m << ' ' << tau_m2 << ' ' << evee << ' ' << evem << ' ' << evem2 << ' ' << stde << ' ' << stdm << endl;
	}
	file.close();
	return 0;
}