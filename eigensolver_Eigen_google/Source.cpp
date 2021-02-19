#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

# define pai 3.1415926
#include <iostream>
#include <cmath>
using namespace Eigen;
using namespace std;

int main() {
	int i, j, N;
	double U, a, L, V0;
	N = 100;              
	a = 1.0;
	V0 = 1.0;
	L = 50.0;
	MatrixXd H(N, N);
	MatrixXd S(N, N);
	VectorXd V(N);
	
	// initialize U, Matrix H, S
	U = 10;
	for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) { 
			if (i == 0 && j == 0) {
				S(i, j) = 1;
				H(i, j) = -V0 / L * a;
			}
			else {
				if (i == j) {				
					S(i, j) = 1;				
					H(i, j) = i * i * pai * pai / L / L - V0 / L *(a + sin(2.0 * i * pai / L * a) / 2.0 * i * pai / L);
				}			
				else {				
					S(i, j) = 0;			
					H(i, j) = -V0 / L *(sin((i - j) * pai / L * a) / ((i - j) * pai / L) + sin((i + j) * pai / L * a) / ((i + j) * pai / L));
				}	
			}			
		}
	}

	// solve for eigenvalues and return them as vector V
	EigenSolver<MatrixXd> es(H, false);
	V = es.eigenvalues().real();
	for (i = 0; i < N; i++) {
		if (V(i) < U) U = V(i);
	}
	
	cout << U << endl;
	cin.get();
}