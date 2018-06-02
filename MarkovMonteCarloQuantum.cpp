#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <time.h>

using namespace std;



const double delta = 2.0;
const double K = 1.0;
const double m = 1.0;
const double omega = 1.0;
const double hbar = 1.0;
const int rangeOfTTests = 10;
const int iterations = 1000000;

ofstream oFile;
ofstream eFile;

void computeChainStep(int &n1, int &n2, int d1, int d2, double T);

main() {
	int n1 = 0;
	int n2 = 0;
	double Etot = 0.0;
	double ESqTot = 0.0;
	double E = 0.0;
	srand (time(NULL));
	double d;
	int intD;
	double T;
	oFile.open("data.txt", ios::app);
	eFile.open("error.txt", ios::app);
	for(double i=0; i < rangeOfTTests; i = i+.1) {
		T = i;
		n1 = 0;
		n2 = 0;
		Etot = 0.0;
		ESqTot = 0.0;
		for(int i=0; i < iterations; i++) {
			d = ((double)(rand())/(double)(RAND_MAX));
			if(d < (1.0/3.0)) {
				intD = -1;
			}
			else if(d < (2.0/3.0)){
				intD = 0;

			}
			else {
				intD = 1;
			}

			computeChainStep(n1, n2, intD, 0, T);

			d = ((double)(rand())/(double)(RAND_MAX));
			if(d < (1.0/3.0)) {
				intD = -1;
			}
			else if(d < (2.0/3.0)){
				intD = 0;

			}
			else {
				intD = 1;
			}

			computeChainStep(n1, n2, 0, intD, T);


			//cout << i << "	" << n << endl;

			E = (n1+.5)*hbar*omega + (n2+.5)*hbar*omega;
			Etot += E;
			ESqTot += E*E;



		}
		oFile << T << "	" << Etot / (double)iterations << endl;
		eFile << T << "	" << pow((ESqTot / (double)iterations - pow(Etot / (double)iterations,2)),.5)/pow(iterations,.5) << endl;

	}
	oFile.close();
	eFile.close();
}



void computeChainStep(int &n1, int &n2, int d1, int d2, double T) {
	double pi = exp(-1*((n1+d1+.5)*hbar*omega+(n2+d2+.5)*hbar*omega)/(K*T)) / exp(-1*((n1+.5)*hbar*omega+(n2+.5)*hbar*omega)/(K*T));
	if(n1+d1 >= 0 && n2+d2 >= 0 && pi > (double)(rand())/(double)(RAND_MAX)) {
		n1 = n1 + d1;
		n2 = n2 + d2;
		//cout << "accepted for T: " << T << endl;
	}

}