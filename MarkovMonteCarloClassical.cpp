#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <time.h>
#include <thread>
#include <mutex>
#include <vector>

using namespace std;

/*const double x0 = 0.0;
const double y0 = 0.0;
const double pX0 = 0.0;
const double pY0 = 0.0;
const double xMin = -4.0;
const double xMax = 4.0;
const double yMin = -4.0;
const double yMax = 4.0;
const double pXMin = -4.0;
const double pXMax = 4.0;
const double pYMin = -4.0;
const double pYMAx = 4.0;*/

const double delta = 1.0;
//const double T = 1.0;
const double K = 1.0;
const double m = 1.0;
const double omega = 1.0;
const int rangeOfTTests = 100;
const int iterations = 1000000;

mutex mtx;

ofstream oFile;
ofstream eFile;

double V(double x, double y);
void computeChainStepPos(double &x, double &y, double dx, double dy, double T);
void computeChainStepMom(double &p, double dP, double T);
void runSim(double T, int id, vector<double> &Envec, vector<double> &Ervec);




main() {
	
	vector<thread> threadList;
	vector<double> resultList(rangeOfTTests);
	vector<double> errorList(rangeOfTTests);
	for(int i=0; i < rangeOfTTests; i++) {
		thread th(runSim, (double)i/10.0, i, ref(resultList), ref(errorList));
		threadList.push_back(move(th));
	}

	for(int i=0; i < rangeOfTTests; i++) {
		threadList.at(i).join();
	}
	oFile.open("data.txt", ios::app);
	eFile.open("error.txt", ios::app);
	for(int i=0; i < rangeOfTTests; i++) {
		oFile << (double)i/10.0 << "	" << resultList.at(i) << endl;
		eFile << (double)i/10.0 << "	" << errorList.at(i) << endl;
	}
	oFile.close();
	eFile.close();
}

void runSim(double T, int id, vector<double> &Envec, vector<double> &Ervec) {
	//mtx.lock();
	/*mtx.lock();
	cout << "Started thread: " << T << endl;
	mtx.unlock();*/
	double x = 0.0;
	double y = 0.0;
	double px = 0.0;
	double py = 0.0;
	double E = 0.0;
	double Etot = 0.0;
	double ESqtot = 0.0;
	srand (time(NULL));
	double d;
	for(int i=0; i < iterations; i++) {
		d = ((double)(rand())/(double)(RAND_MAX)*2*delta) - delta;
		computeChainStepPos(x,y,d,0,T);

		d = ((double)(rand())/(double)(RAND_MAX)*2*delta) - delta;
		computeChainStepPos(x,y,0,d,T);

		d = ((double)(rand())/(double)(RAND_MAX)*2*delta) - delta;
		computeChainStepMom(px, d, T);

		d = ((double)(rand())/(double)(RAND_MAX)*2*delta) - delta;
		computeChainStepMom(py, d, T);

		//cout << (pow(px,2) + pow(py,2))/(2*m) + .5 * m * pow(omega,2) * (pow(x,2) + pow(y,2)) << endl;

		E = ((pow(px,2) + pow(py,2))/(2*m) + .5 * m * pow(omega,2) * (pow(x,2) + pow(y,2)))/1.0;
		Etot += E;
		ESqtot += E*E;




	}
	Envec.at(id) = Etot / (double)iterations;
	//mtx.lock();
	//cout << "ID: " << id << " <E^2>: " << ESqtot/ (double)iterations << " <E>: " << Etot/ (double)iterations << endl;
	//mtx.unlock();
	Ervec.at(id) = pow((ESqtot / (double)iterations - pow(Etot / (double)iterations,2)),.5)/pow(iterations,.5);

	//mtx.unlock();
}

double V(double x, double y) {
	double r =  pow(pow(x,2) + pow(y,2), .5);
	return(.5*m*pow(omega,2)*pow(r,2));

}

void computeChainStepPos(double &x, double &y, double dx, double dy, double T) {
	double pi = exp(-1*V(x+dx,y+dy)/(K*T))/exp(-1*V(x,y)/(K*T));
	if( pi > (double)(rand())/(double)(RAND_MAX)) {
		x = x + dx;
		y = y + dy;
	}
}

void computeChainStepMom(double &p, double dP, double T) {
	double pi = exp(-1*pow(p + dP, 2)/(2*m*K*T))/exp(-1*pow(p, 2)/(2*m*K*T));
	if( pi > (double)(rand())/(double)(RAND_MAX)) {
		p = p + dP;
	}

}