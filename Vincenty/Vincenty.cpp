// Vincenty.cpp: Calcolare LAT e LONG a partire da un punto conosciuto (sempre espresso in lat e long) e la relativa distanza e azimuth fra i due punti
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <vector>

using namespace std;


//------------------------------------------------------------------

//Variables and constants: should be read by a txt config file, like ros parameter server

const double pi = 3.1415927;
const int nauticalMile = 1852;
const float earthEquatorialRad = 6378137.0;
const double earthPolarRad = 6356752.3142;
const double  earthEccentricity = 6.73949674227e-3;
const double flattng = 0.0033528106647;


double LAT1 = 0, LONG1 = 0, A_12 = 0, RANGE = 0;
double LAT2 = 0, LONG2 = 0, A21 = 0;

/*
input:
  LAT1: Known Latitude (double type)
  LONG1: Known Longitude (double type)
  RANGE: Distanza in Miglia Nautiche fra P1-2 (float type)
  A12: Azimuth Relativo fra P1-2

Output:
  LAT2: Lat calcolata del P2
  LONG2: Long del P2
  A21: Azimuth che collega 2 a 1.

*/


 //----------------------------------------------------------------------------------------

// Functions Definitions:






inline char chooseVincenty() {
	char tmp = NULL;
	cout << "Hi! You wanna solve direct (D) or inverse (I) Vincenty?" << endl;
	cin >> tmp;
	fflush(stdin);
	if (tmp == 'D' || tmp == 'd' || tmp == 'I' || tmp == 'i') {
		return tmp;		
	}
	else {
		cout << "Wrong value! Please insert again!" << endl;
		chooseVincenty();
	}
}
vector <double>  getCoordinate(char &mode) {
	
	double value = 0.0;
	vector <double> tmpCoordinate;
	
	if (mode == 'D' || mode == 'd') {
		cout << "Please insert initial Latitude" << endl;
		cin >> value;
		tmpCoordinate.push_back(value);
		fflush(stdin);
		cout << "Please insert initial Longitude" << endl;
		cin >> value;
		tmpCoordinate.push_back(value);
		fflush(stdin);
		cout << "Please insert initial Azimuth" << endl;
		cin >> value;
		tmpCoordinate.push_back(value);
		fflush(stdin);
		cout << "Please insert initial Distance" << endl;
		cin >> value;
		tmpCoordinate.push_back(value);
		fflush(stdin);
	}else {

		cout << "Please insert Point 1 Latitude" << endl;
		cin >> value;
		tmpCoordinate.push_back(value);
		fflush(stdin);
		cout << "Please insert Point 1 Longitude" << endl;
		cin >> value;
		tmpCoordinate.push_back(value);
		fflush(stdin);
		cout << "Please insert Point 2 Latitude" << endl;
		cin >> value;
		tmpCoordinate.push_back(value);
		fflush(stdin);
		cout << "Please insert Point 2 Longitude" << endl;
		cin >> value;
		tmpCoordinate.push_back(value);
		fflush(stdin);

	}

	return tmpCoordinate;
}
bool dataVerifier(vector <double> &Coord, char& mode) {
	if (mode == 'D' || mode == 'd') cout << "These are saved data: " << endl << "LATITUDE: " << Coord[0] << endl << "LONGITUDE: " << Coord[1] << endl << "AZIMUTH: " << Coord[2] << endl << "RANGE: " << Coord[3] << endl;
	else cout << "These are saved data: " << endl << "LAT1: " << Coord[0] << endl << "LONG1: " << Coord[1] << endl << "LAT2: " << Coord[2] << endl << "LONG2: " << Coord[3] << endl;

	cout << "Do you confirm? Y/N" << endl;
	char confirm = NULL;
	cin >> confirm;
	fflush(stdin);

	if (confirm == 'N' || confirm == 'n') return false;
	else if (confirm == 'Y' || confirm == 'y') return true;
	else dataVerifier(Coord, mode);
}

//common Math Function
double toRad(double coordinates) {  //Returns a given degree angle transformed to Radians
	return (coordinates / 180) * pi;
}
double nmToMeters(double nm) {  //Returns a given nautical miles value to meters.
	return nm * nauticalMile;
}


double Afunct(double u){ return (1 + (u / 16384) * (4096 + u * (-768 + u * (320 - 175 * u))));}
double Bfunct(double u){ return ((u / 1024) * (256 + u * (-128 + u * (74 - 47 * u))));}
double Cfunct(double cos){ return (flattng / 16) * (cos) * (4 + flattng * (4 - 3 * (cos)));}
double deltaSigma(double &BBig, double &sigma, double &cos2SigmaM) { return (BBig * sin(sigma) * (cos2SigmaM + (BBig / 4) * (cos(sigma) * (-1 + 2 * pow(cos2SigmaM, 2)) - (BBig / 6) * cos2SigmaM * (-3 + 4 * pow(sin(sigma), 2) * (-3 + 4 * pow(cos2SigmaM, 2)))))); }
double dL(double& bigC, double& sinA, double& S, double& cos2S) { return (1 - bigC) * flattng * sinA * (S + bigC * sin(S) * (cos2S + bigC * cos(S) * (-1 + 2 * pow(cos2S, 2)))); }



void directVincenty(vector <double>& Coord) {
	//Una volta che ricevo l'input inizio a trasformare le variabili da Deg a Rad, da NM a M

	double LAT1Rad = toRad(Coord[0]); 
	double LONG1Rad = toRad(Coord[1]);
	double A12Rad = toRad(Coord[2]);
	double meterRange = nmToMeters(Coord[3]);

	double U1 = atan((1 - flattng) * tan(LAT1Rad));
	double sigma1 = atan2(tan(U1), cos(A12Rad));
	double sinAlpha = cos(U1) * sin(A12Rad);
	double sqCosAlpha = 1 - sinAlpha;

	double sqU = sqCosAlpha * earthEccentricity;

	double ABig = Afunct(sqU);
	double BBig = Bfunct(sqU);

	double sigmaInit = meterRange / (earthPolarRad * ABig);
	double sigma = sigmaInit;
	double oldSigma = 1000;
	int cnt = 0;
	double cos2SigmaM = 0.0;



	while (abs(sigma - oldSigma) > 1e-12 and cnt < 30) {

		cnt += 1;

		cos2SigmaM = cos(2 * sigma1 + sigma);
		double dSigma = deltaSigma(BBig, sigma, cos2SigmaM);

		oldSigma = sigma;
		sigma = sigmaInit + dSigma;

	}

	double denTanPhi2 = (1 - flattng) * pow(pow(sinAlpha, 2) + (pow(sin(U1) * sin(sigma) - cos(U1) * cos(sigma) * cos(A12Rad), 2)), 0.5);
	double numTanPhi2 = (sin(U1) * cos(sigma) + cos(U1) * sin(sigma) * cos(A12Rad));

	double LAT2Rad = 0.0;
	double LONG2Rad = 0.0;


	if ((denTanPhi2 == 0) && (numTanPhi2 > 0)) LAT2Rad = pi / 2;
	else if ((denTanPhi2 == 0) && (numTanPhi2 < 0))  LAT2Rad = -pi / 2;
	else {

		double tanphi2 = numTanPhi2 / denTanPhi2;
		double C = Cfunct(sqCosAlpha);  
		double lambda = atan2((sin(sigma) * sin(A12Rad)), (cos(U1) * cos(sigma) - sin(U1) * sin(sigma) * cos(A12Rad)));
		double omega = lambda - dL(C, sinAlpha, sigma, cos2SigmaM); 
		
		LONG2Rad = LONG1Rad + omega;
		LAT2Rad = atan(tanphi2);
	}



	LAT2 = LAT2Rad * 180 / pi;
	LONG2 = LONG2Rad * 180 / pi;
	A21 = atan2(sinAlpha, (-sin(U1) * sin(sigma) + cos(U1) * cos(sigma) * cos(A12Rad))) * 180 / pi + 180;

	cout << LAT2 << " " << LONG2 << " " << A21 << endl;
}

void inverseVincenty(vector <double>& Coord) {

	LAT1 = Coord[0];
	LONG1 = Coord[1];
	LAT2 = Coord[2];
	LONG2 = Coord[3];

	double LAT1Rad = toRad(LAT1);
	double LAT2Rad = toRad(LAT2);
	double LONG1Rad = toRad(LONG1);
	double LONG2Rad = toRad(LONG2);


	double U1 = atan(tan(LAT1Rad)*(1 - flattng));
	double U2 = atan(tan(LAT2Rad)*(1 - flattng));
	double L = LONG2Rad - LONG1Rad;

	double lambda = L;

	if(abs(LAT1Rad - LAT2Rad) < 1e-15 && abs(LONG2Rad - LONG1Rad) < 1e-15) return; //Verify that points are not overlapping
















}


int main() {	

	char dOrI = chooseVincenty();
	vector <double> Coord = getCoordinate(dOrI);
	bool coordOk = dataVerifier(Coord, dOrI);


	//check coordinates
	if (coordOk) {
		cout << "I begin computation" << endl;

		if (dOrI == 'D' or dOrI == 'd') directVincenty(Coord);
		else if (dOrI == 'I' or dOrI == 'i') inverseVincenty(Coord);
	}
	else getCoordinate(dOrI);	

	













































	
}

	



//TODO : Reverse Bearing e LAT2 sono COMPLETAMENTE sbagliati, la latitudine punta a un punto che in realtà sta a 3 volte le nautical miles della realtà, mentre bearing restituisce numeri senza senso
//(tipo azimuth =10 ritorna reverse bearing = 181).