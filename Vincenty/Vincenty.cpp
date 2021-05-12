// Vincenty.cpp: Calcolare LAT e LONG a partire da un punto conosciuto (sempre espresso in lat e long) e la relativa distanza e azimuth fra i due punti
#include <stdio.h>
#include <math.h>
#include <iostream>

using namespace std;


//------------------------------------------------------------------

//Variables and constants:

const float pi = 3.1415926535897932384626433832795028841971693993751;
const int nauticalMile = 1852;
double LAT1 = 45.881433, LONG1 = 8.677776, A_12 = 10, RANGE = 10;
const float earthEquatorialRad = 6378137.0;
const float earthPolarRad = 6356752.3142;
const double  earthEccentricity = 6.73949674227e-3;
const double flattng = 0.0033528106647;

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

double toRad(double coordinates) {  //Returns a given degree angle transformed to Radians
  return (coordinates/180)*pi;
}

double nmToMeters(double nm) {  //Returns a given nautical miles value to meters.
	return nm * nauticalMile;
}






int main() {

	//Leggo da stdin la quaterna di dati che rappresentano tutto l' INPUT:

	//scanf_s("%lf %lf %lf %lf", &LAT1, &LONG1, &RANGE, &A_12);

	//Una volta che ricevo l'input inizio a trasformare le variabili da Deg a Rad, da NM a M

	double LAT1Rad = toRad(LAT1);
	double LONG1Rad = toRad(LONG1);
	double A12Rad = toRad(A_12);
	double meterRange = nmToMeters(RANGE);

	double LAT2 = 0;
	double LONG2 = 0;
	double A21 = 0;

	double U1 = atan((1 - flattng) * tan(LAT1Rad));
	double sigma1 = atan(tan(U1) / cos(A12Rad));
	double sinAlpha = pow(cos(U1) * sin(A12Rad), 2);
	double sqCosAlpha = 1 - sinAlpha;

	double sqU = sqCosAlpha * earthEccentricity;

	double ABig = 1 + (sqU / 16384) * (4096 + sqU * (-768 + sqU * (320 - 175 * sqU)));
	double BBig = (sqU / 1024) * (256 + sqU * (-128 + sqU * (74 - 47 * sqU)));

	double sigmaInit = meterRange / (earthPolarRad * ABig);
	double sigma = sigmaInit;
	double oldSigma = 1000;
	int cnt = 0;
	double cos2SigmaM = 0.0;



	while (abs(sigma - oldSigma) > 1e-12 and cnt < 150) {

		cnt += 1;

		cos2SigmaM = cos(2 * sigma1 + sigma);
		double deltaSigma = BBig * sin(sigma) * (cos2SigmaM + (BBig / 4) * (cos(sigma) * (-1 + 2 * pow(cos2SigmaM, 2)) - (BBig / 6) * cos2SigmaM * (-3 + 4 * pow(sin(sigma), 2) * (-3 + 4 * pow(cos2SigmaM, 2)))));

		oldSigma = sigma;
		sigma = sigmaInit + deltaSigma;

	}

	double denTanPhi2 = (pow((1 - flattng) * (pow(sinAlpha, 2) + pow(sin(U1) * sin(sigma) - cos(U1) * cos(sigma) * cos(A12Rad), 2)), 0.5));
	double numTanPhi2 = (sin(U1) * cos(sigma) + cos(U1) * sin(sigma) * cos(A12Rad));


	double LAT2Rad = 0.0;
	double LONG2Rad = 0.0;


	if ((denTanPhi2 == 0) && (numTanPhi2 > 0)) LAT2Rad = pi / 2;
	else if ((denTanPhi2 == 0) && (numTanPhi2 < 0))  LAT2Rad = -pi / 2;
	else {

		double tanphi2 = numTanPhi2 / denTanPhi2;

		double C = (flattng / 16) * (sqCosAlpha) * (4 + flattng * (4 - 3 * (sqCosAlpha)));

		double lambda = atan2((sin(sigma) * sin(A12Rad)), (cos(U1) * cos(sigma) - sin(U1) * sin(sigma) * cos(A12Rad)));

		double omega = lambda - (1 - C) * flattng * sinAlpha * (sigma + C * sin(sigma) * (cos2SigmaM + C * cos(sigma) * (-1 + 2 * pow(cos2SigmaM, 2))));

		LONG2Rad = LONG1Rad + omega;

		LAT2Rad = atan(tanphi2);
	}

	LAT2 = LAT2Rad * 180 / pi;

	LONG2 = LONG2Rad * 180 / pi;

	A21 = atan2(sinAlpha, (-sin(U1) * sin(sigma) + cos(U1) * cos(sigma) * cos(A12Rad))) * 180 / pi + 180;

	cout << LAT2 << " " << LONG2 << " " << A21 << endl;

	return 0;
}