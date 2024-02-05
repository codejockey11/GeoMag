#ifndef _CLegendreFunction
#define _CLegendreFunction

#include "standard.h"

#include "CCoordSpherical.h"

class CLegendreFunction
{
public:

	int NumTerms;

	double* Pcup;	/* Legendre Function */
	double* dPcup;	/* Derivative of Legendre fcn */

	CLegendreFunction();

	CLegendreFunction(int n);

	~CLegendreFunction();

	void AssociatedLegendreFunction(CCoordSpherical* CoordSpherical, int nMax);

	void PcupLow(double x, int nMax);

	void PcupHigh(double x, int nMax);

};

#endif