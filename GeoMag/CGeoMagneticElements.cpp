#include "CGeoMagneticElements.h"

CGeoMagneticElements::CGeoMagneticElements()
{
	memset(this, 0x00, sizeof(CGeoMagneticElements));
}

CGeoMagneticElements::~CGeoMagneticElements()
{
	if (LegendreFunction)
	{
		delete LegendreFunction;
	}

	if (SphVariables)
	{
		delete SphVariables;
	}

	if (MagneticResultsSph)
	{
		delete MagneticResultsSph;
	}

	if (MagneticResultsSphVar)
	{
		delete MagneticResultsSphVar;
	}

	if (MagneticResultsGeo)
	{
		delete MagneticResultsGeo;
	}
	
	if (MagneticResultsGeoVar)
	{
		delete MagneticResultsGeoVar;
	}
}

void CGeoMagneticElements::CalculateFieldElements(CEllipsoid* Ellip, CCoordSpherical* CoordSpherical, CCoordGeodetic* CoordGeodetic, CMagneticModel* TimedMagneticModel)
{
	NumTerms = ((TimedMagneticModel->nMax + 1) * (TimedMagneticModel->nMax + 2) / 2);

	LegendreFunction = new CLegendreFunction(NumTerms);

	SphVariables = new CSphericalHarmonicVariables(TimedMagneticModel->nMax);

	SphVariables->ComputeSphericalHarmonicVariables(Ellip, CoordSpherical, TimedMagneticModel->nMax);

	LegendreFunction->AssociatedLegendreFunction(CoordSpherical, TimedMagneticModel->nMax);

	MagneticResultsSph = new CMagneticResults();

	MagneticResultsSph->Summation(LegendreFunction, TimedMagneticModel, SphVariables, CoordSpherical);

	MagneticResultsSphVar = new CMagneticResults();

	MagneticResultsSphVar->SecVarSummation(LegendreFunction, TimedMagneticModel, SphVariables, CoordSpherical);

	MagneticResultsGeo = new CMagneticResults();

	MagneticResultsGeo->RotateMagneticVector(CoordSpherical, CoordGeodetic, MagneticResultsSph);

	MagneticResultsGeoVar = new CMagneticResults();

	MagneticResultsGeoVar->RotateMagneticVector(CoordSpherical, CoordGeodetic, MagneticResultsSphVar);


	CGeoMagneticElements::CalculateGeoMagneticElements(MagneticResultsGeo);

	CGeoMagneticElements::CalculateSecularVariationElements(MagneticResultsGeoVar);
}

void CGeoMagneticElements::CalculateGeoMagneticElements(CMagneticResults* MagneticResultsGeo)
{
	X = MagneticResultsGeo->Bx;
	Y = MagneticResultsGeo->By;
	Z = MagneticResultsGeo->Bz;

	H = sqrt(MagneticResultsGeo->Bx * MagneticResultsGeo->Bx + MagneticResultsGeo->By * MagneticResultsGeo->By);
	
	F = sqrt(H * H + MagneticResultsGeo->Bz * MagneticResultsGeo->Bz);
	
	Decl = RAD2DEG(atan2(Y, X));
	Incl = RAD2DEG(atan2(Z, H));
}

void CGeoMagneticElements::CalculateSecularVariationElements(CMagneticResults* MagneticVariation)
{
	Xdot = MagneticVariation->Bx;
	Ydot = MagneticVariation->By;
	Zdot = MagneticVariation->Bz;
	Hdot = (X * Xdot + Y * Ydot) / H; /* See equation 19 in the WMM technical report */
	Fdot = (X * Xdot + Y * Ydot + Z * Zdot) / F;
	
	Decldot = 180.0 / M_PI * (X * Ydot - Y * Xdot) / (H * H);
	Incldot = 180.0 / M_PI * (H * Zdot - Z * Hdot) / (F * F);
	
	GVdot = Decldot;
}