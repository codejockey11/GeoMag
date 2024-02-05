#ifndef _CMagneticResults
#define _CMagneticResults

#include "standard.h"

#include "CLegendreFunction.h"
#include "CMagneticModel.h"
#include "CSphericalHarmonicVariables.h"
#include "CCoordSpherical.h"

class CMagneticResults
{
public:

	double Bx; /* North */
	double By; /* East */
	double Bz; /* Down */

	CMagneticResults();

	~CMagneticResults();

	void Summation(CLegendreFunction* LegendreFunction, CMagneticModel* MagneticModel, CSphericalHarmonicVariables* SphVariables, CCoordSpherical* CoordSpherical);

	void SummationSpecial(CMagneticModel* MagneticModel, CSphericalHarmonicVariables* SphVariables, CCoordSpherical* CoordSpherical);

	void SecVarSummation(CLegendreFunction* LegendreFunction, CMagneticModel* MagneticModel, CSphericalHarmonicVariables* SphVariables, CCoordSpherical* CoordSpherical);

	void SecVarSummationSpecial(CMagneticModel* MagneticModel, CSphericalHarmonicVariables* SphVariables, CCoordSpherical* CoordSpherical);

	void RotateMagneticVector(CCoordSpherical* CoordSpherical, CCoordGeodetic* CoordGeodetic, CMagneticResults* MagneticResults);

	void SphericalToCartesian(CCoordSpherical* CoordSpherical, double* x, double* y, double* z);
};

#endif