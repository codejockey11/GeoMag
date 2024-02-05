#ifndef _CSphericalHarmonicVariables
#define _CSphericalHarmonicVariables

#include "standard.h"

#include "CEllipsoid.h"
#include "CCoordSpherical.h"

class CSphericalHarmonicVariables
{
public:

	int nMax;

	double* RelativeRadiusPower;	/* [earth_reference_radius_km / sph. radius ]^n */
	double* cos_mlambda;			/* cp(m)  - cosine of (m*spherical coord. longitude) */
	double* sin_mlambda;			/* sp(m)  - sine of (m*spherical coord. longitude) */

	CSphericalHarmonicVariables();

	CSphericalHarmonicVariables(int n);

	~CSphericalHarmonicVariables();

	void ComputeSphericalHarmonicVariables(CEllipsoid* Ellip, CCoordSpherical* CoordSpherical, int nMax);
};

#endif