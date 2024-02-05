#include "CSphericalHarmonicVariables.h"

CSphericalHarmonicVariables::CSphericalHarmonicVariables()
{
	memset(this, 0x00, sizeof(CSphericalHarmonicVariables));
}

CSphericalHarmonicVariables::CSphericalHarmonicVariables(int n)
{
	memset(this, 0x00, sizeof(CSphericalHarmonicVariables));

	nMax = n;

	RelativeRadiusPower = (double*)calloc((nMax + 1), sizeof(double));
	
	cos_mlambda = (double*)calloc((nMax + 1), sizeof(double));
	
	sin_mlambda = (double*)calloc((nMax + 1), sizeof(double));
}

CSphericalHarmonicVariables::~CSphericalHarmonicVariables()
{
	if (RelativeRadiusPower)
	{
		free(RelativeRadiusPower);
	}
		
	if (cos_mlambda)
	{
		free(cos_mlambda);
	}
		
	if (sin_mlambda)
	{
		free(sin_mlambda);
	}
}

void CSphericalHarmonicVariables::ComputeSphericalHarmonicVariables(CEllipsoid* Ellip, CCoordSpherical* CoordSpherical, int nMax)
{
    double cos_lambda, sin_lambda;
    int m, n;
    
    cos_lambda = cos(DEG2RAD(CoordSpherical->lambda));
    sin_lambda = sin(DEG2RAD(CoordSpherical->lambda));
    
    /* for n = 0 ... model_order, compute (Radius of Earth / Spherical radius r)^(n+2)
    for n  1..nMax-1 (this is much faster than calling pow MAX_N+1 times).      */
    RelativeRadiusPower[0] = (Ellip->re / CoordSpherical->r) * (Ellip->re / CoordSpherical->r);
    
    for (n = 1; n <= nMax; n++)
    {
        RelativeRadiusPower[n] = RelativeRadiusPower[n - 1] * (Ellip->re / CoordSpherical->r);
    }

    /*
     Compute cos(m*lambda), sin(m*lambda) for m = 0 ... nMax
           cos(a + b) = cos(a)*cos(b) - sin(a)*sin(b)
           sin(a + b) = cos(a)*sin(b) + sin(a)*cos(b)
     */
    cos_mlambda[0] = 1.0;
    sin_mlambda[0] = 0.0;

    cos_mlambda[1] = cos_lambda;
    sin_mlambda[1] = sin_lambda;
    
    for (m = 2; m <= nMax; m++)
    {
        cos_mlambda[m] = cos_mlambda[m - 1] * cos_lambda - sin_mlambda[m - 1] * sin_lambda;
        sin_mlambda[m] = cos_mlambda[m - 1] * sin_lambda + sin_mlambda[m - 1] * cos_lambda;
    }
}