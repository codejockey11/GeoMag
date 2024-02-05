#ifndef _CCoordGeodetic
#define _CCoordGeodetic

#include "standard.h"

#include "CGeoid.h"

class CCoordGeodetic
{
public:

	double lambda; /* longitude */
	double phi; /* geodetic latitude */
	
	double HeightAboveEllipsoid; /* height above the ellipsoid (HaE) */
	double HeightAboveGeoid; /* (height above the EGM96 geoid model ) */
	
	int UseGeoid;

	CCoordGeodetic();

	CCoordGeodetic(CGeoid* Geoid, double lon, double lat, double alt);

	~CCoordGeodetic();
};

#endif