#ifndef _CCoordSpherical
#define _CCoordSpherical

#include "standard.h"

#include "CCoordGeodetic.h"
#include "CEllipsoid.h"

class CCoordSpherical
{
public:

	double lambda;		/* longitude */
	double phig;		/* geocentric latitude */
	double r;			/* distance from the center of the ellipsoid */

	CCoordSpherical();

	~CCoordSpherical();

	void GeodeticToSperical(CEllipsoid* Ellip, CCoordGeodetic* CoordGeodetic);

private:

	double CosLat, SinLat, rc, xp, zp;
};

#endif