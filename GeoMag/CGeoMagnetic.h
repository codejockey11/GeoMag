#ifndef _CGeoMagnetic
#define _CGeoMagnetic

#include "standard.h"

#include "CMagneticModel.h"
#include "CGeoid.h"
#include "CEllipsoid.h"
#include "CCoordGeodetic.h"
#include "CDate.h"
#include "CCoordSpherical.h"
#include "CGeoMagneticElements.h"

class CGeoMagnetic
{
public:

	CMagneticModel* MagneticModel;
	CMagneticModel* TimedMagneticModel;
	CGeoid* Geoid;
	CEllipsoid* Ellip;
	CCoordGeodetic* CoordGeodetic;
	CDate* UserDate;
	CCoordSpherical* CoordSpherical;
	CGeoMagneticElements* GeoMagneticElements;

	CGeoMagnetic();

	CGeoMagnetic(double sdate, const char projection);

	~CGeoMagnetic();

};

#endif