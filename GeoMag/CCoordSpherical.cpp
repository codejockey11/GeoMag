#include "CCoordSpherical.h"

CCoordSpherical::CCoordSpherical()
{
	memset(this, 0x00, sizeof(CCoordSpherical));
}

CCoordSpherical::~CCoordSpherical()
{

}

void CCoordSpherical::GeodeticToSperical(CEllipsoid* Ellip, CCoordGeodetic* CoordGeodetic)
{
/*
 ** Convert geodetic coordinates, (defined by the WGS-84
 ** reference ellipsoid), to Earth Centered Earth Fixed Cartesian
 ** coordinates, and then to spherical coordinates.
 */

    CosLat = cos(DEG2RAD(CoordGeodetic->phi));
    SinLat = sin(DEG2RAD(CoordGeodetic->phi));

    /* compute the local radius of curvature on the WGS-84 reference ellipsoid */

    rc = Ellip->a / sqrt(1.0 - Ellip->epssq * SinLat * SinLat);

    /* compute ECEF Cartesian coordinates of specified point (for longitude=0) */

    xp = (rc + CoordGeodetic->HeightAboveEllipsoid) * CosLat;
    zp = (rc * (1.0 - Ellip->epssq) + CoordGeodetic->HeightAboveEllipsoid) * SinLat;

    /* compute spherical radius and angle lambda and phi of specified point */

    r = sqrt(xp * xp + zp * zp);
    phig = RAD2DEG(asin(zp / r)); /* geocentric latitude */
    lambda = CoordGeodetic->lambda; /* longitude */
}