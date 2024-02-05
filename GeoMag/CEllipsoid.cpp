#include "CEllipsoid.h"

CEllipsoid::CEllipsoid()
{
	memset(this, 0x00, sizeof(CEllipsoid));

    /* Sets WGS-84 parameters */

    a = 6378.137; /*semi-major axis of the ellipsoid in */
    b = 6356.7523142; /*semi-minor axis of the ellipsoid in */
    fla = 1 / 298.257223563; /* flattening */
    eps = sqrt(1 - (b * b) / (a * a)); /*first eccentricity */
    epssq = (eps * eps); /*first eccentricity squared */
    re = 6371.2; /* Earth's radius */
}

CEllipsoid::~CEllipsoid()
{

}