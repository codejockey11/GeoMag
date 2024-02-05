#ifndef _CEllipsoid
#define _CEllipsoid

#include "standard.h"

class CEllipsoid
{
public:

    double a; /*semi-major axis of the ellipsoid*/
    double b; /*semi-minor axis of the ellipsoid*/
    double fla; /* flattening */
    double epssq; /*first eccentricity squared */
    double eps; /* first eccentricity */
    double re; /* mean radius of  ellipsoid*/

	CEllipsoid();

	~CEllipsoid();
};

#endif