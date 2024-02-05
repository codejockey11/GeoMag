#include "CCoordGeodetic.h"

CCoordGeodetic::CCoordGeodetic()
{
	memset(this, 0x00, sizeof(CCoordGeodetic));
}

CCoordGeodetic::CCoordGeodetic(CGeoid* Geoid, double lat, double lon, double alt)
{
	memset(this, 0x00, sizeof(CCoordGeodetic));

    lambda = lon;
    phi = lat;
    HeightAboveGeoid = alt;

    /* This converts the height above mean sea level to height above the WGS-84 ellipsoid */
    if (Geoid->UseGeoid == 1)
    {
        double DeltaHeight = Geoid->GetGeoidHeight(lat, lon);

        HeightAboveEllipsoid = HeightAboveGeoid + DeltaHeight / 1000;
    }
    else
    {
        HeightAboveEllipsoid = HeightAboveGeoid;
    }

}

CCoordGeodetic::~CCoordGeodetic()
{

}