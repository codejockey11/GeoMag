#include "CGeoid.h"

#include "EGM9615.h"

CGeoid::CGeoid()
{
	memset(this, 0x00, sizeof(CGeoid));
}

CGeoid::CGeoid(const char p)
{

    /* Sets EGM-96 model file parameters */

    NumbGeoidCols = 1441; /* 360 degrees of longitude at 15 minute spacing */
    NumbGeoidRows = 721; /* 180 degrees of latitude  at 15 minute spacing */
    NumbHeaderItems = 6; /* min, max lat, min, max long, lat, long spacing*/

    ScaleFactor = 4; /* 4 grid cells per degree at 15 minute spacing  */

    NumbGeoidElevs = NumbGeoidCols * NumbGeoidRows;

    Geoid_Initialized = 1; /* TODO: Geoid will be initialized only if this is set to zero */


    /* If needed modify height referencing */
    switch (p)
    {
    case 'E':
    {
        break;
    }
    case 'M':
    {
        UseGeoid = 1; /* height above MSL */

        break;
    }
    default:
    {
        printf("\nProjection should be E or M");

        return;
    }
    }

}

CGeoid::~CGeoid()
{

}

double CGeoid::GetGeoidHeight(double lat, double lon)
{
    /* Latitude out of range */
    if ((lat < -90) || (lat > 90))
    {
        return -2;
    }

    OffsetY = (90.0 - lat) * ScaleFactor;

    /* Longitude out of range */
    if ((lon < -180) || (lon > 360))
    {
        return -1;
    }

    if (lon < 0.0)
    {
        OffsetX = (lon + 360.0) * ScaleFactor;
    }
    else
    {
        OffsetX = lon * ScaleFactor;
    }

    /*  Find Four Nearest Geoid Height Cells for specified Latitude, Longitude;   */
    /*  Assumes that (0,0) of Geoid Height Array is at Northwest corner:          */

    PostX = floor(OffsetX);
    
    if ((PostX + 1) == NumbGeoidCols)
    {
        PostX--;
    }
    
    PostY = floor(OffsetY);
    
    if ((PostY + 1) == NumbGeoidRows)
    {
        PostY--;
    }

    Index = (long)(PostY * NumbGeoidCols + PostX);
    
    ElevationNW = (double)GeoidHeightBuffer[Index];
    ElevationNE = (double)GeoidHeightBuffer[Index + 1];

    Index = (long)((PostY + 1) * NumbGeoidCols + PostX);
    
    ElevationSW = (double)GeoidHeightBuffer[Index];
    ElevationSE = (double)GeoidHeightBuffer[Index + 1];

    /*  Perform Bi-Linear Interpolation to compute Height above Ellipsoid:        */

    DeltaX = OffsetX - PostX;
    DeltaY = OffsetY - PostY;

    UpperY = ElevationNW + DeltaX * (ElevationNE - ElevationNW);
    LowerY = ElevationSW + DeltaX * (ElevationSE - ElevationSW);

    DeltaHeight = UpperY + DeltaY * (LowerY - UpperY);

    return DeltaHeight;
}