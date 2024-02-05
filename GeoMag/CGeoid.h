#ifndef _CGeoid
#define _CGeoid

#include "standard.h"

class CGeoid
{
public:

    int NumbGeoidCols; /* 360 degrees of longitude at 15 minute spacing */
    int NumbGeoidRows; /* 180 degrees of latitude  at 15 minute spacing */
    int NumbHeaderItems; /* min, max lat, min, max long, lat, long spacing*/
    
    int ScaleFactor; /* 4 grid cells per degree at 15 minute spacing  */
    
    // included in EGM9615.h
    //float* GeoidHeightBuffer;
    
    int NumbGeoidElevs;
    int Geoid_Initialized; /* indicates successful initialization */
    int UseGeoid; /* Is the Geoid being used? */

	CGeoid();

    CGeoid(const char p);

	~CGeoid();

    double GetGeoidHeight(double lon, double lat);

    void ConvertGeoidToEllipsoidHeight();

private:

    long Index;

    double DeltaX, DeltaY;
    double ElevationSE, ElevationSW, ElevationNE, ElevationNW;
    double OffsetX, OffsetY;
    double PostX, PostY;
    double UpperY, LowerY;
    double DeltaHeight;
};

#endif