#include "standard.h"

#include "CGeoMagnetic.h"
#include "CCoordGeodetic.h"

#define MODEL_RELEASE_DATE "04 Feb 2019"
#define VERSIONDATE_LARGE "$Date: 2019-02-04 10:40:43 -0700 (Mon, 04 Feb 2019) $"
#define VERSIONDATE_SMALL "2019-02-04"

CGeoMagnetic* GeoMagnetic = nullptr;
CCoordGeodetic* CoordGeodetic = nullptr;

FILE* infile = nullptr;
errno_t infileError;

double sdate;

char projection;

double alt;

double longitude;
double latitude;

char* GetArgcValue(int c, char* v[], const char* var)
{
    for (int i = 0; i < c; i++)
    {
        if (strcmp(v[i], var) == 0)
        {
            return v[i + 1];
        }
    }

    return nullptr;
}

int CheckForArgSwitch(int c, char* v[], const char* var)
{
    for (int i = 0; i < c; i++)
    {
        if (strcmp(v[i], var) == 0)
        {
            return i;
        }
    }

    return -1;
}

int main(int argc, char* argv[])
{
    printf("\n-----------------------------------------------");
    printf("\nWMM File processing program %s", VERSIONDATE_SMALL);
    printf("\n-----------------------------------------------");

    if (CheckForArgSwitch(argc, argv, "-h") > -1)
    {
        printf("\nWelcome to the World Magnetic Model (WMM) C-Program");
        printf("\nof the US National Geophysical Data Center");
        printf("\n               --- File Processing Program ----");
        printf("\n             --- Model Release Date: %s ---", MODEL_RELEASE_DATE);
        printf("\n           --- Software Release Date: %s ---", VERSIONDATE_SMALL);
        printf("\nUSAGE:");
        printf("\nFor example: WMM_file -f input_file output_file");
        printf("\nThis screen: WMM_file -h ");
        printf("\n");
        printf("\nThe input file may have any number of entries but they must follow");
        printf("\nthe following format");
        printf("\nDate and location Formats: ");
        printf("\n   Date: xxxx.xxx for decimal  (2015.15)");
        printf("\n   Altitude: M - Above mean sea level: E above WGS84 Ellipsoid ");
        printf("\n   Altitude: Kxxxxxx.xxx for kilometers  (K1000.13)");
        printf("\n             Mxxxxxx.xxx for meters  (m1389.24)");
        printf("\n             Fxxxxxx.xxx for feet  (F192133.73)");
        printf("\n   Lat/Lon: xxx.xxx in decimal  (-76.53)");
        printf("\n            (Lat and Lon must be specified in the same format.)");
        printf("\n   Date and altitude must fit model.");
        printf("\n   Lat: -90 to 90 (Use - to denote Southern latitude.)");
        printf("\n   Lon: -180 to 180 (Use - to denote Western longitude.)");
        printf("\n   Date: 1990.1f to 1995.1f");
        printf("\n   An example of an entry in input file");
        printf("\n   2015.15 E F30000 -70.3 -30.8 ");

        return 0;
    }

    if (CheckForArgSwitch(argc, argv, "-f") > -1)
    {
        printf("\n'f' switch: converting file with multiple locations.");
        printf("\n    The first five output columns repeat the input coordinates.");
        printf("\n    Then follows D, I, H, X, Y, Z, and F.");
        printf("\n    Finally the SV: dD, dI, dH, dX, dY, dZ, and dF");

        char* f = GetArgcValue(argc, argv, "-f");

        if (f == nullptr)
        {
            printf("\nInput file name not provided");

            return -1;
        }
        else
        {
            infileError = fopen_s(&infile, f, "rt");

            if (infileError)
            {
                printf("\nInput file not found:%s", f);

                return -1;
            }
        }
    }

    if (argc != 6)
    {
        printf("\nArguments not correctly provided");

        return -1;
    }

    printf("\nArguments as provided:");

    for (int i = 1; i < argc; i++)
    {
        printf("%s ", argv[i]);
    }

    
    
    // 2015.15 M F1200 39.342776 -81.440772
    sdate = atof(argv[1]);
    
    projection = argv[2][0];
    
    alt = atof(&argv[3][1]);

    /* Convert altitude to km */
    switch (argv[3][0])
    {
    case 'M':
    {
        alt *= 0.001;

        break;
    }
    case 'F':
    {
        alt /= 3280.0839895;

        break;
    }
    default:
    {
        printf("\nAltitude unit should be K F or M:%c", argv[3][0]);
    }
    }
    
    latitude = atof(argv[4]);
    longitude = atof(argv[5]);


    GeoMagnetic = new CGeoMagnetic(sdate, projection);


    CoordGeodetic = new CCoordGeodetic(GeoMagnetic->Geoid, latitude, longitude, alt);

    GeoMagnetic->CoordSpherical->GeodeticToSperical(GeoMagnetic->Ellip, CoordGeodetic);

    GeoMagnetic->GeoMagneticElements->CalculateFieldElements(GeoMagnetic->Ellip, GeoMagnetic->CoordSpherical, CoordGeodetic, GeoMagnetic->TimedMagneticModel);

    delete CoordGeodetic;


    delete GeoMagnetic;
    
    return 0;
}