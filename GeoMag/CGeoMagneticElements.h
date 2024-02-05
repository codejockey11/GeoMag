#ifndef _CGeoMagneticElements
#define _CGeoMagneticElements

#include "standard.h"

#include "CEllipsoid.h"
#include "CCoordSpherical.h"
#include "CCoordGeodetic.h"
#include "CMagneticModel.h"

#include "CLegendreFunction.h"
#include "CSphericalHarmonicVariables.h"
#include "CMagneticResults.h"

class CGeoMagneticElements
{
public:

    double Decl;    /*  1. Angle between the magnetic field vector and true north, positive east */
    double Incl;    /*  2. Angle between the magnetic field vector and the horizontal plane, positive down */
    double F;       /*  3. Magnetic Field Strength */
    double H;       /*  4. Horizontal Magnetic Field Strength */
    double X;       /*  5. Northern component of the magnetic field vector */
    double Y;       /*  6. Eastern component of the magnetic field vector */
    double Z;       /*  7. Downward component of the magnetic field vector */
    double GV;      /*  8. The Grid Variation */
    double Decldot; /*  9. Yearly Rate of change in declination */
    double Incldot; /* 10. Yearly Rate of change in inclination */
    double Fdot;    /* 11. Yearly rate of change in Magnetic field strength */
    double Hdot;    /* 12. Yearly rate of change in horizontal field strength */
    double Xdot;    /* 13. Yearly rate of change in the northern component */
    double Ydot;    /* 14. Yearly rate of change in the eastern component */
    double Zdot;    /* 15. Yearly rate of change in the downward component */
    double GVdot;   /* 16. Yearly rate of change in grid variation */

	CGeoMagneticElements();

	~CGeoMagneticElements();

    void CalculateFieldElements(CEllipsoid* Ellip, CCoordSpherical* CoordSpherical, CCoordGeodetic* CoordGeodetic, CMagneticModel* TimedMagneticModel);

    void CalculateGeoMagneticElements(CMagneticResults* MagneticResultsGeo);

    void CalculateSecularVariationElements(CMagneticResults* MagneticVariation);

private:

    int NumTerms;

    CLegendreFunction* LegendreFunction;
    CSphericalHarmonicVariables* SphVariables;
    
    CMagneticResults* MagneticResultsSph;
    CMagneticResults* MagneticResultsGeo;
    CMagneticResults* MagneticResultsSphVar;
    CMagneticResults* MagneticResultsGeoVar;

};

#endif