#include "CGeoMagnetic.h"

CGeoMagnetic::CGeoMagnetic()
{
	memset(this, 0x00, sizeof(CGeoMagnetic));
}

CGeoMagnetic::~CGeoMagnetic()
{
	delete UserDate;

	delete Geoid;

	delete Ellip;

	delete CoordSpherical;

	delete GeoMagneticElements;

	delete MagneticModel;

	delete TimedMagneticModel;
}

CGeoMagnetic::CGeoMagnetic(double sdate, const char projection)
{
	memset(this, 0x00, sizeof(CGeoMagnetic));

	UserDate = new CDate(sdate);

	Geoid = new CGeoid(projection);

	Ellip = new CEllipsoid();

	CoordSpherical = new CCoordSpherical();

	GeoMagneticElements = new CGeoMagneticElements();

	MagneticModel = new CMagneticModel("WMM.COF");

	TimedMagneticModel = new CMagneticModel(MagneticModel->num_terms);

	TimedMagneticModel->TimelyModifyMagneticModel(UserDate, MagneticModel);
}