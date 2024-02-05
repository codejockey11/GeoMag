#ifndef _CMagneticModel
#define _CMagneticModel

#include "standard.h"

#include "CDate.h"

#define MAXLINELENGTH 1024

#define CALCULATE_NUMTERMS(N)    (N * ( N + 1 ) / 2 + N)

class CMagneticModel
{
public:

    double EditionDate;
    double epoch; /*Base time of Geomagnetic model epoch (yrs)*/
    char ModelName[32];

    double* Main_Field_Coeff_G; /* C - Gauss coefficients of main geomagnetic model (nT) Index is (n * (n + 1) / 2 + m) */
    double* Main_Field_Coeff_H; /* C - Gauss coefficients of main geomagnetic model (nT) */
    double* Secular_Var_Coeff_G; /* CD - Gauss coefficients of secular geomagnetic model (nT/yr) */
    double* Secular_Var_Coeff_H; /* CD - Gauss coefficients of secular geomagnetic model (nT/yr) */

    int nMax; /* Maximum degree of spherical harmonic model */
    int nMaxSecVar; /* Maximum degree of spherical harmonic secular model */
    int SecularVariationUsed; /* Whether or not the magnetic secular variation vector will be needed by program*/

    double CoefficientFileEndDate;

    int num_terms;

	CMagneticModel();
	
	~CMagneticModel();

    CMagneticModel(int n);

	CMagneticModel(const char* f);

    void TimelyModifyMagneticModel(CDate* UserDate, CMagneticModel* MagneticModel);

private:

	FILE* infile;

	errno_t err;

    char line[MAXLINELENGTH];

    char* eof;

    int a;
    int n;
    int m;
    int index;

    double gnm, hnm, dgnm, dhnm;
};

#endif
