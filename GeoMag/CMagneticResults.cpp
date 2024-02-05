#include "CMagneticResults.h"

CMagneticResults::CMagneticResults()
{
	memset(this, 0x00, sizeof(CMagneticResults));
}

CMagneticResults::~CMagneticResults()
{

}

void CMagneticResults::Summation(CLegendreFunction* LegendreFunction, CMagneticModel* MagneticModel, CSphericalHarmonicVariables* SphVariables, CCoordSpherical* CoordSpherical)
{
    int m, n, index;
    double cos_phi;
    
    Bz = 0.0;
    By = 0.0;
    Bx = 0.0;
    for (n = 1; n <= MagneticModel->nMax; n++)
    {
        for (m = 0; m <= n; m++)
        {
            index = (n * (n + 1) / 2 + m);

            /*		    nMax  	(n+2) 	  n     m            m           m
                    Bz =   -SUM (a/r)   (n+1) SUM  [g cos(m p) + h sin(m p)] P (sin(phi))
                                    n=1      	      m=0   n            n           n  */
                                    /* Equation 12 in the WMM Technical report.  Derivative with respect to radius.*/
            Bz -= SphVariables->RelativeRadiusPower[n] *
                (MagneticModel->Main_Field_Coeff_G[index] * SphVariables->cos_mlambda[m] +
                    MagneticModel->Main_Field_Coeff_H[index] * SphVariables->sin_mlambda[m])
                * (double)(n + 1) * LegendreFunction->Pcup[index];

            /*		  1 nMax  (n+2)    n     m            m           m
                    By =    SUM (a/r) (m)  SUM  [g cos(m p) + h sin(m p)] dP (sin(phi))
                               n=1             m=0   n            n           n  */
                               /* Equation 11 in the WMM Technical report. Derivative with respect to longitude, divided by radius. */
            By += SphVariables->RelativeRadiusPower[n] *
                (MagneticModel->Main_Field_Coeff_G[index] * SphVariables->sin_mlambda[m] -
                    MagneticModel->Main_Field_Coeff_H[index] * SphVariables->cos_mlambda[m])
                * (double)(m)*LegendreFunction->Pcup[index];
            /*		   nMax  (n+2) n     m            m           m
                    Bx = - SUM (a/r)   SUM  [g cos(m p) + h sin(m p)] dP (sin(phi))
                               n=1         m=0   n            n           n  */
                               /* Equation 10  in the WMM Technical report. Derivative with respect to latitude, divided by radius. */

            Bx -= SphVariables->RelativeRadiusPower[n] *
                (MagneticModel->Main_Field_Coeff_G[index] * SphVariables->cos_mlambda[m] +
                    MagneticModel->Main_Field_Coeff_H[index] * SphVariables->sin_mlambda[m])
                * LegendreFunction->dPcup[index];



        }
    }

    cos_phi = cos(DEG2RAD(CoordSpherical->phig));
    if (fabs(cos_phi) > 1.0e-10)
    {
        By = By / cos_phi;
    }
    else
        /* Special calculation for component - By - at Geographic poles.
         * If the user wants to avoid using this function,  please make sure that
         * the latitude is not exactly +/-90. An option is to make use the function
         * MAG_CheckGeographicPoles.
         */
    {
        CMagneticResults::SummationSpecial(MagneticModel, SphVariables, CoordSpherical);
    }
}

void CMagneticResults::SummationSpecial(CMagneticModel* MagneticModel, CSphericalHarmonicVariables* SphVariables, CCoordSpherical* CoordSpherical)
{
    int n, index;
    double k, sin_phi, * PcupS, schmidtQuasiNorm1, schmidtQuasiNorm2, schmidtQuasiNorm3;

    PcupS = (double*)malloc((MagneticModel->nMax + 1) * sizeof(double));

    PcupS[0] = 1;
    schmidtQuasiNorm1 = 1.0;

    By = 0.0;
    sin_phi = sin(DEG2RAD(CoordSpherical->phig));

    for (n = 1; n <= MagneticModel->nMax; n++)
    {

        /*Compute the ration between the Gauss-normalized associated Legendre
  functions and the Schmidt quasi-normalized version. This is equivalent to
  sqrt((m==0?1:2)*(n-m)!/(n+m!))*(2n-1)!!/(n-m)!  */

        index = (n * (n + 1) / 2 + 1);
        schmidtQuasiNorm2 = schmidtQuasiNorm1 * (double)(2 * n - 1) / (double)n;
        schmidtQuasiNorm3 = schmidtQuasiNorm2 * sqrt((double)(n * 2) / (double)(n + 1));
        schmidtQuasiNorm1 = schmidtQuasiNorm2;
        if (n == 1)
        {
            PcupS[n] = PcupS[n - 1];
        }
        else
        {
            k = (double)(((n - 1) * (n - 1)) - 1) / (double)((2 * n - 1) * (2 * n - 3));
            PcupS[n] = sin_phi * PcupS[n - 1] - k * PcupS[n - 2];
        }

        /*		  1 nMax  (n+2)    n     m            m           m
                By =    SUM (a/r) (m)  SUM  [g cos(m p) + h sin(m p)] dP (sin(phi))
                           n=1             m=0   n            n           n  */
                           /* Equation 11 in the WMM Technical report. Derivative with respect to longitude, divided by radius. */

        By += SphVariables->RelativeRadiusPower[n] *
            (MagneticModel->Main_Field_Coeff_G[index] * SphVariables->sin_mlambda[1] -
                MagneticModel->Main_Field_Coeff_H[index] * SphVariables->cos_mlambda[1])
            * PcupS[n] * schmidtQuasiNorm3;
    }

    if (PcupS)
    {
        free(PcupS);
    }

}

void CMagneticResults::SecVarSummation(CLegendreFunction* LegendreFunction, CMagneticModel* MagneticModel, CSphericalHarmonicVariables* SphVariables, CCoordSpherical* CoordSpherical)
{
    /*This Function sums the secular variation coefficients to get the secular variation of the Magnetic vector.
    INPUT :  LegendreFunction
                    MagneticModel
                    SphVariables
                    CoordSpherical
    OUTPUT : MagneticResults

    CALLS : MAG_SecVarSummationSpecial

     */
    int m, n, index;
    double cos_phi;
    
    //MagneticModel->SecularVariationUsed = TRUE;
    MagneticModel->SecularVariationUsed = 1;
    
    Bz = 0.0;
    By = 0.0;
    Bx = 0.0;
    
    for (n = 1; n <= MagneticModel->nMaxSecVar; n++)
    {
        for (m = 0; m <= n; m++)
        {
            index = (n * (n + 1) / 2 + m);

            /*		    nMax  	(n+2) 	  n     m            m           m
                    Bz =   -SUM (a/r)   (n+1) SUM  [g cos(m p) + h sin(m p)] P (sin(phi))
                                    n=1      	      m=0   n            n           n  */
                                    /*  Derivative with respect to radius.*/
            Bz -= SphVariables->RelativeRadiusPower[n] *
                (MagneticModel->Secular_Var_Coeff_G[index] * SphVariables->cos_mlambda[m] +
                    MagneticModel->Secular_Var_Coeff_H[index] * SphVariables->sin_mlambda[m])
                * (double)(n + 1) * LegendreFunction->Pcup[index];

            /*		  1 nMax  (n+2)    n     m            m           m
                    By =    SUM (a/r) (m)  SUM  [g cos(m p) + h sin(m p)] dP (sin(phi))
                               n=1             m=0   n            n           n  */
                               /* Derivative with respect to longitude, divided by radius. */
            By += SphVariables->RelativeRadiusPower[n] *
                (MagneticModel->Secular_Var_Coeff_G[index] * SphVariables->sin_mlambda[m] -
                    MagneticModel->Secular_Var_Coeff_H[index] * SphVariables->cos_mlambda[m])
                * (double)(m)*LegendreFunction->Pcup[index];
            /*		   nMax  (n+2) n     m            m           m
                    Bx = - SUM (a/r)   SUM  [g cos(m p) + h sin(m p)] dP (sin(phi))
                               n=1         m=0   n            n           n  */
                               /* Derivative with respect to latitude, divided by radius. */

            Bx -= SphVariables->RelativeRadiusPower[n] *
                (MagneticModel->Secular_Var_Coeff_G[index] * SphVariables->cos_mlambda[m] +
                    MagneticModel->Secular_Var_Coeff_H[index] * SphVariables->sin_mlambda[m])
                * LegendreFunction->dPcup[index];
        }
    }
    cos_phi = cos(DEG2RAD(CoordSpherical->phig));
    if (fabs(cos_phi) > 1.0e-10)
    {
        By = By / cos_phi;
    }
    else
        /* Special calculation for component By at Geographic poles */
    {
        CMagneticResults::SecVarSummationSpecial(MagneticModel, SphVariables, CoordSpherical);
    }

} /*MAG_SecVarSummation*/

void CMagneticResults::SecVarSummationSpecial(CMagneticModel* MagneticModel, CSphericalHarmonicVariables* SphVariables, CCoordSpherical* CoordSpherical)
{
    /*Special calculation for the secular variation summation at the poles.


    INPUT: MagneticModel
               SphVariables
               CoordSpherical
    OUTPUT: MagneticResults
    CALLS : none


     */
    int n, index;
    double k, sin_phi, * PcupS, schmidtQuasiNorm1, schmidtQuasiNorm2, schmidtQuasiNorm3;

    PcupS = (double*)malloc((MagneticModel->nMaxSecVar + 1) * sizeof(double));

    PcupS[0] = 1;
    schmidtQuasiNorm1 = 1.0;

    By = 0.0;
    sin_phi = sin(DEG2RAD(CoordSpherical->phig));

    for (n = 1; n <= MagneticModel->nMaxSecVar; n++)
    {
        index = (n * (n + 1) / 2 + 1);
        schmidtQuasiNorm2 = schmidtQuasiNorm1 * (double)(2 * n - 1) / (double)n;
        schmidtQuasiNorm3 = schmidtQuasiNorm2 * sqrt((double)(n * 2) / (double)(n + 1));
        schmidtQuasiNorm1 = schmidtQuasiNorm2;
        if (n == 1)
        {
            PcupS[n] = PcupS[n - 1];
        }
        else
        {
            k = (double)(((n - 1) * (n - 1)) - 1) / (double)((2 * n - 1) * (2 * n - 3));
            PcupS[n] = sin_phi * PcupS[n - 1] - k * PcupS[n - 2];
        }

        /*		  1 nMax  (n+2)    n     m            m           m
                By =    SUM (a/r) (m)  SUM  [g cos(m p) + h sin(m p)] dP (sin(phi))
                           n=1             m=0   n            n           n  */
                           /* Derivative with respect to longitude, divided by radius. */

        By += SphVariables->RelativeRadiusPower[n] *
            (MagneticModel->Secular_Var_Coeff_G[index] * SphVariables->sin_mlambda[1] -
                MagneticModel->Secular_Var_Coeff_H[index] * SphVariables->cos_mlambda[1])
            * PcupS[n] * schmidtQuasiNorm3;
    }

    if (PcupS)
    {
        free(PcupS);
    }

}/*SecVarSummationSpecial*/

void CMagneticResults::RotateMagneticVector(CCoordSpherical* CoordSpherical, CCoordGeodetic* CoordGeodetic, CMagneticResults* MagneticResults)
{
    double Psi;
    /* Difference between the spherical and Geodetic latitudes */
    Psi = (M_PI / 180) * (CoordSpherical->phig - CoordGeodetic->phi);

    /* Rotate spherical field components to the Geodetic system */
    Bz = MagneticResults->Bx * sin(Psi) + MagneticResults->Bz * cos(Psi);
    Bx = MagneticResults->Bx * cos(Psi) - MagneticResults->Bz * sin(Psi);
    By = MagneticResults->By;

} /*MAG_RotateMagneticVector*/

void CMagneticResults::SphericalToCartesian(CCoordSpherical* CoordSpherical, double* x, double* y, double* z)
{
    double radphi;
    double radlambda;

    radphi = CoordSpherical->phig * (M_PI / 180);
    radlambda = CoordSpherical->lambda * (M_PI / 180);

    *x = CoordSpherical->r * cos(radphi) * cos(radlambda);
    *y = CoordSpherical->r * cos(radphi) * sin(radlambda);
    *z = CoordSpherical->r * sin(radphi);
}