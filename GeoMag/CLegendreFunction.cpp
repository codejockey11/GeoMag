#include "CLegendreFunction.h"

CLegendreFunction::CLegendreFunction()
{
	memset(this, 0x00, sizeof(CLegendreFunction));
}

CLegendreFunction::CLegendreFunction(int n)
{
	memset(this, 0x00, sizeof(CLegendreFunction));

	NumTerms = n;

    Pcup = (double*)calloc((NumTerms + 1), sizeof(double));
    
	dPcup = (double*)calloc((NumTerms + 1), sizeof(double));
}

CLegendreFunction::~CLegendreFunction()
{
	if (Pcup)
	{
		free(Pcup);
	}

	if (dPcup)
	{
		free(dPcup);
	}
}

void CLegendreFunction::AssociatedLegendreFunction(CCoordSpherical* CoordSpherical, int nMax)
{
	double sin_phi = sin(DEG2RAD(CoordSpherical->phig)); /* sin  (geocentric latitude) */

	/* If nMax is less tha 16 or at the poles */
	if (nMax <= 16 || (1 - fabs(sin_phi)) < 1.0e-10)
	{
        CLegendreFunction::PcupLow(sin_phi, nMax);
	}
	else
	{
        CLegendreFunction::PcupHigh(sin_phi, nMax);
	}
}

void CLegendreFunction::PcupLow(double x, int nMax)
{
    int n, m, index, index1, index2, NumTerms;
    
    double k, z;
    
    double* schmidtQuasiNorm;
    
    Pcup[0] = 1.0;
    dPcup[0] = 0.0;
    
    /*sin (geocentric latitude) - sin_phi */
    z = sqrt((1.0 - x) * (1.0 + x));

    NumTerms = ((nMax + 1) * (nMax + 2) / 2);
    schmidtQuasiNorm = (double*)malloc((NumTerms + 1) * sizeof(double));

    /*	 First,	Compute the Gauss-normalized associated Legendre  functions*/
    for (n = 1; n <= nMax; n++)
    {
        for (m = 0; m <= n; m++)
        {
            index = (n * (n + 1) / 2 + m);
            if (n == m)
            {
                index1 = (n - 1) * n / 2 + m - 1;
                Pcup[index] = z * Pcup[index1];
                dPcup[index] = z * dPcup[index1] + x * Pcup[index1];
            }
            else if (n == 1 && m == 0)
            {
                index1 = (n - 1) * n / 2 + m;
                Pcup[index] = x * Pcup[index1];
                dPcup[index] = x * dPcup[index1] - z * Pcup[index1];
            }
            else if (n > 1 && n != m)
            {
                index1 = (n - 2) * (n - 1) / 2 + m;
                index2 = (n - 1) * n / 2 + m;
                if (m > n - 2)
                {
                    Pcup[index] = x * Pcup[index2];
                    dPcup[index] = x * dPcup[index2] - z * Pcup[index2];
                }
                else
                {
                    k = (double)(((n - 1) * (n - 1)) - (m * m)) / (double)((2 * n - 1) * (2 * n - 3));
                    Pcup[index] = x * Pcup[index2] - k * Pcup[index1];
                    dPcup[index] = x * dPcup[index2] - z * Pcup[index2] - k * dPcup[index1];
                }
            }
        }
    }
    /* Compute the ration between the the Schmidt quasi-normalized associated Legendre
     * functions and the Gauss-normalized version. */

    schmidtQuasiNorm[0] = 1.0;
    for (n = 1; n <= nMax; n++)
    {
        index = (n * (n + 1) / 2);
        index1 = (n - 1) * n / 2;
        /* for m = 0 */
        schmidtQuasiNorm[index] = schmidtQuasiNorm[index1] * (double)(2 * n - 1) / (double)n;

        for (m = 1; m <= n; m++)
        {
            index = (n * (n + 1) / 2 + m);
            index1 = (n * (n + 1) / 2 + m - 1);
            schmidtQuasiNorm[index] = schmidtQuasiNorm[index1] * sqrt((double)((n - m + 1) * (m == 1 ? 2 : 1)) / (double)(n + m));
        }
    }

    /* Converts the  Gauss-normalized associated Legendre
              functions to the Schmidt quasi-normalized version using pre-computed
              relation stored in the variable schmidtQuasiNorm */

    for (n = 1; n <= nMax; n++)
    {
        for (m = 0; m <= n; m++)
        {
            index = (n * (n + 1) / 2 + m);
            Pcup[index] = Pcup[index] * schmidtQuasiNorm[index];
            dPcup[index] = -dPcup[index] * schmidtQuasiNorm[index];
            /* The sign is changed since the new WMM routines use derivative with respect to latitude
            insted of co-latitude */
        }
    }

    if (schmidtQuasiNorm)
    {
        free(schmidtQuasiNorm);
    }
}

void CLegendreFunction::PcupHigh(double x, int nMax)
{
    double pm2, pm1, pmm, plm, rescalem, z, scalef;
    double* f1, * f2, * PreSqr;
    int k, kstart, m, n, NumTerms;

    NumTerms = ((nMax + 1) * (nMax + 2) / 2);


    if (fabs(x) == 1.0)
    {
        printf("Error in PcupHigh: derivative cannot be calculated at poles\n");
        //return FALSE;
        return;
    }


    f1 = (double*)malloc((NumTerms + 1) * sizeof(double));

    PreSqr = (double*)malloc((NumTerms + 1) * sizeof(double));

    f2 = (double*)malloc((NumTerms + 1) * sizeof(double));

    scalef = 1.0e-280;

    for (n = 0; n <= 2 * nMax + 1; ++n)
    {
        PreSqr[n] = sqrt((double)(n));
    }

    k = 2;

    for (n = 2; n <= nMax; n++)
    {
        k = k + 1;
        f1[k] = (double)(2 * n - 1) / (double)(n);
        f2[k] = (double)(n - 1) / (double)(n);
        for (m = 1; m <= n - 2; m++)
        {
            k = k + 1;
            f1[k] = (double)(2 * n - 1) / PreSqr[n + m] / PreSqr[n - m];
            f2[k] = PreSqr[n - m - 1] * PreSqr[n + m - 1] / PreSqr[n + m] / PreSqr[n - m];
        }
        k = k + 2;
    }

    /*z = sin (geocentric latitude) */
    z = sqrt((1.0 - x) * (1.0 + x));
    pm2 = 1.0;
    Pcup[0] = 1.0;
    dPcup[0] = 0.0;
    
    if (nMax == 0)
        //return FALSE;
        return;
    
    pm1 = x;
    Pcup[1] = pm1;
    dPcup[1] = z;
    k = 1;

    for (n = 2; n <= nMax; n++)
    {
        k = k + n;
        plm = f1[k] * x * pm1 - f2[k] * pm2;
        Pcup[k] = plm;
        dPcup[k] = (double)(n) * (pm1 - x * plm) / z;
        pm2 = pm1;
        pm1 = plm;
    }

    pmm = PreSqr[2] * scalef;
    rescalem = 1.0 / scalef;
    kstart = 0;

    for (m = 1; m <= nMax - 1; ++m)
    {
        rescalem = rescalem * z;

        /* Calculate Pcup(m,m)*/
        kstart = kstart + m + 1;
        pmm = pmm * PreSqr[2 * m + 1] / PreSqr[2 * m];
        Pcup[kstart] = pmm * rescalem / PreSqr[2 * m + 1];
        dPcup[kstart] = -((double)(m)*x * Pcup[kstart] / z);
        pm2 = pmm / PreSqr[2 * m + 1];
        /* Calculate Pcup(m+1,m)*/
        k = kstart + m + 1;
        pm1 = x * PreSqr[2 * m + 1] * pm2;
        Pcup[k] = pm1 * rescalem;
        dPcup[k] = ((pm2 * rescalem) * PreSqr[2 * m + 1] - x * (double)(m + 1) * Pcup[k]) / z;
        /* Calculate Pcup(n,m)*/
        for (n = m + 2; n <= nMax; ++n)
        {
            k = k + n;
            plm = x * f1[k] * pm1 - f2[k] * pm2;
            Pcup[k] = plm * rescalem;
            dPcup[k] = (PreSqr[n + m] * PreSqr[n - m] * (pm1 * rescalem) - (double)(n)*x * Pcup[k]) / z;
            pm2 = pm1;
            pm1 = plm;
        }
    }

    /* Calculate Pcup(nMax,nMax)*/
    rescalem = rescalem * z;
    kstart = kstart + m + 1;
    pmm = pmm / PreSqr[2 * nMax];
    
    Pcup[kstart] = pmm * rescalem;
    dPcup[kstart] = -(double)(nMax)*x * Pcup[kstart] / z;
    
    free(f1);
    free(PreSqr);
    free(f2);
}