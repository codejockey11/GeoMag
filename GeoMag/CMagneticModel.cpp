#include "CMagneticModel.h"

CMagneticModel::CMagneticModel()
{
	memset(this, 0x00, sizeof(CMagneticModel));

}

CMagneticModel::CMagneticModel(int n)
{
	memset(this, 0x00, sizeof(CMagneticModel));

	num_terms = n;

	Main_Field_Coeff_G = (double*)calloc((num_terms + 1), sizeof(double));
	Main_Field_Coeff_H = (double*)calloc((num_terms + 1), sizeof(double));
	
	Secular_Var_Coeff_G = (double*)calloc((num_terms + 1), sizeof(double));
	Secular_Var_Coeff_H = (double*)calloc((num_terms + 1), sizeof(double));
}

CMagneticModel::CMagneticModel(const char* f)
{
	memset(this, 0x00, sizeof(CMagneticModel));

	if (f == nullptr)
	{
		printf("\nMagnetic Model file name not provided");

		return;
	}

	err = fopen_s(&infile, f, "rt");

	if (err)
	{
		printf("\nMagnetic Model file not found:%s", f);

		return;
	}

	// skipping the header while checking for end-of-file
	if (fgets(line, MAXLINELENGTH, infile) == eof)
	{
		printf("\nMagnetic Model file is empty");

		return;
	}

	eof = fgets(line, MAXLINELENGTH, infile);

	while (eof != nullptr)
	{
		// checking first value
		a = sscanf_s(line, "%d", &n);

		if (n > 12)
		{
			break;
		}

		if (n > nMax)
		{
			nMax = n;
		}

		eof = fgets(line, MAXLINELENGTH, infile);
	}

	fclose(infile);

	nMaxSecVar = nMax;

	num_terms = CALCULATE_NUMTERMS(nMax);

	Main_Field_Coeff_G = (double*)calloc((num_terms + 1), sizeof(double));
	Main_Field_Coeff_H = (double*)calloc((num_terms + 1), sizeof(double));
	
	Secular_Var_Coeff_G = (double*)calloc((num_terms + 1), sizeof(double));
	Secular_Var_Coeff_H = (double*)calloc((num_terms + 1), sizeof(double));

	CoefficientFileEndDate = epoch + 5;



	err = fopen_s(&infile, f, "rt");

	Main_Field_Coeff_H[0] = 0.0;
	Main_Field_Coeff_G[0] = 0.0;
	Secular_Var_Coeff_H[0] = 0.0;
	Secular_Var_Coeff_G[0] = 0.0;

	eof = fgets(line, MAXLINELENGTH, infile);
	
	sscanf(line, "%lf%s", &epoch, ModelName);
	//MagneticModel->epoch = epoch;

	eof = fgets(line, MAXLINELENGTH, infile);

	while (eof != nullptr)
	{
		if (line[0] == '9')
		{
			break;
		}

		sscanf(line, "%d%d%lf%lf%lf%lf", &n, &m, &gnm, &hnm, &dgnm, &dhnm);
		
		if (m <= n)
		{
			index = (n * (n + 1) / 2 + m);
			
			Main_Field_Coeff_G[index] = gnm;
			Secular_Var_Coeff_G[index] = dgnm;
			Main_Field_Coeff_H[index] = hnm;
			Secular_Var_Coeff_H[index] = dhnm;
		}

		eof = fgets(line, MAXLINELENGTH, infile);
	}

	CoefficientFileEndDate = epoch + 5;

	fclose(infile);
}

CMagneticModel::~CMagneticModel()
{
	if (Main_Field_Coeff_G)
	{
		free(Main_Field_Coeff_G);
	}
	
	if (Main_Field_Coeff_H)
	{
		free(Main_Field_Coeff_H);
	}
	
	if (Secular_Var_Coeff_G)
	{
		free(Secular_Var_Coeff_G);
	}
	
	if (Secular_Var_Coeff_H)
	{
		free(Secular_Var_Coeff_H);
	}
}

void CMagneticModel::TimelyModifyMagneticModel(CDate* UserDate, CMagneticModel* MagneticModel)
{
	int n, m, index, a, b;
	
	EditionDate = MagneticModel->EditionDate;
	
	epoch = MagneticModel->epoch;
	
	nMax = MagneticModel->nMax;
	
	nMaxSecVar = MagneticModel->nMaxSecVar;
	
	a = nMaxSecVar;
	
	b = (a * (a + 1) / 2 + a);
	
	strcpy_s(ModelName, 32, MagneticModel->ModelName);
	
	for (n = 1; n <= MagneticModel->nMax; n++)
	{
		for (m = 0; m <= n; m++)
		{
			index = (n * (n + 1) / 2 + m);
			
			if (index <= b)
			{
				Main_Field_Coeff_H[index] = MagneticModel->Main_Field_Coeff_H[index] + (UserDate->DecimalYear - MagneticModel->epoch) * MagneticModel->Secular_Var_Coeff_H[index];
				Main_Field_Coeff_G[index] = MagneticModel->Main_Field_Coeff_G[index] + (UserDate->DecimalYear - MagneticModel->epoch) * MagneticModel->Secular_Var_Coeff_G[index];
				
				/* We need a copy of the secular var coef to calculate secular change */
				Secular_Var_Coeff_H[index] = MagneticModel->Secular_Var_Coeff_H[index];
				Secular_Var_Coeff_G[index] = MagneticModel->Secular_Var_Coeff_G[index];
			}
			else
			{
				Main_Field_Coeff_H[index] = MagneticModel->Main_Field_Coeff_H[index];
				Main_Field_Coeff_G[index] = MagneticModel->Main_Field_Coeff_G[index];
			}
		}
	}
}