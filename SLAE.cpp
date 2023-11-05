#include "slae.h"


void SLAE::Input(FILE *paramf, FILE* iaf, FILE* jaf, FILE* alf, FILE* auf, FILE*dif, FILE *bf) {
	fscanf_s(paramf, "%d", &n);
	fscanf_s(paramf, "%d", &maxiter);
	fscanf_s(paramf, "%lf", &eps);

	ia = new int[n + 1];
	for (int i = 0; i <= n; i++)
		fscanf_s(iaf, "%d", &ia[i]);
	nProfile = ia[n] - ia[0];
	AllocateMemory();
	for (int i = 0; i < nProfile; i++)
		fscanf_s(jaf, "%d", &ja[i]);
	if (ia[0]) {
		for (int i = 0; i <= n; i++)
			ia[i]--;
		for (int i = 0; i < nProfile; i++)
			ja[i]--;
	}

	for (int i = 0; i < nProfile; i++)
		fscanf_s(alf, "%lf", &al[i]);

	for (int i = 0; i < nProfile; i++)
		fscanf_s(auf, "%lf", &au[i]);

	for (int i = 0; i < n; i++)
		fscanf_s(dif, "%lf", &di[i]);

	for (int i = 0; i < n; i++)
		fscanf_s(bf, "%lf", &b[i]);

	//for (int i = 0; i < n; i++)
	//	x0[i] = 0;
}

void SLAE::OutputDense()
{
	int flagfound = 0;
	for (int i = 0; i < n; i++)
	{	
		int k = ia[i + 1] - ia[i];
		if (k == 0)
		{
			for (int j = 0; j < i; j++)
			{
				printf(REALOUTD, 0.0);
			}
		}
		else
		{
			int lastj = 0;
			for (int j = ia[i]; j < ia[i + 1]; j++) //Вот это 100 проц правильно.
			{
				for (int p = lastj; p < ja[j]; p++) //Вот это 100 проц правильно.
				{
					printf(REALOUTD, 0.0);
				}
				printf(REALOUTD, al[j]);
				lastj = ja[j] + 1;
			}
			for (int j = lastj; j < i ; j++) //??
			{
				printf(REALOUTD, 0.0);
			}
		}

		printf(REALOUTD, di[i]);

		for (int j = i + 1; j < n; j++)
		{
			k = ia[j + 1] - ia[j];
			if (k == 0) {
				printf(REALOUTD, 0.0);
			}
			else
			{
				flagfound = 0;
				for (k = ia[j]; k < ia[j + 1]; k++)
				{

					if (ja[k] == i)
					{
						printf(REALOUTD, au[k]);
						flagfound = 1;
						break;				
					}
				}
				if(flagfound == 0)
					printf(REALOUTD, 0.0);
			}
		}
		printf("\n");
	}

}

void SLAE::MatrixVectorMultiplication(double*vectorMult, double* vectorOut)
{
	VectorZeroing(vectorOut);
	for (int i = 0; i < n; i++)
	{
		vectorOut[i] = di[i] * vectorMult[i];
		for (int j = ia[i]; j < ia[i + 1]; j++)
		{
			vectorOut[i] += al[j] * vectorMult[ja[j]];
			vectorOut[ja[j]] += au[j] * vectorMult[i];
		}
	}
}

void SLAE::InitialApproximation()
{
	MatrixVectorMultiplication(x0, tmp);
	VecotorSubtract(b, tmp);
	VectorCopy(tmp, r);
	VectorCopy(r, z);
}

void SLAE::Calculate_ak()
{
	MatrixVectorMultiplication(z, tmp);
	ak = VectorScalarProduction(r) / VectorScalarProduction(tmp, z);
}

double SLAE::CalculateRelativeDiscrepancy()
{
	return VectorNorm(r)/VectorNorm(b);
}

void SLAE::Calculate_xk()
{
	for (int i = 0; i < n; i++)
	{
		x[i] = x[i] + ak * z[i];
	}
}

void SLAE::Calculate_rk()
{
	MatrixVectorMultiplication(z, tmp); //Ну вообще в tmp и так лежит уже A*z
	for (int i = 0; i < n; i++)
	{
		r[i] = r[i] - ak * tmp[i];;
	}
}

void SLAE::Calculate_zk()
{
	for (int i = 0; i < n; i++)
	{
		z[i] = r[i] + bk * z[i];;
	}
}

void SLAE::MethodOfConjugateGradients()
{
	InitialApproximation();
	double RelDiscrepancy = CalculateRelativeDiscrepancy();
	double rk_rkPrev = 0;
	double rk_rkCur = 0;
	for (int curIt = 0; curIt < maxiter and RelDiscrepancy > eps; curIt++)
	{
		printf("Iteration: %d ", curIt + 1);
		Calculate_ak();
		Calculate_xk();
		rk_rkPrev = VectorScalarProduction(r);
		Calculate_rk();
		rk_rkCur = VectorScalarProduction(r);
		bk = rk_rkCur / rk_rkPrev;
		Calculate_zk();
		RelDiscrepancy = CalculateRelativeDiscrepancy();
		printf("RelDiscrepancy: %.15lf\n", RelDiscrepancy);
	}
}

void SLAE::VecotorSubtract(double* first, double* second)
{
	for (int i = 0; i < n; i++)
	{
		tmp[i] = first[i] - second[i];
	}
}

void SLAE::VectorCopy(double* first, double* second)
{
	for (int i = 0; i < n; i++)
	{
		second[i] = first[i];
	}
}

double SLAE::VectorScalarProduction(double* vector1, double* vector2)
{
	double prod = 0;
	for (int i = 0; i < n; i++)
	{
		prod += vector1[i] * vector2[i];
	}
	return prod;
}

double SLAE::VectorScalarProduction(double* vector)
{
	double prod = 0;
	for (int i = 0; i < n; i++)
	{
		prod += vector[i] * vector[i];
	}
	return prod;
}

double SLAE::VectorNorm(double* vector)
{
	double norm = 0;
	for (int i = 0; i < n; i++)
	{
		norm += vector[i] * vector[i];
	}
	return sqrt(norm);
}

void SLAE::VectorZeroing(double *vector)
{
	for (int i = 0; i < n; i++)
	{
		vector[i] = 0;
	}
}



void SLAE::TransposedMatrixVectorMultiplication(double* vectorMult, double* vectorOut)
{
	VectorZeroing(vectorOut);
	for (int i = 0; i < n; i++)
	{
		vectorOut[i] = di[i] * vectorMult[i];
		for (int j = ia[i]; j < ia[i + 1]; j++)
		{
			vectorOut[i] += au[j] * vectorMult[ja[j]];
			vectorOut[ja[j]] += al[j] * vectorMult[i];
		}
	}
}

void SLAE::VectorOutputSolution(FILE *out)
{
	for (int i = 0; i < n; i++)
	{
		fprintf(out, "%lf\n", x[i]);
	}
	printf("\n");
	ClearMemory();
}

void SLAE::AllocateMemory() 
{
	al = new double[nProfile]();
	au = new double[nProfile]();
	di = new double[n]();
	ja = new int[nProfile]();
	b = new double[n]();
	x = new double[n]();
	x0 = new double[n]();
	r = new double[n]();
	z = new double[n]();
	tmp = new double[n]();
}

void SLAE::ClearMemory()
{
	delete [] al;
	delete[] au;
	delete[] di;
	delete[] ja;
	delete[] b;
	delete[] x;
	delete[] x0;
	delete[] z;
	delete[] tmp;
}