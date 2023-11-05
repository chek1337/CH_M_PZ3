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

//**********************************************************************
//void SLAE::MethodOfConjugateGradients()
//{
//	MatrixVectorMultiplication(x0, tmp);
//	VecotorSubtract(b, tmp);
//	VectorCopy(tmp, r);
//	VectorCopy(r, z);
//
//	double normB = VectorNorm(b);
//	double ak = 0, bk = 0;
//	double RelDiscrepancy = CalculateRelativeDiscrepancy(normB);
//	
//	for (int curIt = 0; curIt < maxiter and RelDiscrepancy > eps; curIt++)
//	{
//		printf("Iteration: %d ", curIt + 1);
//		
//		double r_rPrev = VectorScalarProduction(r); // ( r(k-1) , r(k-1) )
//		MatrixVectorMultiplication(z, tmp); // A*z(k-1)
//		double Az_zPrev = VectorScalarProduction(tmp, z); //( A*z(k-1) , z(k-1) )
//		ak = r_rPrev / Az_zPrev;
//
//		for (int i = 0; i < n; i++) // xk = x(k-1) + ak*z(k-1)
//			x[i] = x[i] + ak * z[i];
//
//		for (int i = 0; i < n; i++) // rk = r(k-1) - ak*A*z(k-1)
//			r[i] = r[i] - ak * tmp[i];
//
//		double r_rCur = VectorScalarProduction(r); // ( rk , rk )
//		bk = r_rCur / r_rPrev;
//
//		for (int i = 0; i < n; i++) // zk = M^(-1)*rk + bk* z(k-1)
//			z[i] = r[i] + bk * z[i];
//		RelDiscrepancy = CalculateRelativeDiscrepancy(normB);
//		printf("RelDiscrepancy: %.15lf\n", RelDiscrepancy);
//	}
//}

void SLAE::MethodOfConjugateGradients()
{
	MatrixVectorMultiplication(x0, tmp);
	VecotorSubtract(b, tmp);
	VectorCopy(tmp, r);
	VectorCopy(r, z);

	double normB = VectorNorm(b);
	double ak = 0, bk = 0;
	double RelDiscrepancy = CalculateRelativeDiscrepancy(normB);

	for (int curIt = 0; curIt < maxiter and RelDiscrepancy > eps; curIt++)
	{
		printf("Iteration: %d ", curIt + 1);

		double r_rPrev = VectorScalarProduction(r); // ( r(k-1) , r(k-1) )
		MatrixVectorMultiplication(z, tmp); // A*z(k-1)
		double Az_zPrev = VectorScalarProduction(tmp, z); //( A*z(k-1) , z(k-1) )
		ak = r_rPrev / Az_zPrev;

		for (int i = 0; i < n; i++) // xk = x(k-1) + ak*z(k-1)
			x[i] = x[i] + ak * z[i];

		for (int i = 0; i < n; i++) // rk = r(k-1) - ak*A*z(k-1)
			r[i] = r[i] - ak * tmp[i];

		double r_rCur = VectorScalarProduction(r); // ( rk , rk )
		bk = r_rCur / r_rPrev;

		for (int i = 0; i < n; i++) // zk = M^(-1)*rk + bk* z(k-1)
			z[i] = r[i] + bk * z[i];
		RelDiscrepancy = CalculateRelativeDiscrepancy(normB);
		printf("RelDiscrepancy: %.15lf\n", RelDiscrepancy);
	}
}

double SLAE::CalculateRelativeDiscrepancy(double norm)
{
	return VectorNorm(r) / norm;
}

//**********************************************************************
void SLAE::MethodOfConjugateGradientsWithDiagP()
{
	MatrixVectorMultiplication(x0, tmp);
	VecotorSubtract(b, tmp);
	VectorCopy(tmp, r);
	VectorConditionalityDiagP(r, z);

	double normB = VectorNorm(b);
	double RelDiscrepancy = CalculateRelativeDiscrepancy(normB);
	double ak=0, bk=0;
	for (int curIt = 0; curIt < maxiter and RelDiscrepancy > eps; curIt++)
	{
		printf("Iteration: %d ", curIt + 1);
		VectorConditionalityDiagP(r, tmp); //M^(-1)*r(k-1)
		double Mr_rPrev = VectorScalarProduction(tmp, r); //( M^(-1)*r(k-1) , r(k-1) )
		MatrixVectorMultiplication(z, tmp); //A*z(k-1)
		double Az_zPrev = VectorScalarProduction(tmp, z); // ( A*z(k-1) , z(k-1) )
		ak = Mr_rPrev / Az_zPrev;

		for (int i = 0; i < n; i++) // xk = x(k-1) + ak*z(k-1)
			x[i] = x[i] + ak * z[i];
		
		for (int i = 0; i < n; i++) // rk = r(k-1) - ak*A*z(k-1)
			r[i] = r[i] - ak * tmp[i];

		MatrixVectorMultiplication(r, tmp); // M^(-1)*rk
		double Mr_rCur = VectorScalarProduction(tmp, r); //( M^(-1)*rk , rk )
		bk = Mr_rCur / Mr_rPrev;

		for (int i = 0; i < n; i++) // zk = M^(-1)*rk + bk* z(k-1)
			z[i] = tmp[i] + bk * z[i];
		
		RelDiscrepancy = CalculateRelativeDiscrepancy(normB);
		printf("RelDiscrepancy: %.15lf\n", RelDiscrepancy);
	}
}

void SLAE::VectorConditionalityDiagP(double *vectorIn, double *vectorOut)
{
	for (int i = 0; i < n; i++)
	{
		vectorOut[i] = vectorIn[i]/ sqrt(di[i]);
	}
}

//**********************************************************************

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