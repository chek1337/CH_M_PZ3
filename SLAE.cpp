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
			for (int j = ia[i]; j < ia[i + 1]; j++) //��� ��� 100 ���� ���������.
			{
				for (int p = lastj; p < ja[j]; p++) //��� ��� 100 ���� ���������.
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

void SLAE::MatrixVectorMultiplication(double* vectorMult, double* vectorOut)
{
	VectorZeroing(vectorOut);
	for (int i = 0; i < n; i++)
	{
		vectorOut[i] = di[i] * vectorMult[i];
		for (int k = ia[i]; k < ia[i + 1]; k++)
		{
			int j = ja[k];
			vectorOut[i] += al[k] * vectorMult[j];
			vectorOut[j] += au[k] * vectorMult[i];
		}
	}
}

void SLAE::TransposedMatrixVectorMultiplication(double* vectorMult, double* vectorOut)
{
	VectorZeroing(vectorOut);
	for (int i = 0; i < n; i++)
	{
		vectorOut[i] = di[i] * vectorMult[i];
		for (int k = ia[i]; k < ia[i + 1]; k++)
		{
			int j = ja[k];
			vectorOut[i] += au[k] * vectorMult[j];
			vectorOut[j] += al[k] * vectorMult[i];
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

void SLAE::MethodOfConjugateGradientsForSymMatrix()
{
	MatrixVectorMultiplication(x0, tmp1);
	VectorSubtract(b, tmp1, tmp1);
	VectorCopy(tmp1, r);
	VectorCopy(r, z);

	double normB = VectorNorm(b);
	double ak = 0, bk = 0;
	double RelDiscrepancy = CalculateRelativeDiscrepancy(normB);

	for (int curIt = 0; curIt < maxiter and RelDiscrepancy > eps; curIt++)
	{
		printf("Iteration: %d ", curIt + 1);

		double r_rPrev = VectorScalarProduction(r); // ( r(k-1) , r(k-1) )
		MatrixVectorMultiplication(z, tmp1); // A*z(k-1)
		double Az_zPrev = VectorScalarProduction(tmp1, z); //( A*z(k-1) , z(k-1) )
		ak = r_rPrev / Az_zPrev;

		for (int i = 0; i < n; i++) // xk = x(k-1) + ak*z(k-1)
			x[i] = x[i] + ak * z[i];

		for (int i = 0; i < n; i++) // rk = r(k-1) - ak*A*z(k-1)
			r[i] = r[i] - ak * tmp1[i];

		double r_rCur = VectorScalarProduction(r); // ( rk , rk )
		bk = r_rCur / r_rPrev;

		for (int i = 0; i < n; i++) // zk = M^(-1)*rk + bk* z(k-1)
			z[i] = r[i] + bk * z[i];
		RelDiscrepancy = CalculateRelativeDiscrepancy(normB);
		printf("RelDiscrepancy: %.15lf\n", RelDiscrepancy);
	}
}
/**********************************************************/

void SLAE::MethodOfConjugateGradientsForNonSymMatrix()
{
	MatrixVectorMultiplication(x0, tmp1);
	VectorSubtract(b, tmp1, tmp1);
	SolveForward(al, tmp1, tmp2); //+
	SolveBackward(al, tmp2, tmp1); //+
	TransposedMatrixVectorMultiplication(tmp1, tmp2); //+
	SolveForward(au, tmp2, tmp1); //+
	VectorCopy(tmp1, r);
	VectorCopy(r, z);
	MatrixUVectorMultiplication(x0, x);
	double r_rPrev, r_rCur, Newz_zPrev, ak, bk;
	double normB = VectorNorm(b);
	double RelDiscrepancy = CalculateRelativeDiscrepancy(normB);

	for (int curIt = 0; curIt < maxiter and RelDiscrepancy > eps; curIt++)
	{
		printf("Iteration: %d ", curIt + 1);

		r_rPrev = VectorScalarProduction(r); // ( r(k-1) , r(k-1) )
		VectorCopy(z, tmp1); //���� ���, ��� ��� CalculateZ �� ����� ������ �������
		CalculateZ(tmp2); //U^(-T)*A^T*L^(-T)*L^(-1)*A*U^(-1)*z(k-1)       //������-�� �������� z, ����� ����� ������
		Newz_zPrev = VectorScalarProduction(tmp2, z); //( A*z(k-1) , z(k-1) )
		ak = r_rPrev / Newz_zPrev;

		for (int i = 0; i < n; i++) // xk = x(k-1) + ak*z(k-1)
			x[i] = x[i] + ak * z[i];

		for (int i = 0; i < n; i++) // rk = r(k-1) - ak*U^(-T)*A^T*L^(-T)*L^(-1)*A*U^(-1)*z(k-1)
			r[i] = r[i] - ak * tmp2[i];

		r_rCur = VectorScalarProduction(r); // ( rk , rk )
		bk = r_rCur / r_rPrev;

		for (int i = 0; i < n; i++) // zk = M^(-1)*rk + bk* z(k-1)
			z[i] = r[i] + bk * z[i];
		RelDiscrepancy = CalculateRelativeDiscrepancy(normB);
		printf("RelDiscrepancy: %.15lf\n", RelDiscrepancy);
	}
	SolveBackward(au, x, tmp1);
	VectorCopy(tmp1, x);
}


void SLAE::SolveForward(double *lowerTringMat, double *rightVector, double *vectorX) {
	for (int i0, i1, i = 0; i < n; i++)
	{
		double sum = 0;
		i0 = ia[i];
		i1 = ia[i + 1];
		//int j = i - (i1 - i0);
		for ( int k = i0; k < i1; k++)
		{
			int j = ja[k];
			sum += lowerTringMat[k] * vectorX[j];
		}
		vectorX[i] = (rightVector[i] - sum) / di[i];
	}
}

void SLAE::SolveBackward(double* upperTringMat, double* rightVector, double* vectorX) { //�� ����� ������ ���������� ��-�� -=
	for (int i0, i1, i = n - 1; i >= 0; i--)
	{
		i0 = ia[i];
		i1 = ia[i + 1];
		vectorX[i] = rightVector[i]/di[i];
		for (int j, k = i0; k < i1; k++)
		{
			j = ja[k];
			rightVector[j] -= upperTringMat[k] * vectorX[i];
		}
	}
}

void SLAE::MatrixUVectorMultiplication(double *vectorMult, double *vectorOut)
{
	for (int i = 0; i < n; i++)
	{
		vectorOut[i] = vectorMult[i] * di[i];
		for (int j, k = ia[i]; k < ia[i+1]; k++)
		{
			j = ja[k];
			vectorOut[j] += au[k] * vectorMult[i];
		}
	}
}

void SLAE::CalculateZ(double* vectorOut)
{
	SolveBackward(au, tmp1, tmp1);
	MatrixVectorMultiplication(tmp1, tmp2);
	SolveForward(al, tmp2, tmp1);
	SolveBackward(al, tmp1, tmp2);
	TransposedMatrixVectorMultiplication(tmp2, tmp1);
	SolveForward(au, tmp1, vectorOut);
}

double SLAE::CalculateRelativeDiscrepancy(double norm)
{
	return VectorNorm(r) / norm;
}

//**********************************************************************
void SLAE::MethodOfConjugateGradientsWithDiagP()
{
	MatrixVectorMultiplication(x0, tmp1);
	VectorSubtract(b, tmp1, tmp1);
	VectorCopy(tmp1, r);
	VectorConditionalityDiagP(r, z);

	double normB = VectorNorm(b);
	double RelDiscrepancy = CalculateRelativeDiscrepancy(normB);
	double ak=0, bk=0;
	for (int curIt = 0; curIt < maxiter and RelDiscrepancy > eps; curIt++)
	{
		printf("Iteration: %d ", curIt + 1);
		VectorConditionalityDiagP(r, tmp1); //M^(-1)*r(k-1)
		double Mr_rPrev = VectorScalarProduction(tmp1, r); //( M^(-1)*r(k-1) , r(k-1) )
		MatrixVectorMultiplication(z, tmp1); //A*z(k-1)
		double Az_zPrev = VectorScalarProduction(tmp1, z); // ( A*z(k-1) , z(k-1) )
		ak = Mr_rPrev / Az_zPrev;

		for (int i = 0; i < n; i++) // xk = x(k-1) + ak*z(k-1)
			x[i] = x[i] + ak * z[i];
		
		for (int i = 0; i < n; i++) // rk = r(k-1) - ak*A*z(k-1)
			r[i] = r[i] - ak * tmp1[i];

		VectorConditionalityDiagP(r, tmp1); // M^(-1)*rk
		double Mr_rCur = VectorScalarProduction(tmp1, r); //( M^(-1)*rk , rk )
		bk = Mr_rCur / Mr_rPrev;

		for (int i = 0; i < n; i++) // zk = M^(-1)*rk + bk* z(k-1)
			z[i] = tmp1[i] + bk * z[i];
		
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

void SLAE::VectorSubtract(double* first, double* second, double* result)
{
	for (int i = 0; i < n; i++)
	{
		result[i] = first[i] - second[i];
	}
}

void SLAE::VectorCopy(double* from, double* to)
{
	for (int i = 0; i < n; i++)
	{
		to[i] = from[i];
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
	tmp1 = new double[n]();
	tmp2 = new double[n]();
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
	delete[] tmp1;
	delete[] tmp2;
}