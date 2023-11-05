#pragma once
#pragma once
#include <cstdio>
#include <math.h>

#define REALOUTD "%.0f\t"

class SLAE {
public:
	int n, maxiter;
	double eps, ak, bk;
	double* al, * au, * di;
	double *x,*x0, * b, *r, *z, *tmp;
	int* ia, *ja, nProfile;

	void Input(FILE* paramf, FILE* iaf, FILE* jaf, FILE* alf, FILE* auf, FILE* dif, FILE* bf);

	void MatrixVectorMultiplication(double* vectorMult, double* vectorOut);
	void TransposedMatrixVectorMultiplication(double* vectorMult, double* vectorOut);
	void MethodOfConjugateGradients();
	void InitialApproximation();
	double CalculateRelativeDiscrepancy();

	
	void Calculate_ak();
	void Calculate_xk();
	void Calculate_rk();
	void Calculate_zk();

	void VecotorSubtract(double* first, double* second);
	double VectorScalarProduction(double* vector);
	double VectorScalarProduction(double* vector1, double* vector2);
	double VectorNorm(double* vector);
	void VectorCopy(double* first, double* second);
	void VectorZeroing(double* vector);

	void OutputDense();

	//void GenerateHilbertMatrix();

	void VectorOutputSolution(FILE* out);

protected:
	void AllocateMemory();
	void ClearMemory();
};