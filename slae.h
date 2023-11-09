#pragma once
#pragma once
#include <cstdio>
#include <math.h>

#define REALOUTD "%.0f\t"

class SLAE {
public:
	int n, maxiter;
	double eps;
	double* al, * au, * di;
	double* alLU, * auLU, * diLU;
	double* x, * x0, * b;
	double *r, *z, *tmp1;
	int* ia, *ja, nProfile;

	void Input(FILE* paramf, FILE* iaf, FILE* jaf, FILE* alf, FILE* auf, FILE* dif, FILE* bf);

	void MatrixVectorMultiplication(double* vectorMult, double* vectorOut);
	void TransposedMatrixVectorMultiplication(double* vectorMult, double* vectorOut);
	double CalculateRelativeDiscrepancy(double norm);

	void MethodOfConjugateGradientsForSymMatrix();
	void MethodOfConjugateGradientsForNonSymMatrix();

	void MethodOfConjugateGradientsForSymMatrixWithDiagP();
	void MethodOfConjugateGradientsForNonSymMatrixWithDiagP();
	void VectorConditionalityDiagP(double* vectorIn, double* vectorOut);

	void CalculateLU();

	void SolveForward(double* lowerTringMat, double* rightVector, double* vectorX);
	void SolveBackward(double* upperTringMat, double* rightVector, double* vectorX);
	void MatrixUVectorMultiplication(double* vectorMult, double* vectorOut);
	void CalculateZ(double* vectorOut);

	void VectorSubtract(double* first, double* second, double* result);
	double VectorScalarProduction(double* vector);
	double VectorScalarProduction(double* vector1, double* vector2);
	double VectorNorm(double* vector);
	void VectorCopy(double* first, double* second);
	void VectorZeroing(double* vector);

	void OutputDense();
	void OutputLUDense();

	//void GenerateHilbertMatrix();

	void VectorOutputSolution(FILE* out);

protected:
	void AllocateMemory();
	void ClearMemory();
};