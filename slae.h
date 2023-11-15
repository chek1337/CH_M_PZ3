#pragma once
#pragma once
#include <cstdio>
#include <math.h>

#define REALOUTD "%.7f\t"

class SLAE {
public:
	int n, maxiter = 10000, nProfile = 0;
	double  eps = 1e-13;
	double* al, * au, * di;
	double* alLU, * auLU, * diLU;
	double* x, * x0, * b, * xtrue;
	double* r, * z, * tmp1;
	int* ia, * ja;

	void Input(FILE* paramf, FILE* iaf, FILE* jaf, FILE* alf, FILE* auf, FILE* dif, FILE* bf);

	void MatrixVectorMultiplication(double* vectorMult, double* vectorOut);
	void TransposedMatrixVectorMultiplication(double* vectorMult, double* vectorOut);
	double CalculateRelativeDiscrepancyWithR(double norm);
	double CalculateRelativeDiscrepancy(double norm);

	void MethodOfConjugateGradientsForSymMatrix();
	void MethodOfConjugateGradientsForNonSymMatrixAtA();
	//void MethodOfConjugateGradientsForNonSymMatrix();

	void MethodOfConjugateGradientsForSymMatrixWithDiagP();
	//void MethodOfConjugateGradientsForNonSymMatrixWithDiagP();
	void MethodOfConjugateGradientsForNonSymMatrixAtAWithDiagP();

	void MethodOfConjugateGradientsForSymMatrixWithLuP();
	//void MethodOfConjugateGradientsForNonSymMatrixWithLuP();
	void MethodOfConjugateGradientsForNonSymMatrixAtAWithLuP();

	void VectorConditionalityForSymMatrixDiagP(double* vectorIn, double* vectorOut);
	void VectorConditionalityForNonSymMatrixDiagP(double* vectorIn, double* vectorOut);

	void CalculateLU();

	void GenerateHilbertMatrix(int size);

	void SolveForward(double* lowerTringMat, double* diag, double* rightVector, double* vectorX);
	void SolveBackward(double* upperTringMat, double* diag, double* rightVector, double* vectorX);


	void SolveForwardLU(double* lowerTringMat, double* rightVector, double* vectorX);
	void SolveBackwardLU(double* upperTringMat, double* rightVector, double* vectorX);
	void SolveForwardLU(double* lowerTringMat, double* diag, double* rightVector, double* vectorX);
	void SolveBackwardLU(double* upperTringMat, double* diag, double* rightVector, double* vectorX);

	void MatrixUVectorMultiplicationLU(double* U, double* vectorMult, double* vectorOut);
	void CalculateZ(double* vectorOut);

	void VectorSubtract(double* first, double* second, double* result);
	double VectorScalarProduction(double* vector);
	double VectorScalarProduction(double* vector1, double* vector2);
	double VectorNorm(double* vector);
	void VectorCopy(double* first, double* second);
	void VectorZeroing(double* vector);

	void OutputDense();
	void OutputLUDense();


	void VectorOutputSolution(FILE* out);

protected:
	void AllocateMemory();
	void ClearMemory();
};