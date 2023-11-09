#include <iostream>
#include "slae.h"

int main()
{
    SLAE slae;
    FILE* paramf, * iaf, * jaf, * alf, * auf, * dif, * bf;
    fopen_s(&paramf, "paramLU.txt", "r");
    fopen_s(&iaf, "iaLU.txt", "r");
    fopen_s(&jaf, "jaLU.txt", "r");
    fopen_s(&alf, "alLU.txt", "r");
    fopen_s(&auf, "auLU.txt", "r");
    fopen_s(&dif, "diLU.txt", "r");
    fopen_s(&bf, "b12.txt", "r");
    FILE* out;
    fopen_s(&out, "out.txt", "w");
    

    slae.Input(paramf, iaf, jaf, alf, auf, dif, bf);
    slae.OutputDense();
    slae.CalculateLU();
    slae.OutputDense();
  

    //double* temp;
    //temp = new double[slae.n];
    //for (int i = 0; i < slae.n; i++)
    //{
    //    temp[i] = i + 1;
    //}
    //slae.TransposedMatrixVectorMultiplication(temp, slae.x);
    //slae.MethodOfConjugateGradientsForNonSymMatrix();
   // slae.MethodOfConjugateGradientsForNonSymMatrixWithDiagP();
    //slae.VectorOutputSolution(out);
}