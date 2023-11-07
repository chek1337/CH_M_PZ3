#include <iostream>
#include "slae.h"

int main()
{
    SLAE slae;
    FILE* paramf, * iaf, * jaf, * alf, * auf, * dif, * bf;
    fopen_s(&paramf, "param12.txt", "r");
    fopen_s(&iaf, "ia12.txt", "r");
    fopen_s(&jaf, "ja12.txt", "r");
    fopen_s(&alf, "al12.txt", "r");
    fopen_s(&auf, "au12.txt", "r");
    fopen_s(&dif, "di12.txt", "r");
    fopen_s(&bf, "b12.txt", "r");
    FILE* out;
    fopen_s(&out, "out.txt", "w");
    

    slae.Input(paramf, iaf, jaf, alf, auf, dif, bf);
    slae.OutputDense();
    //double* temp;
    //temp = new double[slae.n];
    //for (int i = 0; i < slae.n; i++)
    //{
    //    temp[i] = i + 1;
    //}
    //slae.TransposedMatrixVectorMultiplication(temp, slae.x);
    slae.MethodOfConjugateGradientsForNonSymMatrix();
    slae.VectorOutputSolution(out);
}