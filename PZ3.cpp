#include <iostream>
#include "slae.h"

int main()
{
    SLAE slae;
    FILE* paramf, * iaf, * jaf, * alf, * auf, * dif, * bf;
    fopen_s(&paramf, "param.txt", "r");
    fopen_s(&iaf, "ia.txt", "r");
    fopen_s(&jaf, "ja.txt", "r");
    fopen_s(&alf, "al.txt", "r");
    fopen_s(&auf, "au.txt", "r");
    fopen_s(&dif, "di.txt", "r");
    fopen_s(&bf, "b.txt", "r");
    FILE* out;
    fopen_s(&out, "out.txt", "w");
    

    slae.Input(paramf, iaf, jaf, alf, auf, dif, bf);
    slae.OutputDense();
    slae.CalculateLU();
    slae.OutputLUDense();
    slae.MethodOfConjugateGradientsForNonSymMatrix();
    //slae.MethodOfConjugateGradientsForNonSymMatrixWithLuP();
    slae.VectorOutputSolution(out);
  

    //double* temp;
    //temp = new double[slae.n];
    //for (int i = 0; i < slae.n; i++)
    //{
    //    temp[i] = i + 1;
    //}
}