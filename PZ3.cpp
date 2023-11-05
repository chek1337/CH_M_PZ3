#include <iostream>
#include "slae.h"

int main()
{
    SLAE slae;
    FILE* paramf, * iaf, * jaf, * alf, * auf, * dif, * bf;
    fopen_s(&paramf, "param8.txt", "r");
    fopen_s(&iaf, "ia8.txt", "r");
    fopen_s(&jaf, "ja8.txt", "r");
    fopen_s(&alf, "al8.txt", "r");
    fopen_s(&auf, "au8.txt", "r");
    fopen_s(&dif, "di8.txt", "r");
    fopen_s(&bf, "b8.txt", "r");
    FILE* out;
    fopen_s(&out, "out.txt", "w");
    

    slae.Input(paramf, iaf, jaf, alf, auf, dif, bf);
    slae.OutputDense();
    slae.MethodOfConjugateGradients();
    slae.VectorOutputSolution(out);
}