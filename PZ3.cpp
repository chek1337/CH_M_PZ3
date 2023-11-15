#include <iostream>
#include "slae.h"

//Не забыть вынести сделать функцию для подсчета относитльной невязки (F-Ax)


int main()
{
    setlocale(LC_ALL, "Russian");
    int size;
    SLAE slae;
    FILE* paramf, * iaf, * jaf, * alf, * auf, * dif, * bf;
    FILE* out;
    fopen_s(&out, "out.txt", "w");

    int menu;
    printf("Введите пункт меню:\n");
    printf("1) МСГ для симметричной матрицы без обуславливания\n");
    printf("2) МСГ для симметричной матрицы с диагональным обуславливанием\n");
    printf("3) МСГ для симметричной матрицы с LU обуславливанием\n\n");

    printf("4) МСГ для несимметричной матрицы без обуславливания (решение AtA)\n");
    printf("5) МСГ для несимметричной матрицы с диагональным обуславливанием (решение AtA)\n");
    printf("6) МСГ для несимметричной матрицы с LU обуславливанием (решение AtA)\n\n");
    scanf_s("%d", &menu);
    switch (menu)
    {
    case 1:
        fopen_s(&paramf, "param8.txt", "r");
        fopen_s(&iaf, "ia8.txt", "r");
        fopen_s(&jaf, "ja8.txt", "r");
        fopen_s(&alf, "al8.txt", "r");
        fopen_s(&auf, "au8.txt", "r");
        fopen_s(&dif, "di8.txt", "r");
        fopen_s(&bf, "b8.txt", "r");
        slae.Input(paramf, iaf, jaf, alf, auf, dif, bf);
        slae.OutputDense();
        slae.MethodOfConjugateGradientsForSymMatrix();
        slae.VectorOutputSolution(out);
        break;
    case 2:
        fopen_s(&paramf, "param8.txt", "r");
        fopen_s(&iaf, "ia8.txt", "r");
        fopen_s(&jaf, "ja8.txt", "r");
        fopen_s(&alf, "al8.txt", "r");
        fopen_s(&auf, "au8.txt", "r");
        fopen_s(&dif, "di8.txt", "r");
        fopen_s(&bf, "b8.txt", "r");
        slae.Input(paramf, iaf, jaf, alf, auf, dif, bf);
        slae.OutputDense();
        slae.MethodOfConjugateGradientsForSymMatrixWithDiagP();
        slae.VectorOutputSolution(out);
        break;
    case 3:
        fopen_s(&paramf, "param8.txt", "r");
        fopen_s(&iaf, "ia8.txt", "r");
        fopen_s(&jaf, "ja8.txt", "r");
        fopen_s(&alf, "al8.txt", "r");
        fopen_s(&auf, "au8.txt", "r");
        fopen_s(&dif, "di8.txt", "r");
        fopen_s(&bf, "b8.txt", "r");
        slae.Input(paramf, iaf, jaf, alf, auf, dif, bf);
        slae.OutputDense();
        slae.CalculateLU();
        slae.MethodOfConjugateGradientsForSymMatrixWithLuP();
        slae.VectorOutputSolution(out);
        break;
    case 4:
        fopen_s(&paramf, "param12.txt", "r");
        fopen_s(&iaf, "ia12.txt", "r");
        fopen_s(&jaf, "ja12.txt", "r");
        fopen_s(&alf, "al12.txt", "r");
        fopen_s(&auf, "au12.txt", "r");
        fopen_s(&dif, "di12.txt", "r");
        fopen_s(&bf, "b12.txt", "r");
        slae.Input(paramf, iaf, jaf, alf, auf, dif, bf);
        slae.OutputDense();
        slae.MethodOfConjugateGradientsForNonSymMatrixAtA();
        slae.VectorOutputSolution(out);
        break;
    case 5:
        fopen_s(&paramf, "param12.txt", "r");
        fopen_s(&iaf, "ia12.txt", "r");
        fopen_s(&jaf, "ja12.txt", "r");
        fopen_s(&alf, "al12.txt", "r");
        fopen_s(&auf, "au12.txt", "r");
        fopen_s(&dif, "di12.txt", "r");
        fopen_s(&bf, "b12.txt", "r");
        slae.Input(paramf, iaf, jaf, alf, auf, dif, bf);
        slae.OutputDense();
        slae.MethodOfConjugateGradientsForNonSymMatrixAtAWithDiagP();
        slae.VectorOutputSolution(out);
        break;
    case 6:
        fopen_s(&paramf, "param12.txt", "r");
        fopen_s(&iaf, "ia12.txt", "r");
        fopen_s(&jaf, "ja12.txt", "r");
        fopen_s(&alf, "al12.txt", "r");
        fopen_s(&auf, "au12.txt", "r");
        fopen_s(&dif, "di12.txt", "r");
        fopen_s(&bf, "b12.txt", "r");
        slae.Input(paramf, iaf, jaf, alf, auf, dif, bf);
        slae.OutputDense();
        slae.CalculateLU();
        slae.OutputLUDense();
        slae.MethodOfConjugateGradientsForNonSymMatrixAtAWithLuP();
        slae.VectorOutputSolution(out);
        break;
    default:
        printf("Неправильный пункт меню!");
        break;
    }

    //fopen_s(&paramf, "param12.txt", "r");
    //fopen_s(&iaf, "ia12.txt", "r");
    //fopen_s(&jaf, "ja12.txt", "r");
    //fopen_s(&alf, "al12.txt", "r");
    //fopen_s(&auf, "au12.txt", "r");
    //fopen_s(&dif, "di12.txt", "r");
    //fopen_s(&bf, "b12.txt", "r");

    //slae.Input(paramf, iaf, jaf, alf, auf, dif, bf);
    //slae.MethodOfConjugateGradientsForNonSymMatrixAtA();
    //slae.MethodOfConjugateGradientsForNonSymMatrixAtAWithDiagP();
    //slae.MethodOfConjugateGradientsForNonSymMatrixAtA();
    //slae.VectorOutputSolution(out);

    //fopen_s(&paramf, "paramLU.txt", "r");
    //fopen_s(&iaf, "iaLU.txt", "r");
    //fopen_s(&jaf, "jaLU.txt", "r");
    //fopen_s(&alf, "alLU.txt", "r");
    //fopen_s(&auf, "auLU.txt", "r");
    //fopen_s(&dif, "diLU.txt", "r");
    //fopen_s(&bf, "bLU.txt", "r");
    //slae.Input(paramf, iaf, jaf, alf, auf, dif, bf);
    //slae.OutputDense();
    //slae.CalculateLU();
    //slae.OutputLUDense();
}