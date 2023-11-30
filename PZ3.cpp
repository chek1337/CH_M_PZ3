#include <iostream>
#include <chrono>
#include "slae.h"

using namespace std::chrono;
using namespace std;
//Не забыть вынести сделать функцию для подсчета относитльной невязки (F-Ax)


int main()
{
    setlocale(LC_ALL, "Russian");
    
    int size;
    SLAE slae;
    FILE* paramf, * iaf, * jaf, * alf, * auf, * dif, * bf;
    FILE* out;
    fopen_s(&out, "out.txt", "w");
    double seconds;
    int menu;
    auto start_time = steady_clock::now();
    auto end_time = steady_clock::now();

    printf("Введите пункт меню:\n");
    printf("1) МСГ для симметричной матрицы без обуславливания\n");
    printf("2) МСГ для симметричной матрицы с диагональным обуславливанием\n");
    printf("3) МСГ для симметричной матрицы с LU обуславливанием\n\n");

    printf("4) МСГ для несимметричной матрицы без обуславливания \n");
    printf("5) МСГ для несимметричной матрицы с диагональным обуславливанием \n");
    printf("6) МСГ для несимметричной матрицы с LU обуславливанием \n");
    printf("7) МСГ для несимметричной матрицы с LU* обуславливанием \n");
    printf("8) МСГ для несимметричной матрицы с LU(sq) обуславливанием \n\n");

    printf("9) МСГ для матрицы Гильберта без обуславливания \n");
    printf("10) МСГ для матрицы Гильберта с диагональным обуславливанием \n");
    printf("11) МСГ для матрицы Гильберта с LU обуславливанием \n");
    printf("12) МСГ для матрицы Гильберта с LU* обуславливанием \n");
    printf("13) МСГ для матрицы Гильберта с LU(sq) обуславливанием \n");
    scanf_s("%d", &menu);
    setlocale(LC_ALL, "US");
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
        fopen_s(&paramf, "param4545.txt", "r");
        fopen_s(&iaf, "ia4545.txt", "r");
        fopen_s(&jaf, "ja4545.txt", "r");
        fopen_s(&alf, "al4545.txt", "r");
        fopen_s(&auf, "au4545.txt", "r");
        fopen_s(&dif, "di4545.txt", "r");
        fopen_s(&bf, "b4545.txt", "r");
        slae.Input(paramf, iaf, jaf, alf, auf, dif, bf);
        //slae.OutputDense();

        start_time = steady_clock::now();
        slae.MethodOfConjugateGradientsForSymMatrixWithDiagP();
        end_time = steady_clock::now();
        cout << duration_cast<microseconds>(end_time - start_time).count();

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

        start_time = steady_clock::now();
        slae.MethodOfConjugateGradientsForSymMatrixWithLuP();
        end_time = steady_clock::now();
        cout << duration_cast<microseconds>(end_time - start_time).count();

        slae.VectorOutputSolution(out);
        break;
    case 4:
        fopen_s(&paramf, "param945.txt", "r");
        fopen_s(&iaf, "ia945.txt", "r");
        fopen_s(&jaf, "ja945.txt", "r");
        fopen_s(&alf, "al945.txt", "r");
        fopen_s(&auf, "au945.txt", "r");
        fopen_s(&dif, "di945.txt", "r");
        fopen_s(&bf, "b945.txt", "r");
        slae.Input(paramf, iaf, jaf, alf, auf, dif, bf);
       //slae.OutputDense();

        start_time = steady_clock::now();
        slae.MethodOfConjugateGradientsForNonSymMatrix();
        end_time = steady_clock::now();
        cout << duration_cast<milliseconds>(end_time - start_time).count();

        slae.VectorOutputSolution(out);
        break;
    case 5:
        fopen_s(&paramf, "param945.txt", "r");
        fopen_s(&iaf, "ia945.txt", "r");
        fopen_s(&jaf, "ja945.txt", "r");
        fopen_s(&alf, "al945.txt", "r");
        fopen_s(&auf, "au945.txt", "r");
        fopen_s(&dif, "di945.txt", "r");
        fopen_s(&bf, "b945.txt", "r");
        slae.Input(paramf, iaf, jaf, alf, auf, dif, bf);
        //slae.OutputDense();

        start_time = steady_clock::now();
        slae.MethodOfConjugateGradientsForNonSymMatrixWithDiagP();
        end_time = steady_clock::now();
        cout << duration_cast<milliseconds>(end_time - start_time).count();

        slae.VectorOutputSolution(out);
        break;
    case 6:
        fopen_s(&paramf, "param6.txt", "r");
        fopen_s(&iaf, "ia6.txt", "r");
        fopen_s(&jaf, "ja6.txt", "r");
        fopen_s(&alf, "al6.txt", "r");
        fopen_s(&auf, "au6.txt", "r");
        fopen_s(&dif, "di6.txt", "r");
        fopen_s(&bf, "b6.txt", "r");
        slae.Input(paramf, iaf, jaf, alf, auf, dif, bf);
        //slae.OutputDense();
        slae.CalculateLU();
        slae.OutputLUDense();

        start_time = steady_clock::now();
        slae.MethodOfConjugateGradientsForNonSymMatrixWithLuP();
        end_time = steady_clock::now();
        cout << duration_cast<milliseconds>(end_time - start_time).count();

        slae.VectorOutputSolution(out);
        break;
    case 7:
        fopen_s(&paramf, "param6.txt", "r");
        fopen_s(&iaf, "ia6.txt", "r");
        fopen_s(&jaf, "ja6.txt", "r");
        fopen_s(&alf, "al6.txt", "r");
        fopen_s(&auf, "au6.txt", "r");
        fopen_s(&dif, "di6.txt", "r");
        fopen_s(&bf, "b6.txt", "r");
        slae.Input(paramf, iaf, jaf, alf, auf, dif, bf);
        slae.OutputDense();
        slae.CalculateLUaster();
        slae.OutputLUDense();

        start_time = steady_clock::now();
        slae.MethodOfConjugateGradientsForNonSymMatrixWithLuAsterP();
        end_time = steady_clock::now();
        cout << duration_cast<milliseconds>(end_time - start_time).count();

        slae.VectorOutputSolution(out);
        break;
    case 8:
        fopen_s(&paramf, "param945.txt", "r");
        fopen_s(&iaf, "ia945.txt", "r");
        fopen_s(&jaf, "ja945.txt", "r");
        fopen_s(&alf, "al945.txt", "r");
        fopen_s(&auf, "au945.txt", "r");
        fopen_s(&dif, "di945.txt", "r");
        fopen_s(&bf, "b945.txt", "r");
        slae.Input(paramf, iaf, jaf, alf, auf, dif, bf);
        //slae.OutputDense();
        slae.CalculateLUsq();
        //slae.OutputLUDense();

        start_time = steady_clock::now();
        slae.MethodOfConjugateGradientsForNonSymMatrixWithLuSqP();
        end_time = steady_clock::now();
        cout << duration_cast<milliseconds>(end_time - start_time).count();

        slae.VectorOutputSolution(out);
        break;

    case 9:

        slae.GenerateHilbertMatrix(50);
        //slae.OutputDense();
        //slae.OutputLUDense();
        start_time = steady_clock::now();
        slae.MethodOfConjugateGradientsForNonSymMatrix();
        end_time = steady_clock::now();
        cout << duration_cast<milliseconds>(end_time - start_time).count();

        slae.VectorOutputSolution(out);
        break;

    case 10:

        slae.GenerateHilbertMatrix(10);
        //slae.OutputDense();
        //slae.OutputLUDense();
        start_time = steady_clock::now();
        slae.MethodOfConjugateGradientsForNonSymMatrixWithDiagP();
        end_time = steady_clock::now();
        cout << duration_cast<milliseconds>(end_time - start_time).count();

        slae.VectorOutputSolution(out);
        break;
    case 11:

        slae.GenerateHilbertMatrix(12);
        //slae.OutputDense();
        //slae.OutputLUDense();

        slae.CalculateLU();
        start_time = steady_clock::now();
        slae.MethodOfConjugateGradientsForNonSymMatrixWithLuP();
        end_time = steady_clock::now();
        cout << duration_cast<milliseconds>(end_time - start_time).count();

        slae.VectorOutputSolution(out);
        break;

    case 12:

        slae.GenerateHilbertMatrix(8);
        slae.OutputDense();
        

        slae.CalculateLUaster();
        slae.OutputLUDense();
        start_time = steady_clock::now();
        slae.MethodOfConjugateGradientsForNonSymMatrixWithLuAsterP();
        end_time = steady_clock::now();
        cout << duration_cast<milliseconds>(end_time - start_time).count();

        slae.VectorOutputSolution(out);
        break;
    case 13:

        slae.GenerateHilbertMatrix(10);
        //slae.OutputDense();
        //slae.OutputLUDense();
        slae.CalculateLUsq();
        start_time = steady_clock::now();
        slae.MethodOfConjugateGradientsForNonSymMatrixWithLuSqP();
        end_time = steady_clock::now();
        cout << duration_cast<milliseconds>(end_time - start_time).count();

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