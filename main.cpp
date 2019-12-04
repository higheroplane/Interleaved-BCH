#include "bch.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <ctime> 
#include <cstdlib>


const int N = 31;
const int M = 5;
const int K = 21;
const int L = 70; 
const int T = (N - K)*L / (L + 1) * 2;

inline int hamm_dist (int * w, int size)
{
    int weight = 0;
    for (int i = 0; i < size; i ++) if (w[i] != 0) weight ++; 
    return weight;
}


const int _MAX_ITER = 10000;

int main ()
{
    BCH code (M, N, K);

    //code.dump ();

    //printf ("<<%d>>\n", T);

    int ** errv = new int* [L];
    for (int i = 0; i < L; i ++) errv [i] = (int*)calloc (N, sizeof (int));
    
    srand (time (NULL));
    int loc [2][T] = {};
    for (int k = 0; k < T; k ++) {loc [0][k] = rand() % N;}
    for (int k = 0; k < T; k ++) {loc [1][k] = rand() % N;}

    qsort(loc, T, sizeof(int), [](const void* a, const void * b) {return *static_cast<const int*>(a) - *static_cast<const int*>(b);});

    for (int i = 0; i < L; i ++)
    {
        for (int j = 0; j < T; j ++)
        {
            errv [i][loc[0][j]]= rand()%2;
        }
    }

    
    for (int i = 0; i < L; i ++)
    {
        int k = 0;
        for (int j = 0; j < N; j ++)
        {
            if (j == loc[i/(L/2)] [k] && k < T) 
            {
                //printf ("\033[31m%d\033[0m", errv [i][j] );
                while (loc[i%(L/2)][k] == loc [i/(L/2)][k + 1]) k ++; k++;
            }
            //else printf ("%d", errv [i][j]);
        }
       // printf ("\n");
    }


    //printf ("result of decoding is %d", code.collaborative_decoder (errv, L));
    
    return 0;
}

int rain ()
{
    BCH code (M, N, K);

    int errv[N] = {};
    int errv1[N] = {};

    errv [8] = errv[1] = errv[2] = errv [7] = errv [3] = 1;

    std::cout << code.syndrome(errv) << "\n";

    //errv [0] = errv[1] = errv[2] = errv [7] = errv [3] = 0;

    errv1 [8] = errv1[1] = errv1[2] = errv1 [5] = errv1 [6] = 1;



    std::cout << code.syndrome(errv1)  << "\n";

    galois::GaloisFieldElement prod [4];
    galois::GaloisFieldElement inner_prod;

    for (int i = 0; i < 4; i ++) {prod [i] = (code.syndrome(errv1).poly[i]^2) * (code.syndrome(errv).poly[i] ^2); inner_prod += prod[i];}

    galois::GaloisFieldPolynomial product (code.g->field(), 3, prod);

    std::cout << product << "\n";

    
    std::cout << inner_prod << "\n";
    

    errv [7] = errv [3] = 0;

    std::cout << code.syndrome(errv) << "\n";

    return 0;

}
