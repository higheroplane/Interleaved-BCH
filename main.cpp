#include "bch.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <ctime> 
#include <cstdlib>

#define N 100
#define M 7
#define K 80

int hamm_dist (int * w, int size)
{
    int weight = 0;
    for (int i = 0; i < size; i ++) if (w[i] != 0) weight ++; 
    return weight;
}


const int _MAX_ITER = 100;

int main ()
{
    BCH code (M, N, K);

    code.dump ();

    int errv[N] = {};

    for (int k = 0; k < 5; k ++) errv[rand()%N] = 1; 

    code.decoder(errv);
    //return 0;

    int stats [N+1][2] = {};

    srand(time(NULL));

    for (int i = 0; i < _MAX_ITER; i ++)
    {
        if (i % (_MAX_ITER / 20) == 0) printf ("%d\n", i);
        
        for (int j = 0; j < N/2; j ++) 
        {
            for (int k = 0; k < N; k ++) errv[k] = 0;
            for (int k = 0; k < j; k ++) errv[rand()%N] = 1; 
            int res = code.decoder(errv);
            int weight = hamm_dist(errv, N);
            stats [weight][res] ++;
        }
        for (int j = N/2; j <= N; j ++) 
        {
            for (int k = 0; k < N; k ++) errv[k] = 1;
            for (int k = 0; k < N - j; k ++) errv[rand()%N] = 0; 
            int res = code.decoder(errv);
            int weight = hamm_dist(errv, N);
            stats [weight][res] ++;
        }
        

    }

    FILE * out = fopen ("res.txt", "wb");
    for (int i = 0; i < N + 1; i ++) fprintf ( out, "%d %lf %lf \n", i,
                                     (double)stats[i][1] / ((double)stats[i][1] + (double)stats[i][0]), 
                                     (double)stats[i][0]/ ((double)stats[i][1] + (double)stats[i][0]));
    
    return 0;
}