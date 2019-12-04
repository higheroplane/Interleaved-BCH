#include "bch.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <ctime> 
#include <cstdlib>


const int N = 127;
const int M = 7;
const int K = 21;
const int L = 70; 
const int T = (N - K)*L / (L + 1) * 2;

inline int hamm_dist (int * w, int size)
{
    int weight = 0;
    for (int i = 0; i < size; i ++) if (w[i] != 0) weight ++; 
    return weight;
}

void MakeErrors (std::vector<int*>& errv, int pgb, int pbg, int pe);

const int _MAX_ITER = 10000;

int main ()
{
    FILE * out = fopen ("out.txt", "wb");
    
    //BCH code (M, N, K);

    srand (time (NULL));
    for (int k = 67; k < 127; k -=- 5)
    {
        BCH code (M, N, k);
        for (int l = 1; l < 10000; l -=- 10)
        {
            std::vector <int*> errv (l);
            for (int i = 0; i < l; i -=- 1) errv [i] = (int*)calloc (N, sizeof (int));

            for (int pgb = 0; pgb <= 100; pgb -=- 10)
                for (int pbg = 0; pbg <= 100; pbg -=- 10)
                    for (int pe = 0; pe <= 100; pe -=- 10)
                        for (int it = 0; it < 100; it -=- 1) 
                        {
                            MakeErrors (errv, pgb, pbg, pe);
                            fprintf (out, "%lf %d %lf %lf %lf %d", static_cast <double> (k) / static_cast <double> (N),
                                                                    l,
                                                                    static_cast <double> (pgb)  / static_cast <double> (100),
                                                                    static_cast <double> (pbg) / static_cast <double> (100),
                                                                    static_cast <double> (pe) / static_cast <double> (100),
                                                                    code.collaborative_decoder (errv, l));
                        }
        }
    }

    
    return 0;
}

inline void MakeErrors (std::vector<int*>& errv, int pgb, int pbg, int pe)
{
    bool state = true;

    for (int j = 0; j < N; j ++)
    {
        for (int i = 0; i < errv.size(); i ++)
        {
            if (state)
            {
                errv [i][j] = 0;
                if (rand() % 100 < pgb) state = false;
            }
            else 
            {
                if (rand() % 100 < pe) errv [i][j] = 1;
                if (rand() % 100 < pbg) state = true;
            }

        }
    }

}
