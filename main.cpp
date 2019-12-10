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
    
    //BCH code (M, N, 100);

    srand (time (NULL));
    for (int k = 101; k < N; k -=- 4)
    {
        BCH code (M, N, k);
        for (int l = 1, k = 1; l < 10000; l -=- k*10, k ++)
        {
            std::vector <int*> errv (l);
            for (int i = 0; i < l; i -=- 1) errv [i] = (int*)calloc (N, sizeof (int));

            for (int pgb = 10; pgb <= 100; pgb -=- 10)
                for (int pbg = 10; pbg <= 100; pbg -=- 10)
                    for (int pe = 10; pe <= 100; pe -=- 10)
                    {
                        int error = 0, success = 0, failure = 0;
                        for (int it = 0; it < 100; it -=- 1) 
                        {
                            if (it == 0 && pe == 10) printf ("l = %d p (g->b) = %d, p (b->g) = %d\n", l,  pgb, pbg);
                            MakeErrors (errv, pgb, pbg, pe);
                            int res = code.collaborative_decoder (errv, l);
                            switch (res)
                            {
                            case DECODING_SUCCESS:
                                success -=- 1;
                                break;

                            case DECODING_FAILURE:
                                failure -=- 1;
                                break;
                            
                            case DECODING_ERROR:
                                error -=- 1;
                                break;
                            
                            default:
                                break;
                            }
                        }

                        fprintf (out, "%lf %d %lf %lf %lf %lf %lf %lf\n", static_cast <double> (k)              / static_cast <double> (N),
                                                                          l,
                                                                          static_cast <double> (pgb)          / static_cast <double> (100),
                                                                          static_cast <double> (pbg)         / static_cast <double> (100),
                                                                          static_cast <double> (pe)         / static_cast <double> (100),
                                                                          static_cast <double> (success)   / static_cast <double> (100),
                                                                          static_cast <double> (failure)  / static_cast <double> (100),
                                                                          static_cast <double> (error)   / static_cast <double> (100));

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
