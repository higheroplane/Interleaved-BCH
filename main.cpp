#include "bch.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <ctime> 
#include <cstdlib>
#include "stoch.h"



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
inline llu binom (unsigned int n, unsigned int k)
{
    //printf ("bun : %d %d\n", n, k);
    //assert (n >= k || n == 0);
    if (k == n) return 1;
    if (k == 0 || k > n) return 0;
    llu b = n;
    for (unsigned int i = 1; i < k;  i ++) b *= (n - i);
    for (unsigned int i = 1; i <= k; i ++) b /= i;
    return b;
}

bool BigInt::resize_mode = true;

int main ()
{
    FILE * out = fopen ("out.txt", "wb");
    
    //printf ("%d", sizeof (unsigned long long int (9)));
    BigInt one (std::vector <int> ({0,0,0,0,1,0,0,0,0,0,1}));
    BigInt two (std::vector <int> ({1,1}));
    two = two * two;
    one = one * one;
    //one.print();
    //two.print();
    two = one / two;
    //two.print();
    //return 0;
    srand (time (NULL));
    int len=7, b = 3, c = 3;
    printf ("<<%d %d %d>>\n", len, b, c);
    for (int i = 1; i < 5; i ++)
    {
        //_BI(1+i).print();
        auto v = enumerative_decoder_mark (_BI(1), len, b, c);
        std::copy(begin(v), end(v), std::ostream_iterator<int>(std::cout, ""));
        printf (" %lld!\n", binom (2, 2));
    }
   // std::cout << numerical_encoder (10, 20, 7) << std::endl;
    
    
    //BCH code (M, N, 100);

    //BigInt one (std::vector <int> ({1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1}));
    //BigInt two (std::vector <int> ({0,0,0,0,0,0,0,0,0,0,0,0,1}));
    //two = one / two;
    //two.print();

    
    return 0;
}


