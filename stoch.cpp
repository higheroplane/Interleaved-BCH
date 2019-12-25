
#include "stoch.h"

void make_hilbert_errors (std::vector<int*>& errv, int pgb, int pbg, int pe, int N)
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

inline int binom (int n, int k)
{
    assert (n >= k || n == 0);
    if (k == 0 || k == n) return 0;
    int b = n;
    for (int i = 1; i < k; i ++) b *= (n-i);
    for (int i = 1; i <= k; i ++) b /= i;
    return b;
}

std::vector<int> numerical_encoder_w (int n, int len, int weight)
{
    std::vector <int> ret (len, 0);
    for (int i = 0, k = 0; i < weight; i ++)
    {
        while (binom (len-k - 1, weight - i) > n) k++;
        n -= binom (len - k - 1, weight - i);
        ret [k] = 1;
        k ++;
    }

    return ret;
}
