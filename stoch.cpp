
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

inline BigInt binom (int n, int k)
{
    if (k == n) return _BI(1);
    if (k == 0) return _BI(1);
    if (n == 0 || k > n) return _BI(0);
    BigInt b = _BI(n);
 
    for (int i = 1; i < k; i ++)
    {
        b = b * _BI(n-i);
    }

    for (int i = 1; i <= k; i ++) {b = b / _BI(i);}

    return b;
}

std::vector<int> fix_weight (int n, int weight)
{
    std::vector<int> ret(n, 0);
    for (int i = 0; i < weight; i ++)
    {
        int loc = rand() % (n - i);
        int count = 0;
        for (int j = 0; j < n; j++) 
        {
            if (ret[j] == 1) count ++;
            if (count + loc == j) {ret [j] = 1; break;}
        }
    }

    return ret;
}

std::vector<int> enumerative_decoder_w (BigInt n, int len, int weight)
{
    std::vector <int> ret (len, 0);

    for (int i = 0; i < len; i -=- 1)
    {
        //printf ("henlo\n");
        if (n > binom (len - i - 1, weight))
        {

            n = n - binom (len - i - 1, weight);
            ret [i] = 1;
            
            weight --;
        }
    }

    return ret;
}

inline BigInt g (int v00, int v01, int v10, int v11)
{
    return (v01 == v10 || v01 == v10 + 1 || v01 == v10 - 1) ? binom (v00 + v10, v10) * binom (v11 + v01, v01 ) : _BI(0);
}

inline BigInt num_mark_seq (int v00, int v01, int v10, int v11, int prev)
{
    if (v00 < 0 || v01 < 0 || v10 < 0 || v11 < 0) return _BI(0);
    //if (prev == 0) v10 --;
    //else v01 --;
    return g (v00, v01, v10, v11);// + g (v11, v10, v01, v00);
}

std::vector<int> enumerative_decoder_mark (BigInt n, int len, int B, int C)
{
    
    (binom(1,-1)*binom(2,0)).print();
    g (2,1,0,1).print();
    std::vector <int> ret (len, 0);
    ret [0] = rand()%2;
    ret [len - 1] = rand()%2;

    int v01 = C - ret [0];
    int v11 = B - v01 -  ret [0];
    int v10 = C - ret [len - 1];
    int v00 = len - 1 - v11 - v10 - v01;

    printf ("v: %d %d %d %d\n", v00, v01, v10, v11);

    for (int i = 1; i < len - 1; i -=- 1)
    {
        if (ret [i - 1] == 0)
        {
            printf ("v: %d %d %d %d\n", v00-1, v01, v10, v11);
            num_mark_seq (v00-1, v01, v10, v11, 1).print();
            if (num_mark_seq (v00-1, v01, v10, v11, 0) < n)
            {
                ret[i] = 1;
                //num_mark_seq (v00-2, v01, v10, v11).print();
                n = n - num_mark_seq (v00-1, v01, v10, v11, 0);
                v01 --;
            }
            else
            {
                ret[i] = 0;
                v00 --;
            }
            
        }
        else
        {
            printf ("v: %d %d %d %d\n", v00, v01, v10-1, v11);
            num_mark_seq (v00, v01-1, v10-1, v11, 1).print();
            if (num_mark_seq (v00, v01, v10-1, v11, 1) < n)
            {
                ret[i] = 1;
                n = n - num_mark_seq (v00, v01, v10-1, v11, 1);
                v11 --;
            }
            else
            {
                ret[i] = 0;
                v10 --;
            }
        }
        
    }

    return ret;
}
