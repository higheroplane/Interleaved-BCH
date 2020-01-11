
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
    //printf ("bun : %d %d = ", n, k);
    //assert (n >= k || n == 0);
    if (k == n) return _BI(1);
    if (k == 0 || k > n) return _BI(0);
    BigInt b = _BI(n);
    //b.print();
 
    for (int i = 1; i < k; i ++)
    {
        //(b/_BI(2)).print();
       // (b * _BI(n-i)).print(); 
        b = b * _BI(n-i);
       // _BI(n-i).print();
        //b.print();
    }
//printf ("gg\n");

    for (int i = 1; i <= k; i ++) {b = b / _BI(i);}

    //std::cout <<b.size() << std::endl;
    return b;
}

std::vector<int> enumerative_decoder_w (BigInt n, int len, int weight)
{
    std::vector <int> ret (len, 0);

    for (int i = 0; i < len; i -=- 1)
    {
        //printf ("henlo\n");
        if (n > binom (len - i - 1, weight))
        {
            //printf ("\thenlo\n");
            //printf ("ret[%d] = %d\n", i, ret[i]);
            n = n - binom (len - i - 1, weight);
            ret [i] = 1;
            
            weight --;
        }
    }

    return ret;
}

inline BigInt g (int v00, int v01, int v10, int v11)
{
    return (v01 == v10 || v01 == v10 + 1) ? binom (v00 + v10 - 1, v10 - 1) * binom (v11 + v01 , v01 ) : _BI(0);
}

inline BigInt num_mark_seq (int v00, int v01, int v10, int v11)
{
    return g (v00, v01, v10, v11) ;//+ g (v11, v10, v01, v00);
}

std::vector<int> enumerative_decoder_mark (BigInt n, int len, int B, int C)
{
    
    std::vector <int> ret (len, 0);
    ret [0] = 0;//rand()%2;
    ret [len - 1] = 0;//rand()%2;

    int v01 = C - ret [0];
    int v11 = B - v01 ;
    int v10 = C - ret [len - 1];
    int v00 = len - 3 - v11 - v10 - v01;

    //printf ("v: %d %d %d %d\n", v00, v01, v10, v11);

    for (int i = 1; i < len - 1; i -=- 1)
    {
        if (ret [i - 1] == 0)
        {
            if (num_mark_seq (v00-1, v01, v10, v11) < n)
            {
                ret[i] = 1;
                n = n - num_mark_seq (v00-1, v01, v10, v11);
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
            if (num_mark_seq (v00, v01, v10-1, v11) < n)
            {
                ret[i] = 1;
                n = n - num_mark_seq (v00, v01, v10-1, v11);
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
