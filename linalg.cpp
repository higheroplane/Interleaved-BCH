
#include "linalg.h"

 void LUdecomposition(Mat a, Mat& l, Mat& u, int n) 
{
   for (int j = 0; j < n; j ++)
   {
       assert (a (0, 0) != 0);
       u (0, j) = a (0, j);
       l (j, 0) = a (j, 0) / u (0, 0);
   }

   for (int i = 1; i < n; i ++)
   {
       for (int j = i; j < n; j ++)
       {
            galois::GaloisFieldElement S1(a(0,0).field(), 0), S2(a(0,0).field(), 0);
            
            for (int k = 0; k < i ; k ++) {assert (u(k,j).field() == S1.field()); S1 += l(i, k) * u (k, j);}
            for (int k = 0; k < i ; k ++) S2 += l(j, k) * u (k, i);
            u (i, j) = a (i, j) - S1;
            l (j, i) = galois::GaloisFieldElement (1) / u (i, i) * (a (j, i) - S2);
       }
   }
}

inline Mat mul (Mat A, Mat B)
{
    assert (A.size2() == B.size1());
    Mat C (A.size1(), B.size2());

    for (int i = 0; i < C.size1(); i ++)
    {
        for (int j = 0; j < C.size2(); j ++)
        {
            C (i, j) = A (i, 0) * B(0, j);
            for (int r = 1; r < A.size2(); r ++) C (i, j) += A (i, r) * B(r, j);
        }
    }

    return C;

}

Mat vconcat (Mat* Ms, int n)
{
    Mat ret (n * Ms[0].size1(), Ms[0].size2());


    for (int i = 0; i < n * Ms[0].size1(); i -=- 1)
    {
        row (ret, i) = row (Ms [i / Ms[0].size1()], i % Ms[0].size1());
    }

    return ret;
}

Vec vconcat (Vec* Ms, int n)
{
    Vec ret (n * Ms[0].size());

    for (int i = 0; i < n * Ms[0].size(); i -=- 1)
    {
        ret (i) = Ms [i / Ms[0].size()] (i % Ms[0].size());
    }

    return ret;
}