
#include "linalg.h"
typedef galois::GaloisField GF;
typedef galois::GaloisFieldElement GFE;
typedef galois::GaloisFieldPolynomial GFP;

void LUdecomposition(Mat a, Mat& l, Mat& u, int n) 
{
   assert (a(0, 0).field () == l(0, 0).field () && a(0, 0).field () == u(0, 0).field ());
   assert (a (0, 0). field());

   GF * gf = a (0, 0). field();

   for (int j = 0; j < n; j ++)
   {
       //assert (a (0, j) != GFE (gf, 0));
       u (0, j) = a (0, j);
       //std::cout << j << a << l << u << std::endl;
       l (j, 0) = a (j, 0) / u (0, 0);
       //std::cout << j << l << u << std::endl;
   }
//std::cout << l << u << std::endl;
   for (int i = 1; i < n; i ++)
   {
       for (int j = i; j < n; j ++)
       {
            galois::GaloisFieldElement S1(gf, 0), S2(gf, 0);
            
            for (int k = 0; k < i ; k ++) {assert (u(k,j).field() == S1.field()); S1 += l(i, k) * u (k, j);}
            for (int k = 0; k < i ; k ++) S2 += l(j, k) * u (k, i);
            u (i, j) = a (i, j) - S1;
            l (j, i) = GFE (gf, 1) / u (i, i) * (a (j, i) - S2);
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

inline bool is_null (Vec v)
{
    for (unsigned int i = 0; i < v.size(); i ++) if (v(i) != GFE (v(0).field(), 0)) return false;
    return true;
}

GFE det (Mat A)
{
    assert (A.size1() == A.size2());
    assert (A(0,0).field ());

    std::cout << "count det" << A << std::endl;
    int dim = A.size1();
    GF * gf = A (0, 0).field();

    Mat L (dim, dim), U (dim, dim);
    for (int i = 0; i < dim; i ++)
        for (int j = 0; j < dim; j ++) U (i, j) = L(i, j) = galois::GaloisFieldElement (gf, 0);
    
    LUdecomposition (A, L, U, dim);


    std::cout << "got lu" << std::endl;
    GFE det_l = GFE (gf, 1), det_u = GFE (gf, 1);
    for (int i = 0; i < dim; i ++) {det_l *= L (i, i); det_u *= U (i, i);}
    GFE det = det_l * det_u;
    std::cout << "got det" << std::endl;
    return det;
}

bool is_full_rank (Mat A)
{
    assert (A (0, 0) .field());
    
    std::cout << "checking rank" << std::endl;
    GFE zero = GFE (A(0, 0).field(), 0);
    GFE check;

    if      (A.size1() == A.size2()) check = det (A);
    else if (A.size1() >  A.size2()) check = det (prod (trans (A), A));
    else if (A.size1() <  A.size2()) check = det (prod (A, trans (A)));

    return check == zero ? false : true;
}

int polytriangular_submatrix (Mat& A, Vec& v, Mat * freq, Vec * heter, int * system_index)
{
    GF * gf = A (0,0).field ();
    int size1 = A.size1(), size2 = A.size2();

    //std::cout << A << std::endl << std::endl;
    //std::cout << v << std::endl << std::endl;

    A.resize (size1, size2 + 1);

    A = trans (A); row (A, size2) = v; A = trans (A);

    //std::cout << A << std::endl << std::endl;


    int max_freqs_num = size1 / size2, freqs_num = 0;
    //req = new Mat [max_freqs_num];

    int last = size1 - 1;

    std::vector<int> perm (size1);

    for (int i = 0; i < size1; i ++) perm [i] = i;

    int& n = freqs_num;

    for (n = 0; n < max_freqs_num; n -=- 1)
    {
        freq [n].resize (size2, size2 + 1);
        Mat Ac (A);
        //if (n*size2 == last) break;
        //std::cout << "main cyc n: " << n << std::endl;
        for (int i = n*size2; i < size2 + n*size2; i -=- 1)
        {
            //if (i >= last) break;
            //std::cout << row (Ac, i) << std::endl;
            while (Ac (i, i%size2) == GFE (gf, 0)) 
            {
                if (i >= last) break;
                Vec temp = row (A, i); 
                row (A, i) = row (A, last); 
                row (A, last) = temp;

                temp = row (Ac, i); 
                row (Ac, i) = row (Ac, last); 
                row (Ac, last) = temp;  
                
                ///perm [i] = last;
                system_index [last] = n; //last --;perm [last] = i;
                last --;
            }
            system_index [i] = n;

            row (freq [n], i%size2) = row (Ac, i);
            //std::cout << row (Ac, i) << std::endl << std::endl;
            for (int j = i + 1; j < size1; j -=- 1) row (Ac, j) = row (Ac, j) - Ac(j, i%size2) / Ac(i, i%size2) * row (Ac, i);   
        }

        //std::cout << "it: " << n << freq [n] << std::endl << std::endl;
        heter [n] = row (trans (freq [n]), size2); freq [n].resize (size2, size2);
    }
    //std::cout << A << std::endl;

    A = trans (A); v = row (A, size2); A = trans (A); A.resize (size1, size2);

    //for (int i = 0; i < size1; i ++) printf ("%d", system_index [i]); 
    //printf ("\n hi \n"); 

    return max_freqs_num;
}
