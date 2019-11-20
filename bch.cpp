#include <math.h>
#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include "sourceCode/galois/GaloisField.h"
#include "sourceCode/galois/GaloisFieldElement.h"
#include "sourceCode/galois/GaloisFieldPolynomial.h"
#include "bch.h"

BCH::BCH (int _m, int _n, int _k): m(_m), n(_n), k(_k), d(n - k + 1) 
{
    
    assert (d < n/2);
    //printf ("About to start creating BCH (%d, %d, %d) code space.\n", n, k, d);
    unsigned int * p = (unsigned int *) calloc (m+1, sizeof (int));

    p[0] = p[m] = 1;

    switch (m)
    {
        case 2:  p[1] = 1;                 break;
        case 3:  p[1] = 1;                 break;
        case 4:  p[1] = 1;                 break;
        case 5:  p[2] = 1;                 break;
        case 6:  p[1] = 1;                 break;
        case 7:  p[1] = 1;                 break;
        case 8:  p[2] = p[3] = p[4] = 1;   break;
        case 9:  p[4] = 1;                 break;
        case 10: p[3] = 1;                 break;
        case 11: p[2] = 1;                 break;
        case 12: p[1] = p[4] = p[6] = 1;   break;
        case 13: p[1] = p[3] = p[4] = 1;   break;
        case 14: p[1] = p[6] = p[8] = 1; break;
        case 15: p[1] = 1;                 break;
        case 16: p[1] = p[3] = p[12] = 1;   break;
        case 17: p[3] = 1;                 break;
        case 18: p[7] = 1;                 break;
        case 19: p[1] = p[2] = p[6] = 1;   break;
        case 20: p[3] = 1;                 break;
        default:
            printf ("Inappropriate size of field.\n");
            break;
    }

    //printf("p(x) = ");

    int l = 1;
	for (int i = 0; i <= m; i++) 
    {
        l *= 2;
	//	printf("%1d", p[i]);
    }

	l = l / 2 - 1;
	int L = (l + 1) / 2 - 1;
    //printf("\n %d %d %d\n", l, n, L);
    assert ((n <= l)&&(n > L));

    gf = new galois::GaloisField (m, p);
    a  = new galois::GaloisFieldElement(gf, gf -> alpha(1));
    g  = new galois::GaloisFieldPolynomial;

    gen_poly ();

    //std::cout << "Generator: " << g -> deg() << "\n";    

    //printf ("BCH (%d, %d, %d) seems to be built", n, k, d);

}


int * BCH::gen_poly()

{
    
    galois::GaloisFieldPolynomial M;
    galois::GaloisFieldPolynomial G;
    M = min_poly (*a);
    *g = M;
    for (int i = 3; i < d; i += 2) *g *= min_poly ((*a)^i);
}

galois::GaloisFieldPolynomial min_poly (galois::GaloisFieldElement a)
{
    galois::GaloisFieldElement temp [2] = {a, galois::GaloisFieldElement(a.field(), 1)};
    galois::GaloisFieldPolynomial M (a.field(), 1, temp);
    galois::GaloisFieldPolynomial G;

    for (int i = 1; (temp[0]^2) != a; i ++) 
    {
        temp [0] = temp[0]^2;
        G = galois::GaloisFieldPolynomial (a.field(), 1, temp);
        M = M * G;
    }

    //std::cout << "M: " << " = " << M.deg() << "\n";
    
    return M ;
}

int BCH::decoder (int * errors)
{
    assert (errors);
    int t = d/2;

    // Reading error vector
    galois::GaloisFieldElement * er = new galois::GaloisFieldElement [n];
    assert (er);
    //printf ("decoding HENLO!%d %d %d!\n", n, k, t);
    for (int i = 0; i < n; i ++) 
    {
        //printf ("it %d init \n", i);
        if (errors [i] == 1) er [i] = galois::GaloisFieldElement(gf, 1);
        else                 er [i] = galois::GaloisFieldElement(gf, 0);
    }
    galois::GaloisFieldPolynomial e (gf, n - 1, er);

    //std::cout << "Recieved word = " << e << "\n";


    // Preparing roots array
   // galois::GaloisFieldElement * roots = new galois::GaloisFieldElement [2*t];
    //for (int i = 1; i <= 2*t; i ++) {roots [i-1] = (*a)^i; std::cout << ((*a)^i).index() << i << "\n";}


    // Calculating syndrome
    galois::GaloisFieldElement * S = new galois::GaloisFieldElement [gf -> size()];
    for (int i = 0; i < gf -> size(); i ++) S[i] = galois::GaloisFieldElement (gf, 0);
    for (int i = 0; i < 2*t; i ++) {S[i] = e ((*a)^(i+1)); 
    /*std::cout << "S[" << roots[i].index() << "] = " << S[roots[i].index()] << "\n";*/}
    galois::GaloisFieldPolynomial Sz  (gf, gf -> size()-1, S);
    Sz.simplify();


    // Preparing p(z) = z^d
    galois::GaloisFieldElement * zd  = new galois::GaloisFieldElement [d+1];
    for (int i = 0; i < d+1 ; i ++) 
      {
        if (i == d-1)      zd[i] = galois::GaloisFieldElement (gf, 1);
        else               zd[i] = galois::GaloisFieldElement (gf, 0);
      }
    galois::GaloisFieldPolynomial z_d (gf, d, zd);


    // Calculating ELP
    galois::GaloisFieldPolynomial sigma(gf, 0);
    sigma = galois::key_equation (Sz, z_d, d);
    //std::cout << sigma << "\n";
    
    // Recovering the word
    int error_num = 0;{}
    //for (int i = 0; i < n; i ++) std::cout << '0';
    //printf ("\n recieved: \n");
    //for (int i = 0; i < n; i ++) std::cout << e [i];
    //printf ("\n decoded: \n");
    for (int i = 0; i < n; i ++)
    {
        if (sigma ((*a)^(gf -> size()  - i)) == galois::GaloisFieldElement (gf, 0)) {/*std::cout << "e";*/ error_num ++;}
        //else std::cout << e [i]; 
    }

    if (error_num != sigma.deg()) {/*printf ("decoding failure");*/ return 1;}
    else return 0;



}

galois::GaloisFieldPolynomial BCH::syndrome (int * errors)
{

    // Reading error vector
    galois::GaloisFieldElement * er = new galois::GaloisFieldElement [n];
    assert (er);
    //printf ("decoding HENLO!%d %d %d!\n", n, k, t);
    for (int i = 0; i < n; i ++) 
    {
        //printf ("it %d init \n", i);
        if (errors [i] == 1) er [i] = galois::GaloisFieldElement(gf, 1);
        else                 er [i] = galois::GaloisFieldElement(gf, 0);
    }
    galois::GaloisFieldPolynomial e (gf, n - 1, er);

    //std::cout << "Recieved word = " << e << "\n";


    // Preparing roots array
   // galois::GaloisFieldElement * roots = new galois::GaloisFieldElement [2*t];
    //for (int i = 1; i <= 2*t; i ++) {roots [i-1] = (*a)^i; std::cout << ((*a)^i).index() << i << "\n";}


    // Calculating syndrome
    galois::GaloisFieldElement * S = new galois::GaloisFieldElement [gf -> size()];
    for (int i = 0; i < gf -> size(); i ++) S[i] = galois::GaloisFieldElement (gf, 0);
    for (int i = 0; i < 2 * (d/2); i ++) {S[i] = e ((*a)^(i+1)); 
    /*std::cout << "S[" << roots[i].index() << "] = " << S[roots[i].index()] << "\n";*/}
    galois::GaloisFieldPolynomial Sz  (gf, gf -> size()-1, S);
    Sz.simplify();

    return Sz;

}

void BCH::dump ()
{
    printf ("This is BCH (%d, %d, %d)\n", n, k, d);

    std::cout << "g(x) =" << *g << "\n";
    std::cout << "alpha =" << *a << "\n";

}



int BCH::collaborative_decoder (int ** error, int l)
{
    int t = l*(n-k) / (l + 1);

    galois::GaloisFieldPolynomial * S = new galois::GaloisFieldPolynomial [l];
    for (int i = 0; i < l; i ++) S[i] = syndrome (error[i]);

    galois::GaloisFieldElement z[1] = {galois::GaloisFieldElement (gf, 0)};
    galois::GaloisFieldElement o[1] = {galois::GaloisFieldElement (gf, 1)};
    galois::GaloisFieldElement x[2] = {galois::GaloisFieldElement (gf, 0), galois::GaloisFieldElement (gf, 1)};
    galois::GaloisFieldPolynomial zero (gf, 0, z);
    galois::GaloisFieldPolynomial one  (gf, 0, o);
    galois::GaloisFieldPolynomial xx  (gf, 1, x);

    boost::numeric::ublas::matrix<galois::GaloisFieldElement> * SS = new boost::numeric::ublas::matrix<galois::GaloisFieldElement> [l];
    boost::numeric::ublas::vector<galois::GaloisFieldElement> * ss = new boost::numeric::ublas::vector<galois::GaloisFieldElement> [l];

    printf ("t = %d\n", t);

    for (int i = 0; i < l; i -=- 1) 
    {
        SS [i] .resize (n - k - t, t);
        ss [i] .resize (n - k - t);
        for (int j = 0; j < SS[i].size1(); j ++)
        {
            if (S[i].deg () >= j + t) ss[i](j) = S[i][j + t];
            else ss[i] (j) = galois::GaloisFieldElement (gf, 0);

            for (int k = 0; k < SS[i].size2(); k ++ )  
                if (S[i].deg () >= j + k) SS[i] (j, k) = S[i][j + k];
                else SS[i] (j, k) = galois::GaloisFieldElement (gf, 0);
        }
    }
    std::cout  <<  ss[0] << std::endl;
    std::cout  <<  SS[0] << std::endl;

    boost::numeric::ublas::vector<galois::GaloisFieldElement> Lambda (t);

    Lambda = boost::numeric::ublas::solve (SS, ss);

}   

int BCH::collaborative_decoder_srl (int ** error, int l)
{
        printf ("decoding henlo!\n");
    galois::GaloisFieldPolynomial * S = new galois::GaloisFieldPolynomial [l];
    for (int i = 0; i < l; i ++) S[i] = syndrome (error[i]);

    galois::GaloisFieldElement z[1] = {galois::GaloisFieldElement (gf, 0)};
    galois::GaloisFieldElement o[1] = {galois::GaloisFieldElement (gf, 1)};
    galois::GaloisFieldElement x[2] = {galois::GaloisFieldElement (gf, 0), galois::GaloisFieldElement (gf, 1)};
    galois::GaloisFieldPolynomial zero (gf, 0, z);
    galois::GaloisFieldPolynomial one  (gf, 0, o);
    galois::GaloisFieldPolynomial xx  (gf, 1, x);

    int te = 0;
    galois::GaloisFieldPolynomial elp(one);
    galois::GaloisFieldElement delta (gf, 0);

    int * t_l = new int [l];
    int * ms = (int*) calloc (l, sizeof (int));
    galois::GaloisFieldPolynomial * elps = new galois::GaloisFieldPolynomial [l];
    galois::GaloisFieldElement * deltas = new galois::GaloisFieldElement [l];
    galois::GaloisFieldElement * xm = new galois::GaloisFieldElement [m];
    for (int i = 0; i < l; i ++) {elps [i] = zero; t_l[i] = 0; deltas[i] = galois::GaloisFieldElement (gf, 1);}
    galois::GaloisFieldPolynomial x_m (one);

    //printf ("cycle starr\n");   
    for (int mi = 0; mi < n - k; mi ++)
    {
        
 //getchar();
        for (int li = 0; li < l; li ++)
        {
            std::cout << "elp:" << elp << " t: " << te << "\n\n";
            
            if (mi - te > 0)
            {
                galois::GaloisFieldElement sum (gf, 0);
                for (int i = 0; i < te; i ++) {if (elp.deg() > i && S[li].deg() > mi - i - 1) sum += elp[i]*S[li][mi - i - 1];}
                delta = S[li][mi] + sum;
                if (delta != galois::GaloisFieldElement (gf, 0))
                {
                    if (mi - ms[li] <= te - t_l[li])
                    {
                        //std::cout << "loook" << (xx^(mi - ms[li])) << "\n";
                        elp = elp - (delta/deltas[li]) * elps[li] * xx^(mi - ms[li]);
                    }
                    else
                    {
                        //std::cout << "loook" << (xx^(mi - ms[li])) << "\n";
                        int tt = te; galois::GaloisFieldPolynomial telp (elp);
                        elp = elp - (delta/deltas[li])*elps[li] * xx^(mi - ms[li]);
                        te = mi + t_l[li] - ms[li];
                        t_l[li] = tt; elps[li] = telp;
                        deltas[li] = delta; ms[li] = mi;

                    }
                }
            }
            //printf ("end\n");
        }
        
        //std::cout << xx << "\n";
    }

    std::cout << "t: " << te << " elp: " << elp << "\n";
    for (int i = 0; i < n; i ++)
    {
        if (elp ((*a)^(gf -> size()  - i)) == galois::GaloisFieldElement (gf, 0)) {std::cout << "e"; }
        else std::cout << "0"; 
    }



}