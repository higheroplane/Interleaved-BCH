
#ifndef BCH_GUARD
#define BCH_GUARD

#include "sourceCode/galois/GaloisField.h"
#include "sourceCode/galois/GaloisFieldElement.h"
#include "sourceCode/galois/GaloisFieldPolynomial.h"


typedef galois::GaloisField GF;

class BCH
{
    private:
    const int m, n, k, d;

    galois::GaloisField * gf;
    galois::GaloisFieldElement * a;
    //galois::GaloisFieldPolynomial * g;

    public:

    galois::GaloisFieldPolynomial * g;

    void dump ();

    BCH (int m, int n, int k);

    
    int decoder (int * error);

    int collaborative_decoder (int** error, int L);
    int collaborative_decoder_srl (int ** error, int l);

    galois::GaloisFieldPolynomial syndrome (int * error);

    int * gen_poly();
    


};

galois::GaloisFieldPolynomial min_poly (galois::GaloisFieldElement a);

#endif





