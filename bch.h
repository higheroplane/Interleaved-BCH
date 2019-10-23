#include "sourceCode/galois/GaloisField.h"
#include "sourceCode/galois/GaloisFieldElement.h"
#include "sourceCode/galois/GaloisFieldPolynomial.h"

class BCH
{
    private:
    int m, n, k, d, t;

    galois::GaloisField * gf;
    galois::GaloisFieldElement * a;
    galois::GaloisFieldPolynomial * g;

    public:

    void dump ();

    BCH (int m, int n, int d);

    
    int decoder (int * error);

    int * gen_poly();
    


};

galois::GaloisFieldPolynomial min_poly (galois::GaloisFieldElement a);

