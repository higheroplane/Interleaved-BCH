
#ifndef LINALG_GUARD
#define LINALG_GUARD

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/tags.hpp>

#include "sourceCode/galois/GaloisField.h"
#include "sourceCode/galois/GaloisFieldElement.h"
#include "sourceCode/galois/GaloisFieldPolynomial.h"

typedef boost::numeric::ublas::matrix<galois::GaloisFieldElement> Mat;
typedef boost::numeric::ublas::vector<galois::GaloisFieldElement> Vec;

const int UPPER = 0;
const int LOWER = 1;

void LUdecomposition(Mat a, Mat& l, Mat& u, int n);
Mat mul (Mat A, Mat B);
Vec mul (Mat A, Vec b);
Vec solve_tr (Mat A, Vec b, int type);
Mat vconcat (Mat* Ms, int n);
Vec vconcat (Vec* Ms, int n);
bool is_null (Vec v);
galois::GaloisFieldElement det (Mat A);
bool is_full_rank (Mat A);
int polytriangular_submatrix (Mat& A, Vec& v, Mat*, Vec*, int*); //of full-rank martix

#endif
