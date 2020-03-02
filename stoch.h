#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <ctime> 
#include <cstdlib>
#include <assert.h>
#include <vector>
#include "Bigint/bigint.h"
#include <boost/math/special_functions/binomial.hpp>

typedef long long unsigned int llu;


void make_hilbert_errors (std::vector<int*>& errv, int pgb, int pbg, int pe, int N);

std::vector<int> enumerative_decoder_w (BigInt k, int len, int weight);
std::vector<int> enumerative_decoder_mark (BigInt n, int len, int B, int C);
std::vector<int> fix_weight (int n, int weigth);

